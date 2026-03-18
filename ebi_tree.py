#!/usr/bin/env python3
"""
V3 mAb Phylogenetic Tree Pipeline
===================================
Reads heavy+light chain sequences from an Excel file, concatenates them
into full Fab sequences, submits them to EBI Clustal Omega for alignment,
sends the alignment to EBI Simple Phylogeny with Distance Correction ON,
and saves a publication-quality PNG of the resulting tree.

EBI REST base URLs:
  https://www.ebi.ac.uk/Tools/services/rest/clustalo
  https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny

No API key required — just supply a valid email address.

Usage
-----
  python ebi_tree.py --xlsx V3_mAb_sequences.xlsx --email you@example.com
  python ebi_tree.py --xlsx V3_mAb_sequences.xlsx --email you@example.com --out tree.png
  python ebi_tree.py --list-params --email you@example.com   # diagnose accepted params
  python ebi_tree.py --help
"""

import argparse
import logging
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests
import openpyxl
from io import StringIO
from Bio import Phylo

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# Suppress SSL verification warnings that appear when verify=False is used.
# This is necessary on institutional networks with corporate SSL proxies.
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

CLUSTALO_URL     = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
PHYLOGENY_URL    = "https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny"
POLL_INTERVAL    = 5     # seconds between status polls
MAX_WAIT_SECONDS = 900   # 15 min hard ceiling per job


# ─────────────────────────────────────────────────────────────────────────────
# EBI REST helpers
# ─────────────────────────────────────────────────────────────────────────────

def _submit(base_url: str, params: dict) -> str:
    """POST to /run; return the plain-text job ID or raise with full response body."""
    resp = requests.post(f"{base_url}/run", data=params, timeout=30, verify=False)
    if not resp.ok:
        raise RuntimeError(
            f"HTTP {resp.status_code} from {base_url}/run\n"
            f"Response body: {resp.text[:800]}"
        )
    job_id = resp.text.strip()
    if not job_id or "<" in job_id:
        raise RuntimeError(
            f"Expected a plain job ID but got HTML from {base_url}/run.\n"
            f"First 300 chars: {resp.text[:300]}"
        )
    log.info("Submitted → job ID: %s", job_id)
    return job_id


def _wait(base_url: str, job_id: str) -> None:
    """Poll /status/<jobId> until FINISHED (or raise on error/timeout)."""
    url      = f"{base_url}/status/{job_id}"
    deadline = time.time() + MAX_WAIT_SECONDS
    while time.time() < deadline:
        resp   = requests.get(url, timeout=15, verify=False)
        resp.raise_for_status()
        status = resp.text.strip()
        log.info("  status: %s", status)
        if status == "FINISHED":
            return
        if status in ("ERROR", "FAILURE", "NOT_FOUND"):
            raise RuntimeError(f"Job {job_id} ended with status: {status}")
        time.sleep(POLL_INTERVAL)
    raise TimeoutError(f"Job {job_id} did not finish within {MAX_WAIT_SECONDS}s")


def _result(base_url: str, job_id: str, result_type: str) -> str:
    """Fetch a named result for a finished job."""
    resp = requests.get(f"{base_url}/result/{job_id}/{result_type}", timeout=60, verify=False)
    resp.raise_for_status()
    return resp.text


# ─────────────────────────────────────────────────────────────────────────────
# Diagnostics
# ─────────────────────────────────────────────────────────────────────────────

def list_phylogeny_params() -> None:
    """
    Fetch /parameters and /parameterdetails/<id> for Simple Phylogeny,
    parse the XML, and print a clean table of every parameter with its
    allowed values and default. Useful for diagnosing 400 errors.
    """
    import xml.etree.ElementTree as ET

    resp = requests.get(f"{PHYLOGENY_URL}/parameters", timeout=20, verify=False)
    resp.raise_for_status()
    root = ET.fromstring(resp.text)
    param_ids = [el.text for el in root.findall("id") if el.text]

    print(f"\nSimple Phylogeny — {len(param_ids)} parameters:\n")
    print(f"  {'Parameter':<16} {'Default':<14} Allowed values")
    print(f"  {'-'*16} {'-'*14} {'-'*44}")

    for pid in param_ids:
        d_resp = requests.get(f"{PHYLOGENY_URL}/parameterdetails/{pid}", timeout=20, verify=False)
        d_resp.raise_for_status()
        d = ET.fromstring(d_resp.text)

        default = d.findtext("defaultValue") or ""
        values  = [v.findtext("value") or "" for v in d.findall(".//values/value")]

        if values:
            print(f"  {pid:<16} {default:<14} {', '.join(values)}")
        else:
            ptype = d.findtext("type") or "string"
            print(f"  {pid:<16} {default:<14} <{ptype}>")


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 – Read sequences from Excel
# ─────────────────────────────────────────────────────────────────────────────

def load_sequences_from_xlsx(path: str) -> dict:
    """
    Read (Sequence Number, Heavy Chain, Light Chain) from the active sheet.
    Returns {seq_label: heavy+light concatenated string}.
    """
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))

    sequences = {}
    for row in rows[1:]:   # skip header
        if row[0] is None:
            continue
        num   = row[0]
        heavy = str(row[1]).strip() if row[1] else ""
        light = str(row[2]).strip() if row[2] else ""
        if heavy:
            sequences[f"Sequence_{num}"] = heavy + light

    log.info("Loaded %d sequences from %s", len(sequences), path)
    return sequences


def _to_fasta(sequences: dict) -> str:
    return "".join(f">{sid}\n{seq}\n" for sid, seq in sequences.items())


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 – EBI Clustal Omega alignment
# ─────────────────────────────────────────────────────────────────────────────

def run_clustalo(sequences: dict, email: str) -> str:
    """
    Submit Fab sequences to EBI Clustal Omega.
    Returns the alignment in Clustal format (required by Simple Phylogeny).
    """
    log.info("=== Submitting to EBI Clustal Omega (%d sequences) ===", len(sequences))
    job_id = _submit(CLUSTALO_URL, {
        "email":    email,
        "sequence": _to_fasta(sequences),
        "stype":    "protein",
        "outfmt":   "clustal",
    })
    _wait(CLUSTALO_URL, job_id)
    alignment = _result(CLUSTALO_URL, job_id, "aln-clustal")
    log.info("Alignment received (%d chars)", len(alignment))
    return alignment


# ─────────────────────────────────────────────────────────────────────────────
# Step 3 – EBI Simple Phylogeny (Distance Correction ON)
# ─────────────────────────────────────────────────────────────────────────────

def run_simple_phylogeny(alignment: str, email: str,
                         clustering: str = "UPGMA") -> str:
    """
    Submit a Clustal-format alignment to EBI Simple Phylogeny.

    Parameters used:
      kimura     = "true"            →  Distance Correction ON (Kimura 2-parameter model)
      tossgaps   = "true"            →  Exclude gap-containing columns
      clustering = "UPGMA"           →  UPGMA produces a rooted, ultrametric tree,
                                        which allows direct depth-threshold clustering.
                                        Use "Neighbour-joining" for the classic unrooted NJ tree.
      tree       = "phylip"          →  Standard Newick output format

    Confirmed accepted values (run --list-params to verify):
      clustering : "UPGMA", "Neighbour-joining"
      tree       : "phylip", "nj", "dist", "nexus"
      kimura     : "true", "false"
      tossgaps   : "true", "false"
    """
    log.info("=== Submitting to EBI Simple Phylogeny (clustering=%s) ===", clustering)
    params = {
        "email":      email,
        "sequence":   alignment,
        "stype":      "protein",
        "kimura":     "true",      # Distance Correction ON
        "tossgaps":   "true",      # Exclude gapped columns
        "clustering": clustering,  # "UPGMA" or "Neighbour-joining"
        "tree":       "phylip",    # Newick output format
    }
    job_id = _submit(PHYLOGENY_URL, params)
    _wait(PHYLOGENY_URL, job_id)
    newick = _result(PHYLOGENY_URL, job_id, "tree")
    log.info("Newick tree received (%d chars)", len(newick))
    return newick


# ─────────────────────────────────────────────────────────────────────────────
# Step 4 – Render tree as PNG
# ─────────────────────────────────────────────────────────────────────────────

def render_tree_png(newick: str, out_path: str, title: str = "V3 mAb Phylogenetic Tree (UPGMA)") -> None:
    """
    Parse the Newick string with Biopython and render a horizontal
    phylogram as a high-resolution PNG.
    """
    tree = Phylo.read(StringIO(newick), "newick")
    n_leaves = tree.count_terminals()

    fig_h = max(8, n_leaves * 0.22)
    fig, ax = plt.subplots(figsize=(16, fig_h))

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )

    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("Evolutionary distance (Kimura-corrected substitutions/site)", fontsize=10)
    ax.set_ylabel("")
    ax.set_yticks([])

    note = mpatches.Patch(
        color="none",
        label="Distance correction: ON  |  Method: Neighbour-Joining  |  Model: Kimura"
    )
    ax.legend(handles=[note], loc="lower right", fontsize=8,
              framealpha=0.6, edgecolor="grey")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Tree PNG saved → %s", out_path)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--xlsx",  default=None,
                   help="Excel file: Sequence Number | Heavy Chain | Light Chain")
    p.add_argument("--email", required=True,
                   help="Your email address (required by EBI REST API)")
    p.add_argument("--out",   default="phylogenetic_tree.png",
                   help="Output PNG path  [default: phylogenetic_tree.png]")
    p.add_argument("--save-intermediates", action="store_true",
                   help="Also save the alignment (.aln) and Newick tree (.nwk)")
    p.add_argument("--clustering", default="UPGMA",
                   choices=["UPGMA", "Neighbour-joining"],
                   help="Tree-building method  [default: UPGMA]")
    p.add_argument("--list-params", action="store_true",
                   help="Print all parameters accepted by Simple Phylogeny API and exit")
    args = p.parse_args()

    if args.list_params:
        list_phylogeny_params()
        return

    if not args.xlsx:
        p.error("--xlsx is required unless using --list-params")

    sequences = load_sequences_from_xlsx(args.xlsx)
    if len(sequences) < 2:
        sys.exit("Need at least 2 sequences to build a tree.")

    alignment = run_clustalo(sequences, args.email)
    newick    = run_simple_phylogeny(alignment, args.email,
                                     clustering=args.clustering)

    if args.save_intermediates:
        base = args.out.replace(".png", "")
        with open(f"{base}.aln", "w") as fh:
            fh.write(alignment)
        with open(f"{base}.nwk", "w") as fh:
            fh.write(newick)
        log.info("Alignment → %s.aln", base)
        log.info("Newick    → %s.nwk", base)

    render_tree_png(newick, args.out)
    print(f"\nDone! Tree image: {args.out}")


if __name__ == "__main__":
    main()
