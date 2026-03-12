#!/usr/bin/env python3
"""
CDR PIM Analysis — Steps 3 & 4 of the HVTN300 representative selection pipeline
=================================================================================
Takes cluster_assignments_CDR.csv (output of scalop.py, which contains a
CDR_concat column with all six CDR loops concatenated) and for each cluster:

  Step 3-1  Align the CDR_concat sequences within each cluster using
            Clustal Omega (local binary).

  Step 3-2  Compute the Percentage Identity Matrix (PIM) from the alignment.

  Step 3-3  Rank sequences by average % identity to all others in the cluster
            (higher = more consensus-like, lower = more divergent).

  Step 4    Pick representative sequences two ways per cluster:
              top_bottom : rank #1 (highest) + rank #last (lowest)
              trisect    : --n-reps evenly-spaced positions spanning the full
                           range (default n=4 → captures high, 2 midpoints,
                           low; n=2 → same as top_bottom)

Outputs (all written to --out directory)
-----------------------------------------
  cluster_N_pim.csv               full PIM matrix for cluster N
  cluster_N_ranked.csv            ranked sequences with Avg_%_Identity + length
  representatives_top_bottom.csv  one file: highest + lowest per cluster
  representatives_trisect.csv     one file: n trisect-point reps per cluster
  PIM_summary.csv                 all clusters combined, with rep flags

Usage
-----
  python PIM.py --csv cluster_assignments_CDR.csv
  python PIM.py --csv cluster_assignments_CDR.csv --n-reps 4 --out results/
  python PIM.py --csv cluster_assignments_CDR.csv --clustalo /usr/bin/clustalo
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from io import StringIO

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Clustal Omega helpers
# ─────────────────────────────────────────────────────────────────────────────

def _find_clustalo(hint: str = None) -> str:
    if hint and os.path.isfile(hint):
        return hint
    path = shutil.which("clustalo") or shutil.which("clustal-omega")
    if path:
        return path
    raise FileNotFoundError(
        "clustalo not found on PATH. Install with:\n"
        "  conda install -c bioconda clustalo   OR\n"
        "  brew install clustal-omega           OR\n"
        "  sudo apt install clustalo\n"
        "or pass --clustalo /path/to/clustalo"
    )


def _to_fasta(seqs: dict) -> str:
    """Dict {name: sequence} → FASTA string."""
    return "".join(f">{name}\n{seq}\n" for name, seq in seqs.items())


def run_clustalo_pim(seqs: dict, clustalo_bin: str) -> pd.DataFrame:
    """
    Align sequences with Clustal Omega and return the PIM as a DataFrame.
    Uses --percent-id --distmat-out to get % identity directly.
    For clusters of size 1 returns a trivial 100% matrix.
    """
    names = list(seqs.keys())

    if len(names) == 1:
        return pd.DataFrame([[100.0]], index=names, columns=names)

    with tempfile.TemporaryDirectory() as tmpdir:
        fa_in   = os.path.join(tmpdir, "in.fa")
        aln_out = os.path.join(tmpdir, "aln.fa")
        pim_out = os.path.join(tmpdir, "pim.txt")

        with open(fa_in, "w") as fh:
            fh.write(_to_fasta(seqs))

        cmd = [
            clustalo_bin,
            "-i", fa_in,
            "-o", aln_out,
            "--distmat-out", pim_out,
            "--percent-id",   # output % identity instead of distance
            "--full",         # full matrix (not mBed approximation)
            "--force",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"clustalo failed:\n{result.stderr}"
            )

        # Parse the PIM: first line = n, then n lines of name + n floats
        with open(pim_out) as fh:
            lines = [l.strip() for l in fh if l.strip()]

        n      = int(lines[0])
        matrix = np.zeros((n, n))
        row_names = []
        for i, line in enumerate(lines[1: n + 1]):
            parts     = line.split()
            row_names.append(parts[0])
            matrix[i] = [float(x) for x in parts[1:]]

    return pd.DataFrame(matrix, index=row_names, columns=row_names)


# ─────────────────────────────────────────────────────────────────────────────
# PIM → ranking
# ─────────────────────────────────────────────────────────────────────────────

def rank_by_avg_pim(pim: pd.DataFrame) -> pd.DataFrame:
    """
    For each sequence compute average % identity to ALL others (excluding self).
    Returns a DataFrame sorted descending: Rank, Seq_ID, Avg_%_Identity.
    """
    n    = len(pim)
    avgs = {}
    for name in pim.index:
        others     = pim.loc[name].drop(name)
        avgs[name] = others.mean() if len(others) > 0 else 100.0

    ranked = (
        pd.DataFrame.from_dict(avgs, orient="index", columns=["Avg_%_Identity"])
        .reset_index()
        .rename(columns={"index": "Seq_ID"})
        .sort_values("Avg_%_Identity", ascending=False)
        .reset_index(drop=True)
    )
    ranked.insert(0, "Rank", ranked.index + 1)
    return ranked


# ─────────────────────────────────────────────────────────────────────────────
# Representative selection
# ─────────────────────────────────────────────────────────────────────────────

def pick_top_bottom(ranked: pd.DataFrame) -> pd.DataFrame:
    """
    Select rank #1 (highest avg % identity) and rank #last (lowest).
    Returns a copy of those rows with a 'Selection' label column.
    """
    top    = ranked.iloc[[0]].copy()
    bottom = ranked.iloc[[-1]].copy()
    top["Selection"]    = "highest"
    bottom["Selection"] = "lowest"
    return pd.concat([top, bottom], ignore_index=True)


def pick_trisect(ranked: pd.DataFrame, n: int = 4) -> pd.DataFrame:
    """
    Select n sequences at evenly-spaced positions spanning the full ranked list.
    Mirrors Ellie's trisect approach:
      n=2 → positions [0, -1]           (same as top_bottom)
      n=4 → positions [0, 1/3, 2/3, -1] of the range

    For a list of length L, positions are computed as:
      np.linspace(0, L-1, n) rounded to nearest integer.
    """
    L       = len(ranked)
    n       = min(n, L)
    indices = np.round(np.linspace(0, L - 1, n)).astype(int)
    # Deduplicate while preserving order (can happen for tiny clusters)
    seen, unique_idx = set(), []
    for i in indices:
        if i not in seen:
            unique_idx.append(i)
            seen.add(i)

    selected = ranked.iloc[unique_idx].copy()
    selected["Selection"] = [f"trisect_{pos+1}_of_{len(unique_idx)}"
                              for pos in range(len(unique_idx))]
    return selected.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# Per-cluster pipeline
# ─────────────────────────────────────────────────────────────────────────────

def process_cluster(cluster_id: int, group: pd.DataFrame,
                    clustalo_bin: str, out_dir: str,
                    n_reps: int = 4) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run the full Step 3–4 pipeline for one cluster.
    Returns (ranked_df, top_bottom_df, trisect_df).
    """
    log.info("Cluster %d: %d sequences", cluster_id, len(group))

    # Build {Seq_ID: CDR_concat} dict
    seqs = dict(zip(group["Seq_ID"], group["CDR_concat"].astype(str)))

    # Step 3-1 & 3-2: align + PIM
    pim = run_clustalo_pim(seqs, clustalo_bin)

    # Save PIM
    pim_path = os.path.join(out_dir, f"cluster_{cluster_id}_pim.csv")
    pim.to_csv(pim_path)
    log.info("  PIM saved → %s", pim_path)

    # Step 3-3: rank
    ranked = rank_by_avg_pim(pim)

    # Attach CDR length for reference
    len_map = group.set_index("Seq_ID")["CDR_concat"].str.len()
    ranked["CDR_length"] = ranked["Seq_ID"].map(len_map)
    ranked.insert(1, "Cluster", cluster_id)

    ranked_path = os.path.join(out_dir, f"cluster_{cluster_id}_ranked.csv")
    ranked.to_csv(ranked_path, index=False)
    log.info("  Ranked table saved → %s", ranked_path)

    # Step 4: pick representatives
    tb       = pick_top_bottom(ranked)
    trisect  = pick_trisect(ranked, n=n_reps)

    tb["Cluster"]      = cluster_id
    trisect["Cluster"] = cluster_id

    return ranked, tb, trisect


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def run(csv_path: str, clustalo_bin: str = None, out_dir: str = ".",
        n_reps: int = 4) -> None:

    os.makedirs(out_dir, exist_ok=True)
    clustalo_bin = _find_clustalo(clustalo_bin)
    log.info("Using clustalo: %s", clustalo_bin)

    df = pd.read_csv(csv_path)
    required = {"Seq_ID", "Cluster", "CDR_concat"}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {missing}\n"
            f"Found: {list(df.columns)}\n"
            "Make sure you pass the output of scalop.py (cluster_assignments_CDR.csv)."
        )

    log.info("Loaded %d sequences across %d clusters",
             len(df), df["Cluster"].nunique())

    all_ranked    = []
    all_tb        = []
    all_trisect   = []

    for cluster_id, group in df.groupby("Cluster"):
        ranked, tb, trisect = process_cluster(
            cluster_id, group, clustalo_bin, out_dir, n_reps
        )
        all_ranked.append(ranked)
        all_tb.append(tb)
        all_trisect.append(trisect)

    # ── Assemble summary ──────────────────────────────────────────────────────
    summary = pd.concat(all_ranked, ignore_index=True)

    tb_ids      = set(pd.concat(all_tb)["Seq_ID"])
    trisect_ids = set(pd.concat(all_trisect)["Seq_ID"])
    summary["rep_top_bottom"] = summary["Seq_ID"].isin(tb_ids)
    summary["rep_trisect"]    = summary["Seq_ID"].isin(trisect_ids)

    summary_path = os.path.join(out_dir, "PIM_summary.csv")
    summary.to_csv(summary_path, index=False)
    log.info("Summary → %s", summary_path)

    # ── Representatives ───────────────────────────────────────────────────────
    tb_df      = pd.concat(all_tb,      ignore_index=True)
    trisect_df = pd.concat(all_trisect, ignore_index=True)

    # Column order: Cluster | Rank | Seq_ID | Avg_%_Identity | CDR_length | Selection
    col_order = ["Cluster", "Rank", "Seq_ID", "Avg_%_Identity", "CDR_length", "Selection"]

    tb_path = os.path.join(out_dir, "representatives_top_bottom.csv")
    tb_df[col_order].to_csv(tb_path, index=False)
    log.info("Top/bottom reps → %s", tb_path)

    trisect_path = os.path.join(out_dir, "representatives_trisect.csv")
    trisect_df[col_order].to_csv(trisect_path, index=False)
    log.info("Trisect reps    → %s", trisect_path)

    # ── Print final tables ────────────────────────────────────────────────────
    print("\n── Top/Bottom representatives (n=2 per cluster) ──")
    print(tb_df[col_order].to_string(index=False))

    print(f"\n── Trisect representatives (n={n_reps} per cluster) ──")
    print(trisect_df[col_order].to_string(index=False))


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--csv",      required=True,
                   help="cluster_assignments_CDR.csv from scalop.py")
    p.add_argument("--out",      default=".",
                   help="Output directory  [default: current directory]")
    p.add_argument("--n-reps",   type=int, default=4,
                   help="Trisect representatives per cluster  [default: 4]")
    p.add_argument("--clustalo", default=None,
                   help="Path to clustalo binary  [default: auto-detect from PATH]")
    args = p.parse_args()

    run(
        csv_path     = args.csv,
        clustalo_bin = args.clustalo,
        out_dir      = args.out,
        n_reps       = args.n_reps,
    )
