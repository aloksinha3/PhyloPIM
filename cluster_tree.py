#!/usr/bin/env python3
"""
Tree Clustering by Depth Threshold (UPGMA)
===========================================
Cuts a rooted UPGMA tree at a depth threshold, producing one cluster per
subtree that crosses that depth — the direct algorithmic equivalent of
drawing vertical brackets on a printed phylogram.

Why UPGMA enables this cleanly
-------------------------------
UPGMA (Unweighted Pair Group Method with Arithmetic mean) produces a rooted,
ULTRAMETRIC tree: all leaves are equidistant from the root. This means depth
from the root is a meaningful, consistent measure across the whole tree.
Cutting at depth d is visually equivalent to drawing a vertical line at
position d on the phylogram — every branch crossing that line becomes a
cluster root.

Neighbour-Joining (NJ) trees are UNROOTED and not ultrametric, so depth
from an arbitrary root is not consistent. Clustering NJ trees required a
detour through a linkage matrix. With UPGMA, we cut the tree directly.

Auto-detection
--------------
If --k is not given, the script finds k by locating the largest gap between
consecutive merge depths. In an ultrametric tree this gap corresponds to the
longest inter-cluster branch — exactly what the eye sees when identifying
major groups on a phylogram.

Usage
-----
  python cluster_tree.py --newick V3_tree_upgma.nwk
  python cluster_tree.py --newick V3_tree_upgma.nwk --k 5
  python cluster_tree.py --newick V3_tree_upgma.nwk --k-min 3 --k-max 10
  python cluster_tree.py --newick V3_tree_upgma.nwk --out results/
"""

import argparse
import logging
import os
import sys
from io import StringIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

PALETTE = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#999999",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
]


# ─────────────────────────────────────────────────────────────────────────────
# Build UPGMA linkage from tree's patristic distances
# ─────────────────────────────────────────────────────────────────────────────

def build_upgma(newick_path: str) -> tuple[str, list, np.ndarray]:
    """
    Load Newick, compute patristic distances, build UPGMA linkage matrix.
    Returns (newick_str, names, Z).

    We always build Z from the tree's own patristic distances so the
    depth thresholds align exactly with the branch lengths in the PNG.
    """
    newick    = open(newick_path).read()
    tree      = Phylo.read(StringIO(newick), "newick")
    terminals = tree.get_terminals()
    names     = [t.name for t in terminals]
    n         = len(names)

    log.info("Building %d×%d patristic distance matrix …", n, n)
    dist = np.array([
        [tree.distance(terminals[i], terminals[j]) for j in range(n)]
        for i in range(n)
    ])
    dist = np.clip(dist, 0, None)

    # UPGMA = average linkage
    Z = linkage(squareform(dist), method="average")
    log.info("UPGMA linkage built — root depth = %.4f", Z[-1, 2])
    return newick, names, Z


# ─────────────────────────────────────────────────────────────────────────────
# Heavy / light chain lookup from Excel
# ─────────────────────────────────────────────────────────────────────────────

def load_heavy_light_from_xlsx(xlsx_path: str) -> dict:
    """
    Read (Sequence Number, Heavy Chain, Light Chain) from the active sheet of
    the original Excel file and return a mapping:

        {'Sequence_1': (heavy_seq, light_seq), ...}

    This lets us annotate cluster assignments with the exact input sequences.
    """
    try:
        import openpyxl
    except ImportError as e:
        raise RuntimeError(
            "openpyxl is required to read the Excel file for heavy/light chains.\n"
            "Install it with `pip install openpyxl` or re-run without --xlsx."
        ) from e

    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))

    mapping: dict[str, tuple[str, str]] = {}
    for row in rows[1:]:  # skip header
        if not row or row[0] is None:
            continue
        num = str(row[0]).strip()
        heavy = str(row[1]).strip() if len(row) > 1 and row[1] else ""
        light = str(row[2]).strip() if len(row) > 2 and row[2] else ""
        if not num:
            continue
        mapping[f"Sequence_{num}"] = (heavy, light)

    log.info("Loaded heavy/light sequences for %d IDs from %s", len(mapping), xlsx_path)
    return mapping


# ─────────────────────────────────────────────────────────────────────────────
# Depth threshold: the core operation
# ─────────────────────────────────────────────────────────────────────────────

def depth_for_k(Z: np.ndarray, k: int) -> float:
    """
    Return the depth threshold that yields exactly k clusters.

    In UPGMA, Z[n-k-1, 2] is the merge height that reduces the tree from
    k+1 clusters to k. Cutting at the midpoint between this height and the
    one above it (Z[n-k, 2]) lands cleanly in the gap between the two levels.

              Z[n-k, 2]   ← merge creating k+1 clusters (one level above)
    cut here ──────────── midpoint
              Z[n-k-1, 2] ← merge creating k clusters
    """
    n       = Z.shape[0] + 1
    h_at_k  = Z[n - k - 1, 2]   # depth at which we have k clusters
    h_above = Z[n - k,     2]   # depth one merge above (k+1 → k+2 clusters)
    return (h_at_k + h_above) / 2.0


def auto_k(Z: np.ndarray, k_min: int = 2, k_max: int = 20) -> tuple[int, float]:
    """
    Find k by locating the largest gap between consecutive UPGMA merge depths
    in the range [k_min, k_max].

    In an ultrametric tree the merge depths are strictly decreasing as k
    increases. The largest gap corresponds to the longest inter-cluster branch
    — the same feature the eye uses to identify major groups on a phylogram.

    Returns (k, gap_size) where k is the number of clusters just BEFORE
    the big gap (i.e. stay split at that level of detail).
    """
    n     = Z.shape[0] + 1
    k_max = min(k_max, n - 1)
    k_min = max(k_min, 2)

    ks      = list(range(k_min, k_max + 1))
    # Depths are the merge heights; higher k = smaller depth (more splits)
    depths  = np.array([Z[n - k - 1, 2] for k in ks])

    # gaps[i] = how much depth decreases going from ks[i] to ks[i+1]
    # A large gap = long branch was just severed = natural cluster boundary
    gaps     = np.diff(depths)          # all negative (depths decrease with k)
    best_idx = int(np.argmin(gaps))     # most negative = largest drop

    # Stay at ks[best_idx] — that's the k just BEFORE the big merge collapses it
    chosen_k = ks[best_idx]
    gap_size = float(abs(gaps[best_idx]))
    return chosen_k, gap_size


def print_depth_table(Z: np.ndarray, chosen_k: int,
                      k_min: int = 2, k_max: int = 20) -> None:
    """Print merge depths for the relevant k range."""
    n     = Z.shape[0] + 1
    k_max = min(k_max, n - 1)
    ks     = list(range(k_min, k_max + 1))
    depths = [Z[n - k - 1, 2] for k in ks]

    print(f"\n  {'k':>4}  {'merge depth':>12}  gap to k+1")
    print(f"  {'-'*4}  {'-'*12}  {'-'*35}")
    for i, (k, d) in enumerate(zip(ks, depths)):
        if i + 1 < len(depths):
            gap    = depths[i] - depths[i + 1]   # positive: depth decreases
            marker = "  ← largest gap, cut here" if k == chosen_k else ""
            print(f"  {k:>4}  {d:12.4f}  {gap:.4f}{marker}")
        else:
            print(f"  {k:>4}  {d:12.4f}")


# ─────────────────────────────────────────────────────────────────────────────
# Assign clusters and build output DataFrame
# ─────────────────────────────────────────────────────────────────────────────

def assign_clusters(names: list, Z: np.ndarray, depth: float) -> pd.DataFrame:
    """
    Cut the UPGMA tree at `depth` using scipy's distance criterion.
    fcluster with criterion='distance' returns every subtree whose root
    is at or above `depth` as a separate cluster — exactly the vertical-
    line-on-phylogram operation.
    Renumber clusters by descending size.
    """
    labels = fcluster(Z, t=depth, criterion="distance")
    df     = pd.DataFrame({"Seq_ID": names, "Cluster": labels.tolist()})

    order = df.groupby("Cluster").size().sort_values(ascending=False).index.tolist()
    remap = {old: new for new, old in enumerate(order, start=1)}
    df["Cluster"] = df["Cluster"].map(remap)
    return df.sort_values(["Cluster", "Seq_ID"]).reset_index(drop=True)


def print_summary(df: pd.DataFrame) -> None:
    print(f"\n  {'Cluster':>8}  {'Size':>6}   Members (first 5)")
    print(f"  {'-'*8}  {'-'*6}   {'-'*45}")
    for cid, grp in df.groupby("Cluster"):
        members = ", ".join(grp["Seq_ID"].tolist()[:5])
        suffix  = " …" if len(grp) > 5 else ""
        print(f"  {cid:>8}  {len(grp):>6}   {members}{suffix}")


# ─────────────────────────────────────────────────────────────────────────────
# Colour-coded PNG
# ─────────────────────────────────────────────────────────────────────────────

def render_png(newick: str, df: pd.DataFrame, depth: float,
               out_path: str, title: str = "") -> None:
    tree     = Phylo.read(StringIO(newick), "newick")
    n_leaves = tree.count_terminals()
    k        = df["Cluster"].nunique()

    colour_map  = {cid: PALETTE[i % len(PALETTE)]
                   for i, cid in enumerate(sorted(df["Cluster"].unique()))}
    leaf_colour = dict(zip(df["Seq_ID"], df["Cluster"].map(colour_map)))

    fig_h = max(10, n_leaves * 0.22)
    fig, ax = plt.subplots(figsize=(17, fig_h))

    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False,
        label_colors=leaf_colour,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )

    if not title:
        title = f"Phylogenetic Tree (UPGMA) — {k} clusters  (depth cut = {depth:.4f})"
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Evolutionary distance (Kimura-corrected substitutions/site)", fontsize=10)
    ax.set_yticks([])

    # Vertical dashed line at the depth cut — this IS the bracket Ellie drew
    ax.axvline(x=depth, color="black", linestyle="--", linewidth=0.9, alpha=0.6)

    handles = [
        mpatches.Patch(color=colour_map[cid],
                       label=f"Cluster {cid}  (n = {len(df[df.Cluster == cid])})")
        for cid in sorted(colour_map)
    ]
    handles.append(plt.Line2D(
        [0], [0], color="black", linestyle="--", linewidth=0.9,
        label=f"depth cut = {depth:.4f}",
    ))
    ax.legend(handles=handles, loc="lower right", fontsize=9,
              framealpha=0.85, edgecolor="grey", title=f"k = {k}")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("PNG → %s", out_path)


# ─────────────────────────────────────────────────────────────────────────────
# Orchestrator
# ─────────────────────────────────────────────────────────────────────────────

def run(newick_path, k=None, k_min=2, k_max=20, out_dir=".", png_path=None,
        xlsx_path=None):
    os.makedirs(out_dir, exist_ok=True)
    newick, names, Z = build_upgma(newick_path)

    heavy_light_map = None
    if xlsx_path:
        try:
            heavy_light_map = load_heavy_light_from_xlsx(xlsx_path)
        except Exception as exc:
            log.error("Failed to load heavy/light sequences from %s: %s",
                      xlsx_path, exc)
            heavy_light_map = None

    if k is not None:
        chosen_k = k
        log.info("Using user-specified k = %d", chosen_k)
    else:
        chosen_k, gap = auto_k(Z, k_min=k_min, k_max=k_max)
        log.info("Auto-detected k = %d  (gap = %.4f)", chosen_k, gap)
        print_depth_table(Z, chosen_k, k_min, k_max)

    depth = depth_for_k(Z, chosen_k)
    log.info("Depth cut = %.4f  (midpoint between k=%d and k=%d merge depths)",
             depth, chosen_k, chosen_k + 1)

    df    = assign_clusters(names, Z, depth)

    # If we have the original Excel, add heavy/light chain sequences
    if heavy_light_map is not None:
        heavy_vals = []
        light_vals = []
        for sid in df["Seq_ID"]:
            h, l = heavy_light_map.get(sid, ("", ""))
            heavy_vals.append(h)
            light_vals.append(l)
        df["Heavy_chain"] = heavy_vals
        df["Light_chain"] = light_vals
        # Ensure column order: Seq_ID, Cluster, Heavy_chain, Light_chain, ...
        base_cols = ["Seq_ID", "Cluster", "Heavy_chain", "Light_chain"]
        other_cols = [c for c in df.columns if c not in base_cols]
        df = df[base_cols + other_cols]
    k_out = df["Cluster"].nunique()
    log.info("Produced %d clusters", k_out)
    print_summary(df)

    csv_path = os.path.join(out_dir, "cluster_assignments.csv")
    df.to_csv(csv_path, index=False)
    log.info("CSV → %s", csv_path)

    if png_path is None:
        png_path = os.path.join(out_dir, "tree_clustered.png")
    render_png(newick, df, depth, png_path)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--newick",  required=True,
                   help="Input Newick file from UPGMA tree (e.g. from ebi_tree.py)")
    p.add_argument("--xlsx", default=None,
                   help="Original Excel file with Sequence Number | Heavy | Light "
                        "to annotate clusters with heavy/light chains")
    p.add_argument("--k",       type=int, default=None,
                   help="Number of clusters (skips auto-detection)")
    p.add_argument("--k-min",   type=int, default=2,
                   help="Min k for auto-detection  [default: 2]")
    p.add_argument("--k-max",   type=int, default=20,
                   help="Max k for auto-detection  [default: 20]")
    p.add_argument("--out",     default=".",
                   help="Output directory  [default: .]")
    p.add_argument("--png",     default=None,
                   help="Output PNG  [default: <out>/tree_clustered.png]")
    args = p.parse_args()
    run(newick_path=args.newick, k=args.k,
        k_min=args.k_min, k_max=args.k_max,
        out_dir=args.out, png_path=args.png,
        xlsx_path=args.xlsx)
