# PhyloPIM

Phylogenetic clustering and CDR representative selection pipeline for antibody sequences.

Given a set of antibody heavy and light chain sequences, PhyloPIM:
1. Builds a UPGMA phylogenetic tree via EBI Clustal Omega + Simple Phylogeny
2. Cuts the tree into k clusters using a depth threshold
3. Extracts CDR loops (CDRH1–CDRL3) for each sequence using SCALOP
4. Computes a Percentage Identity Matrix (PIM) per cluster from the concatenated CDR sequences
5. Ranks sequences by average CDR identity and selects representative sequences

## Installation

```bash
cd Desktop
git clone https://github.com/aloksinha3/PhyloPIM.git
cd PhyloPIM
```

Create and activate the conda environment:

```bash
conda create -n PhyloPIM python=3.10 -y
conda activate PhyloPIM
conda install -c bioconda -c conda-forge --file requirements.txt -y
git clone https://github.com/oxpig/SCALOP
pip install ./SCALOP/
```

## Input format

All input Excel files must follow the format of `V3_mAb_sequences.xlsx`:

| Column 1          | Column 2             | Column 3              |
|-------------------|----------------------|-----------------------|
| Sequence Number   | Heavy Chain sequence | Light Chain sequence  |

## Usage

### Full pipeline

```bash
python main.py \
  --xlsx V3_mAb_sequences.xlsx \
  --email your_email@email_provider.com \
  --k 5 \
  --out results/
```

`--k` sets the number of clusters. If omitted, the pipeline auto-detects k from the largest gap in the UPGMA merge depths, which may not always match the visual groupings in the tree

### Resume from a later step

If the EBI tree has already been built, skip steps 1–2:

```bash
python main.py \
  --xlsx V3_mAb_sequences.xlsx \
  --email your_email@email_provider.com \
  --start-from cluster \
  --newick results/tree.nwk \
  --k 5 \
  --out results/
```

Resume from CDR extraction (clusters already assigned):

```bash
python main.py \
  --start-from scalop \
  --cluster-csv results/cluster_assignments.csv \
  --out results/
```

Resume from PIM only (CDR CSV already exists):

```bash
python main.py \
  --start-from pim \
  --cdr-csv results/cluster_assignments_CDR.csv \
  --out results/
```

### All options

| Flag | Default | Description |
|------|---------|-------------|
| `--xlsx` | — | Input Excel file |
| `--email` | — | Email address for EBI REST API (required for steps 1–2) |
| `--k` | auto | Number of clusters |
| `--k-min` | 2 | Minimum k for auto-detection |
| `--k-max` | 20 | Maximum k for auto-detection |
| `--n-reps` | 4 | Trisect representative sequences per cluster |
| `--clustering-method` | UPGMA | Tree method: `UPGMA` or `Neighbour-joining` |
| `--ncpu` | 1 | Parallel workers for SCALOP |
| `--out` | `results/` | Output directory |
| `--start-from` | `fasta` | Resume from: `fasta`, `tree`, `cluster`, `scalop`, `pim` |
| `--newick` | — | Existing `.nwk` file (required when `--start-from cluster`) |
| `--cluster-csv` | — | Existing `cluster_assignments.csv` (required when `--start-from scalop`) |
| `--cdr-csv` | — | Existing `cluster_assignments_CDR.csv` (required when `--start-from pim`) |

## Outputs

All files are written to `--out` (default: `results/`):

| File | Description |
|------|-------------|
| `tree.nwk` | Newick phylogenetic tree |
| `tree.png` | Phylogenetic tree image |
| `tree_clustered.png` | Colour-coded tree with cluster assignments |
| `cluster_assignments.csv` | Sequence ID, cluster number, heavy/light chains |
| `cluster_assignments_CDR.csv` | Above + CDRH1–CDRL3 sequences, positions, and concatenated CDR string |
| `cluster_N_pim.csv` | Full N×N percentage identity matrix for cluster N |
| `cluster_N_ranked.csv` | Sequences in cluster N ranked by average CDR % identity |
| `PIM_summary.csv` | All clusters combined with rank, avg % identity, and representative flags |
| `representatives_top_bottom.csv` | Highest + lowest ranked sequence per cluster (2 per cluster) |
| `representatives_trisect.csv` | `--n-reps` evenly-spaced sequences spanning the identity range per cluster |

## Pipeline scripts

| Script | Role |
|--------|------|
| `ebi_tree.py` | EBI Clustal Omega + Simple Phylogeny → Newick tree |
| `cluster_tree.py` | UPGMA depth-cut clustering → `cluster_assignments.csv` |
| `extract_cdrs.py` | SCALOP CDR extraction → `cluster_assignments_CDR.csv` |
| `PIM.py` | CDR alignment, PIM, ranking, representative selection |
| `main.py` | Orchestrates all steps end-to-end |

Each script can also be run independently — see the docstring at the top of each file for its own usage instructions.
