import pandas as pd
from pathlib import Path

EXCEL_PATH   = Path("V3_mAb_sequences.xlsx")
ID_COLUMN    = "Ab ID"
SEQIDX_COL   = "signed #"
VH_COL       = "VH AA Sequence"
VL_COL       = "VK/VL AA Sequence"
FASTA_OUT    = Path("V3_mAb_fab.fasta")
IDMAP_OUT    = Path("V3_mAb_fab_id_map.csv")


def main() -> None:
    if not EXCEL_PATH.exists():
        raise SystemExit(f"Excel file not found: {EXCEL_PATH.resolve()}")

    df = pd.read_excel(EXCEL_PATH)

    missing = [c for c in (ID_COLUMN, SEQIDX_COL, VH_COL, VL_COL) if c not in df.columns]
    if missing:
        raise SystemExit(
            f"Missing expected column(s) {missing} in {EXCEL_PATH.name}.\n"
            f"Found columns: {list(df.columns)}"
        )

    df = df.dropna(subset=[ID_COLUMN, SEQIDX_COL, VH_COL, VL_COL])

    records: list[tuple[str, str]] = []
    id_map_rows = []
    # Enumerate rows after filtering so X = 1..N in final FASTA
    for idx, (_, row) in enumerate(df.iterrows(), start=1):
        # Use Sequence_X where X is the sequential index 1..N
        seq_id = f"Sequence_{idx}"

        vh = str(row[VH_COL]).strip().replace(" ", "").replace("\n", "")
        vl = str(row[VL_COL]).strip().replace(" ", "").replace("\n", "")
        if not vh or not vl:
            continue
        fab_seq = vh + vl
        records.append((seq_id, fab_seq))
        id_map_rows.append(
            {
                "Sequence_ID": seq_id,
                "Sequence_Number": idx,
                "Original_Row_Index": row[SEQIDX_COL],
                "Ab_ID": row[ID_COLUMN],
            }
        )

    if not records:
        raise SystemExit("No valid Fab sequences (VH + VL) found in the Excel file.")

    with FASTA_OUT.open("w") as fh:
        for seq_id, seq in records:
            fh.write(f">{seq_id}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # Write mapping table so you can link Sequence_X back to spreadsheet IDs.
    pd.DataFrame(id_map_rows).to_csv(IDMAP_OUT, index=False)

    print(f"Wrote {len(records)} Fab sequences to {FASTA_OUT.resolve()}")
    print(f"Wrote ID mapping table to {IDMAP_OUT.resolve()}")


if __name__ == "__main__":
    main()

