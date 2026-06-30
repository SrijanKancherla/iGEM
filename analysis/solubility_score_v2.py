"""Solubility scoring for the v2 construct-based pipeline."""

import sys
from pathlib import Path
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
# Reuse the scoring function from v1 — only paths change
sys.path.insert(0, str(Path(__file__).parent.parent / "analysis"))
from solubility_score import add_solubility_to_pipeline

INPUT_CSV   = Path("construct_v2/results/6_expression.csv")
OUTPUT_CSV  = Path("construct_v2/results/7_solubility.csv")
WT_FASTA_V2 = Path("construct_v2/data/cleaned/hbsag_construct.fasta")


def main():
    wt_seq = ""
    for rec in SeqIO.parse(WT_FASTA_V2, "fasta"):
        wt_seq = str(rec.seq).upper()
        break
    if not wt_seq:
        raise ValueError(f"Could not read WT sequence from {WT_FASTA_V2}")

    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
    df = add_solubility_to_pipeline(
        input_csv=str(INPUT_CSV),
        output_csv=str(OUTPUT_CSV),
        wt_sequence=wt_seq,
        verbose=False,
    )
    print(f"Solubility scoring complete → {OUTPUT_CSV}")
    print(df["solubility_score"].describe())


if __name__ == "__main__":
    main()
