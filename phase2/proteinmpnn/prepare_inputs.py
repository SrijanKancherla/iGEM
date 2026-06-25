import argparse
import json
from pathlib import Path

import pandas as pd


DEFAULT_INPUT = "results/final_ranked.csv"
DEFAULT_OUTPUT_DIR = Path("phase2/proteinmpnn/inputs")
DEFAULT_PDB = "data/cleaned/hbsag.pdb"


def choose_candidates(df, top_n):
    usable = df[df["recommendation"].isin(["strong_candidate", "candidate", "review_manually"])].copy()
    usable = usable[~usable["flags"].fillna("").str.contains("possible_disulfide_loss|buried_charge|low_solubility")]
    if usable.empty:
        usable = df.copy()
    return usable.sort_values(["final_score", "confidence_score"], ascending=False).head(top_n)


def main():
    parser = argparse.ArgumentParser(description="Prepare Phase 2 ProteinMPNN inputs.")
    parser.add_argument("--ranked", default=DEFAULT_INPUT)
    parser.add_argument("--pdb", default=DEFAULT_PDB)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.ranked)
    candidates = choose_candidates(df, args.top_n)
    candidates.to_csv(out_dir / "proteinmpnn_candidates.csv", index=False)

    mutable_positions = sorted(int(pos) for pos in candidates["pdb_residue"].unique())
    fixed_positions = {
        Path(args.pdb).name: {
            "A": mutable_positions,
            "note": "These are candidate mutable PDB residue numbers. Freeze all other positions in ProteinMPNN.",
        }
    }

    with open(out_dir / "mutable_positions.jsonl", "w") as handle:
        handle.write(json.dumps(fixed_positions) + "\n")

    with open(out_dir / "run_notes.txt", "w") as handle:
        handle.write(
            "ProteinMPNN preparation complete.\n"
            f"PDB: {args.pdb}\n"
            "Candidate CSV: proteinmpnn_candidates.csv\n"
            "Mutable positions JSONL: mutable_positions.jsonl\n\n"
            "Recommended use: constrain/freeze immune epitope residues, antibody-contact residues,\n"
            "critical cysteines, and highly conserved positions unless they are intentionally being tested.\n"
        )

    print(f"ProteinMPNN inputs written to {out_dir}")
    print(f"Candidate count: {len(candidates)}")


if __name__ == "__main__":
    main()
