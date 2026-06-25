import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO


DEFAULT_INPUT = "results/final_ranked.csv"
DEFAULT_FASTA = "data/cleaned/hbsag.fasta"
DEFAULT_OUTPUT_DIR = Path("phase2/md/inputs")


def read_sequence(path):
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq).upper()
    raise ValueError(f"No FASTA sequence found in {path}")


def mutate_sequence(sequence, pos, mut):
    chars = list(sequence)
    chars[int(pos)] = mut
    return "".join(chars)


def choose_candidates(df, top_n):
    usable = df[df["recommendation"].isin(["strong_candidate", "candidate"])].copy()
    if usable.empty:
        usable = df[df["recommendation"].isin(["review_manually"])].copy()
    if usable.empty:
        usable = df.copy()
    return usable.sort_values(["final_score", "confidence_score"], ascending=False).head(top_n)


def main():
    parser = argparse.ArgumentParser(description="Prepare Phase 2 MD candidate inputs.")
    parser.add_argument("--ranked", default=DEFAULT_INPUT)
    parser.add_argument("--fasta", default=DEFAULT_FASTA)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    fasta_dir = out_dir / "mutant_fastas"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    wt_sequence = read_sequence(args.fasta)
    df = pd.read_csv(args.ranked)
    candidates = choose_candidates(df, args.top_n)
    candidates.to_csv(out_dir / "md_candidates.csv", index=False)

    spec_lines = []
    for _, row in candidates.iterrows():
        mutation = row.get("mutation", f"{row['wt']}{row['pdb_residue']}{row['mut']}")
        mut_sequence = mutate_sequence(wt_sequence, int(row["pos"]), row["mut"])
        fasta_path = fasta_dir / f"{mutation}.fasta"
        fasta_path.write_text(f">{mutation}\n{mut_sequence}\n")
        spec_lines.append(
            f"{mutation}\tpos={int(row['pos'])}\tpdb_residue={int(row['pdb_residue'])}"
            f"\twt={row['wt']}\tmut={row['mut']}\tfinal_score={row['final_score']:.3f}"
        )

    (out_dir / "mutation_specs.tsv").write_text("\n".join(spec_lines) + "\n")
    (out_dir / "run_notes.txt").write_text(
        "MD preparation complete.\n"
        "Use these candidates for focused simulation only after Phase 1 ranking.\n\n"
        "Suggested MD readouts:\n"
        "- RMSD against WT/reference structure\n"
        "- RMSF around epitope loop and mutated residue\n"
        "- retention of key contacts/disulfides/salt bridges\n"
        "- gross unfolding or local instability\n"
    )

    print(f"MD inputs written to {out_dir}")
    print(f"Candidate count: {len(candidates)}")


if __name__ == "__main__":
    main()
