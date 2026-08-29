from pathlib import Path

import pandas as pd


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
FASTA_PATH = "data/cleaned/hbsag.fasta"
OUTPUT_PATH = "data/mutations/esm2_mutations.csv"


def read_fasta(path):
    sequence = ""
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()


def heuristic_score(wt, mut):
    score = 0.0
    if mut in {"K", "R", "E", "D"}:
        score -= 0.2
    if mut == "C":
        score -= 1.0
    if mut == "P":
        score -= 0.8
    if wt == "C" and mut != "C":
        score -= 0.6
    return score


def main():
    sequence = read_fasta(FASTA_PATH)
    mutations = []

    for pos, wt in enumerate(sequence):
        for mut in AMINO_ACIDS:
            if mut == wt:
                continue

            mutations.append(
                {
                    "pos": pos,
                    "wt": wt,
                    "mut": mut,
                    "esm_score": heuristic_score(wt, mut),
                }
            )

    df = pd.DataFrame(mutations).sort_values("esm_score", ascending=False)
    Path(OUTPUT_PATH).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)
    print(f"Generated heuristic mutation candidates: {len(df)}")
    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
