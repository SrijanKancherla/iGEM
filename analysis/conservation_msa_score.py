from collections import Counter
from math import log2
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from scoring_utils import conservation_penalty


INPUT_CSV = "results/1_filtered.csv"
OUTPUT_CSV = "results/2_conservation.csv"
MSA_PATH = Path("data/msa/hbsag_msa.fasta")
WT_FASTA = Path("data/cleaned/hbsag.fasta")

AA_GROUPS = [
    {"A", "V", "L", "I", "M"},
    {"F", "Y", "W"},
    {"S", "T", "N", "Q"},
    {"D", "E"},
    {"K", "R", "H"},
    {"G"},
    {"P"},
    {"C"},
]


def read_sequence(path):
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq).upper().replace("-", "")
    raise ValueError(f"No FASTA sequence found in {path}")


def is_conservative(wt, mut):
    return any(wt in group and mut in group for group in AA_GROUPS)


def normalized_entropy(counts):
    total = sum(counts.values())
    if total == 0:
        return 0.0

    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * log2(p)

    return entropy / log2(20)


def load_msa_profile(msa_path=MSA_PATH, wt_fasta=WT_FASTA):
    if not msa_path.exists():
        return None, "fallback_no_msa"

    records = list(SeqIO.parse(msa_path, "fasta"))
    if len(records) < 2:
        return None, "fallback_too_few_msa_sequences"

    aligned = [str(record.seq).upper() for record in records]
    lengths = {len(seq) for seq in aligned}
    if len(lengths) != 1:
        return None, "fallback_unequal_msa_lengths"

    wt_sequence = read_sequence(wt_fasta)
    reference = aligned[0]
    ungapped_pos = -1
    profile = {}

    for aln_col, ref_aa in enumerate(reference):
        if ref_aa == "-":
            continue

        ungapped_pos += 1
        if ungapped_pos >= len(wt_sequence):
            break

        column = [
            seq[aln_col]
            for seq in aligned
            if seq[aln_col] not in {"-", "X", "B", "Z", "J", "U", "O"}
        ]
        counts = Counter(column)
        total = sum(counts.values())
        wt = wt_sequence[ungapped_pos]

        if total == 0:
            wt_frequency = 0.0
            consensus = wt
            entropy = 0.0
        else:
            wt_frequency = counts.get(wt, 0) / total
            consensus = counts.most_common(1)[0][0]
            entropy = normalized_entropy(counts)

        profile[ungapped_pos] = {
            "alignment_column": aln_col,
            "msa_depth": len(records),
            "msa_observed": total,
            "msa_consensus": consensus,
            "wt_frequency": wt_frequency,
            "position_entropy": entropy,
            "conservation_score_raw": 1.0 - entropy,
            "conservation_mode": "msa",
            "counts": counts,
        }

    if not profile:
        return None, "fallback_empty_profile"

    return profile, "msa"


def score_row(row, profile, mode):
    pos = int(row["pos"])
    wt = row["wt"]
    mut = row["mut"]

    if profile is None or pos not in profile:
        penalty = conservation_penalty(pos)
        return pd.Series(
            {
                "alignment_column": None,
                "msa_depth": 0,
                "msa_observed": 0,
                "msa_consensus": None,
                "wt_frequency": None,
                "mut_frequency": None,
                "position_entropy": None,
                "conservation_score_raw": 1.0 if penalty else 0.0,
                "conservation_penalty": penalty,
                "conservation_mode": mode,
            }
        )

    entry = profile[pos]
    counts = entry["counts"]
    observed = max(entry["msa_observed"], 1)
    mut_frequency = counts.get(mut, 0) / observed
    wt_frequency = entry["wt_frequency"]
    entropy = entry["position_entropy"]

    penalty = 0.0
    if wt_frequency >= 0.90:
        penalty += 18
    elif wt_frequency >= 0.70:
        penalty += 10
    elif wt_frequency >= 0.50:
        penalty += 5

    if mut_frequency == 0:
        penalty += 6
    elif mut_frequency < 0.05:
        penalty += 3

    if not is_conservative(wt, mut):
        penalty += 6

    if entropy >= 0.55:
        penalty *= 0.5

    return pd.Series(
        {
            "alignment_column": entry["alignment_column"],
            "msa_depth": entry["msa_depth"],
            "msa_observed": entry["msa_observed"],
            "msa_consensus": entry["msa_consensus"],
            "wt_frequency": wt_frequency,
            "mut_frequency": mut_frequency,
            "position_entropy": entropy,
            "conservation_score_raw": entry["conservation_score_raw"],
            "conservation_penalty": round(penalty, 3),
            "conservation_mode": entry["conservation_mode"],
        }
    )


def main():
    df = pd.read_csv(INPUT_CSV)
    profile, mode = load_msa_profile()
    scored = df.apply(lambda row: score_row(row, profile, mode), axis=1)
    df = pd.concat([df, scored], axis=1)
    Path(OUTPUT_CSV).parent.mkdir(exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Conservation scoring complete: {OUTPUT_CSV}")
    print(f"Conservation mode: {mode}")


if __name__ == "__main__":
    main()
