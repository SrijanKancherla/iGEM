"""
MSA-based conservation scoring for the v2 construct-based pipeline.

Uses construct_v2/data/msa/hbsag_construct_msa.fasta (build with scripts/build_msa_v2.py).
If the v2 MSA is absent, falls back to the existing v1 MSA with a positional remap
so the overlapping 113-aa region still contributes conservation data.
"""

from collections import Counter
from math import log2
from pathlib import Path
import sys

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from scoring_utils import conservation_penalty

INPUT_CSV   = Path("construct_v2/results/1_filtered.csv")
OUTPUT_CSV  = Path("construct_v2/results/2_conservation.csv")
WT_FASTA_V2 = Path("construct_v2/data/cleaned/hbsag_construct.fasta")
MSA_V2      = Path("construct_v2/data/msa/hbsag_construct_msa.fasta")
MSA_V1      = Path("data/msa/hbsag_msa.fasta")          # fallback
WT_FASTA_V1 = Path("data/cleaned/hbsag.fasta")           # fallback reference

AA_GROUPS = [
    {"A", "V", "L", "I", "M"}, {"F", "Y", "W"}, {"S", "T", "N", "Q"},
    {"D", "E"}, {"K", "R", "H"}, {"G"}, {"P"}, {"C"},
]


def read_sequence(path):
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq).upper().replace("-", "")
    raise ValueError(f"No FASTA in {path}")


def is_conservative(wt, mut):
    return any(wt in g and mut in g for g in AA_GROUPS)


def normalized_entropy(counts):
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = -sum((c / total) * log2(c / total) for c in counts.values())
    return entropy / log2(20)


def build_profile(msa_path, wt_seq):
    """Build per-position conservation profile aligned to wt_seq."""
    records = list(SeqIO.parse(msa_path, "fasta"))
    if len(records) < 2:
        return None, "fallback_too_few_msa_sequences"
    aligned = [str(r.seq).upper() for r in records]
    if len(set(len(s) for s in aligned)) != 1:
        return None, "fallback_unequal_msa_lengths"

    reference = aligned[0]
    profile = {}
    ungapped_pos = -1

    for col, ref_aa in enumerate(reference):
        if ref_aa == "-":
            continue
        ungapped_pos += 1
        if ungapped_pos >= len(wt_seq):
            break
        column = [s[col] for s in aligned
                  if s[col] not in {"-", "X", "B", "Z", "J", "U", "O"}]
        counts = Counter(column)
        total = sum(counts.values())
        wt = wt_seq[ungapped_pos]
        wt_frequency = counts.get(wt, 0) / max(total, 1)
        entropy = normalized_entropy(counts)
        profile[ungapped_pos] = {
            "alignment_column": col,
            "msa_depth": len(records),
            "msa_observed": total,
            "msa_consensus": counts.most_common(1)[0][0] if total else wt,
            "wt_frequency": wt_frequency,
            "position_entropy": entropy,
            "conservation_score_raw": 1.0 - entropy,
            "conservation_mode": "msa",
            "counts": counts,
        }
    return (profile, "msa") if profile else (None, "fallback_empty_profile")


def build_v1_remap_profile(wt_v2_seq):
    """
    Use the v1 MSA (113 aa reference) and remap positions to the 226 aa v2 sequence.
    Alignment: v1 WT sequence vs v2 WT sequence, then project MSA columns.
    """
    if not MSA_V1.exists() or not WT_FASTA_V1.exists():
        return None, "fallback_no_v1_msa"

    wt_v1 = read_sequence(WT_FASTA_V1)
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5

    try:
        aln = next(aligner.align(wt_v2_seq, wt_v1))
        gapped_v2 = str(aln[0])
        gapped_v1 = str(aln[1])
    except Exception:
        return None, "fallback_v1_alignment_failed"

    # Map v2 0-indexed positions to v1 0-indexed positions
    v2_to_v1 = {}
    v2_idx = v1_idx = 0
    for a, b in zip(gapped_v2, gapped_v1):
        if a != "-" and b != "-":
            v2_to_v1[v2_idx] = v1_idx
        if a != "-":
            v2_idx += 1
        if b != "-":
            v1_idx += 1

    v1_profile, mode = build_profile(MSA_V1, wt_v1)
    if v1_profile is None:
        return None, mode

    # Reindex: v2 position → v1 conservation data
    profile = {}
    for v2_pos, v1_pos in v2_to_v1.items():
        if v1_pos in v1_profile:
            entry = dict(v1_profile[v1_pos])
            # recalculate wt_frequency using v2 WT at this position
            wt_v2_aa = wt_v2_seq[v2_pos]
            counts = entry["counts"]
            total = max(entry["msa_observed"], 1)
            entry["wt_frequency"] = counts.get(wt_v2_aa, 0) / total
            entry["conservation_mode"] = "msa_v1_remapped"
            profile[v2_pos] = entry

    return (profile, "msa_v1_remapped") if profile else (None, "fallback_v1_remap_empty")


def load_msa_profile():
    wt_seq = read_sequence(WT_FASTA_V2)
    if MSA_V2.exists():
        print(f"Using v2 MSA: {MSA_V2}")
        return build_profile(MSA_V2, wt_seq)
    print(f"v2 MSA not found ({MSA_V2}). Attempting v1 remap from {MSA_V1} ...")
    print("  Run 'python scripts/build_msa_v2.py' for full coverage.")
    profile, mode = build_v1_remap_profile(wt_seq)
    if profile:
        print(f"  v1 remap succeeded: {len(profile)} positions with conservation data.")
    return profile, mode


def score_row(row, profile, mode, wt_seq):
    pos = int(row["pos"])
    wt  = row["wt"]
    mut = row["mut"]

    if profile is None or pos not in profile:
        penalty = conservation_penalty(pos)
        return pd.Series({
            "alignment_column": None, "msa_depth": 0, "msa_observed": 0,
            "msa_consensus": None, "wt_frequency": None, "mut_frequency": None,
            "position_entropy": None, "conservation_score_raw": 1.0 if penalty else 0.0,
            "conservation_penalty": penalty, "conservation_mode": mode,
        })

    entry = profile[pos]
    counts = entry["counts"]
    observed = max(entry["msa_observed"], 1)
    mut_frequency = counts.get(mut, 0) / observed
    wt_frequency  = entry["wt_frequency"]
    entropy       = entry["position_entropy"]

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

    return pd.Series({
        "alignment_column":      entry["alignment_column"],
        "msa_depth":             entry["msa_depth"],
        "msa_observed":          entry["msa_observed"],
        "msa_consensus":         entry["msa_consensus"],
        "wt_frequency":          wt_frequency,
        "mut_frequency":         mut_frequency,
        "position_entropy":      entropy,
        "conservation_score_raw": entry["conservation_score_raw"],
        "conservation_penalty":  round(penalty, 3),
        "conservation_mode":     entry["conservation_mode"],
    })


def main():
    wt_seq = read_sequence(WT_FASTA_V2)
    df = pd.read_csv(INPUT_CSV)
    profile, mode = load_msa_profile()
    scored = df.apply(lambda row: score_row(row, profile, mode, wt_seq), axis=1)
    df = pd.concat([df, scored], axis=1)
    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Conservation scoring complete → {OUTPUT_CSV}  (mode: {mode})")


if __name__ == "__main__":
    main()
