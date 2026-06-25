from pathlib import Path

import pandas as pd
from Bio import SeqIO


INPUT_CSV = "results/3_accessibility.csv"
OUTPUT_CSV = "results/4_epitope_scores.csv"
WT_FASTA = "data/cleaned/hbsag.fasta"

MHR_START = 99
MHR_END = 169
MAJOR_EPITOPE_START = 124
MAJOR_EPITOPE_END = 147

ANTIBODY_CONTACT_RESIDUES = {
    101, 102, 103, 104, 105, 106,
    120, 121, 122, 123, 124, 125, 126,
    127, 128, 129, 130, 131, 132, 133,
    134, 135, 136, 137, 138, 139, 140,
    141, 142, 143, 144, 145, 146, 147,
}

AA_GROUPS = [
    {"A", "V", "L", "I", "M"},
    {"F", "Y", "W"},
    {"S", "T", "N", "Q"},
    {"D", "E"},
    {"K", "R", "H"},
]

HYDROPHOBIC = {"A", "V", "L", "I", "M", "F", "W", "Y"}
AROMATIC = {"F", "W", "Y"}
CHARGED = {"D", "E", "K", "R", "H"}
POLAR = {"S", "T", "N", "Q"}

MHC_I_ANCHORS = {
    "HLA-A02_like": {
        2: {"L", "M", "I", "V", "A"},
        9: {"V", "L", "I", "M", "A"},
    },
    "HLA-A03_like": {
        2: {"A", "V", "S", "T"},
        9: {"K", "R"},
    },
    "HLA-B07_like": {
        2: {"P"},
        9: {"L", "F", "W", "Y", "I", "V"},
    },
}


def read_sequence(path=WT_FASTA):
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq).upper()
    raise ValueError(f"No FASTA sequence found in {path}")


def conservative_change(wt, mut):
    return any(wt in group and mut in group for group in AA_GROUPS)


def mutate_sequence(sequence, pos, mut):
    chars = list(sequence)
    chars[int(pos)] = mut
    return "".join(chars)


def overlapping_windows(sequence, pos, size):
    pos = int(pos)
    start_min = max(0, pos - size + 1)
    start_max = min(pos, len(sequence) - size)
    for start in range(start_min, start_max + 1):
        yield start, sequence[start:start + size]


def score_mhc_i_window(peptide, anchor_spec):
    if len(peptide) != 9:
        return 0.0

    score = 0.0
    for one_indexed_pos, allowed in anchor_spec.items():
        if peptide[one_indexed_pos - 1] in allowed:
            score += 0.45

    if peptide[0] in {"Y", "F", "W", "L", "M"}:
        score += 0.05
    if any(aa in {"P", "G"} for aa in peptide[3:7]):
        score -= 0.10
    if peptide.count("C") > 1:
        score -= 0.10

    return max(0.0, min(score, 1.0))


def max_mhc_i_score(sequence, pos):
    best_score = 0.0
    best_allele = None
    best_peptide = None

    for _, peptide in overlapping_windows(sequence, pos, 9):
        for allele, anchors in MHC_I_ANCHORS.items():
            score = score_mhc_i_window(peptide, anchors)
            if score > best_score:
                best_score = score
                best_allele = allele
                best_peptide = peptide

    return best_score, best_allele, best_peptide


def score_mhc_ii_core(core):
    if len(core) != 9:
        return 0.0

    score = 0.0
    for idx in [0, 3, 5, 8]:
        if core[idx] in HYDROPHOBIC | AROMATIC:
            score += 0.18
        if core[idx] in CHARGED:
            score -= 0.05

    if any(aa in {"P", "G"} for aa in core[2:7]):
        score -= 0.08

    return max(0.0, min(score, 1.0))


def max_mhc_ii_score(sequence, pos):
    best_score = 0.0
    best_peptide = None

    for _, peptide in overlapping_windows(sequence, pos, 15):
        for core_start in range(0, len(peptide) - 8):
            core = peptide[core_start:core_start + 9]
            score = score_mhc_ii_core(core)
            if score > best_score:
                best_score = score
                best_peptide = peptide

    return best_score, best_peptide


def kolaskar_like_antigenicity(sequence, pos):
    scores = {
        "A": 1.064, "C": 1.412, "D": 0.866, "E": 0.851, "F": 1.091,
        "G": 0.874, "H": 1.105, "I": 1.152, "K": 0.930, "L": 1.250,
        "M": 0.826, "N": 0.776, "P": 1.064, "Q": 1.015, "R": 0.873,
        "S": 1.012, "T": 0.909, "V": 1.383, "W": 0.893, "Y": 1.161,
    }
    pos = int(pos)
    start = max(0, pos - 3)
    end = min(len(sequence), pos + 4)
    window = sequence[start:end]
    return sum(scores.get(aa, 1.0) for aa in window) / max(len(window), 1)


def score_row(row, wt_sequence):
    pos = int(row["pos"])
    pdb_residue = int(row.get("pdb_residue", pos))
    wt = row["wt"]
    mut = row["mut"]
    surface = row.get("surface_class", "unknown")

    mut_sequence = mutate_sequence(wt_sequence, pos, mut)

    wt_mhci, wt_mhci_allele, wt_mhci_peptide = max_mhc_i_score(wt_sequence, pos)
    mut_mhci, mut_mhci_allele, mut_mhci_peptide = max_mhc_i_score(mut_sequence, pos)
    wt_mhcii, wt_mhcii_peptide = max_mhc_ii_score(wt_sequence, pos)
    mut_mhcii, mut_mhcii_peptide = max_mhc_ii_score(mut_sequence, pos)

    wt_bcell = kolaskar_like_antigenicity(wt_sequence, pos)
    mut_bcell = kolaskar_like_antigenicity(mut_sequence, pos)
    bcell_delta = mut_bcell - wt_bcell

    in_mhr = MHR_START <= pdb_residue <= MHR_END
    in_major = MAJOR_EPITOPE_START <= pdb_residue <= MAJOR_EPITOPE_END
    in_contact = pdb_residue in ANTIBODY_CONTACT_RESIDUES

    wt_tcell_score = max(wt_mhci, wt_mhcii)
    mut_tcell_score = max(mut_mhci, mut_mhcii)
    tcell_delta = mut_tcell_score - wt_tcell_score

    tcell_penalty = 0.0
    bcell_penalty = 0.0

    if wt_tcell_score >= 0.70 and tcell_delta <= -0.20:
        tcell_penalty += 18
    elif wt_tcell_score >= 0.55 and tcell_delta <= -0.15:
        tcell_penalty += 10

    if mut_tcell_score >= 0.75 and tcell_delta >= 0.25:
        tcell_penalty += 4

    if in_major:
        tcell_penalty += 12
        bcell_penalty += 16
    elif in_mhr:
        bcell_penalty += 8

    if in_contact:
        bcell_penalty += 20

    if surface == "surface" and (in_mhr or in_major):
        bcell_penalty += 5

    if bcell_delta <= -0.15 and (surface == "surface" or in_mhr):
        bcell_penalty += 8

    if conservative_change(wt, mut):
        tcell_penalty *= 0.6
        bcell_penalty *= 0.6

    return pd.Series(
        {
            "in_mhr": in_mhr,
            "in_major_epitope": in_major,
            "in_antibody_contact": in_contact,
            "wt_mhci_score": round(wt_mhci, 3),
            "mut_mhci_score": round(mut_mhci, 3),
            "mhci_delta": round(mut_mhci - wt_mhci, 3),
            "wt_mhci_allele": wt_mhci_allele,
            "mut_mhci_allele": mut_mhci_allele,
            "wt_mhci_peptide": wt_mhci_peptide,
            "mut_mhci_peptide": mut_mhci_peptide,
            "wt_mhcii_score": round(wt_mhcii, 3),
            "mut_mhcii_score": round(mut_mhcii, 3),
            "mhcii_delta": round(mut_mhcii - wt_mhcii, 3),
            "wt_mhcii_peptide": wt_mhcii_peptide,
            "mut_mhcii_peptide": mut_mhcii_peptide,
            "wt_tcell_epitope_score": round(wt_tcell_score, 3),
            "mut_tcell_epitope_score": round(mut_tcell_score, 3),
            "tcell_epitope_delta": round(tcell_delta, 3),
            "wt_bcell_antigenicity": round(wt_bcell, 3),
            "mut_bcell_antigenicity": round(mut_bcell, 3),
            "bcell_antigenicity_delta": round(bcell_delta, 3),
            "tcell_penalty": round(tcell_penalty, 3),
            "bcell_penalty": round(bcell_penalty, 3),
            "epitope_penalty": round(tcell_penalty + bcell_penalty, 3),
            "epitope_prediction_mode": "internal_motif_model",
        }
    )


def main():
    df = pd.read_csv(INPUT_CSV)
    if "pdb_residue" not in df.columns:
        from residue_mapping import residue_map_by_pos

        mapping = residue_map_by_pos()
        df["pdb_residue"] = df["pos"].map(
            lambda pos: mapping.get(int(pos), {}).get("pdb_residue", int(pos))
        )

    wt_sequence = read_sequence()
    scored = df.apply(lambda row: score_row(row, wt_sequence), axis=1)
    df = pd.concat([df, scored], axis=1)
    Path(OUTPUT_CSV).parent.mkdir(exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Epitope prediction scoring complete: {OUTPUT_CSV}")
    print("Epitope prediction mode: internal_motif_model")


if __name__ == "__main__":
    main()
