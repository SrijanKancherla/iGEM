from pathlib import Path

import pandas as pd


INPUT_CSV = "results/3_accessibility.csv"
OUTPUT_CSV = "results/4_epitope_scores.csv"

MHR_START = 99
MHR_END = 169
MAJOR_EPITOPE_START = 124
MAJOR_EPITOPE_END = 147

# Contact-sensitive residues in the 9UBQ HBsAg antigenic loop.
# This is intentionally conservative until a full atom-distance contact map is added.
ANTIBODY_CONTACT_RESIDUES = {
    101, 102, 103, 104, 105, 106,
    120, 121, 122, 123, 124, 125, 126,
    127, 128, 129, 130, 131, 132, 133,
    134, 135, 136, 137, 138, 139, 140,
    141, 142, 143, 144, 145, 146, 147,
}


def conservative_change(wt, mut):
    groups = [
        {"A", "V", "L", "I", "M"},
        {"F", "Y", "W"},
        {"S", "T", "N", "Q"},
        {"D", "E"},
        {"K", "R", "H"},
    ]
    return any(wt in group and mut in group for group in groups)


def score_epitope(row):
    pdb_residue = int(row.get("pdb_residue", row["pos"]))
    wt = row["wt"]
    mut = row["mut"]
    surface = row.get("surface_class", "unknown")

    in_mhr = MHR_START <= pdb_residue <= MHR_END
    in_major = MAJOR_EPITOPE_START <= pdb_residue <= MAJOR_EPITOPE_END
    in_contact = pdb_residue in ANTIBODY_CONTACT_RESIDUES

    tcell_penalty = 0
    bcell_penalty = 0

    if in_major:
        tcell_penalty += 18
        bcell_penalty += 18
    elif in_mhr:
        bcell_penalty += 10

    if in_contact:
        bcell_penalty += 20

    if surface == "surface" and (in_mhr or in_major):
        bcell_penalty += 6

    if conservative_change(wt, mut):
        tcell_penalty *= 0.5
        bcell_penalty *= 0.5

    return pd.Series(
        {
            "in_mhr": in_mhr,
            "in_major_epitope": in_major,
            "in_antibody_contact": in_contact,
            "tcell_penalty": tcell_penalty,
            "bcell_penalty": bcell_penalty,
            "epitope_penalty": tcell_penalty + bcell_penalty,
        }
    )


def main():
    df = pd.read_csv(INPUT_CSV)
    if "pdb_residue" not in df.columns:
        from residue_mapping import residue_map_by_pos

        mapping = residue_map_by_pos()
        df["pdb_residue"] = df["pos"].map(lambda pos: mapping.get(int(pos), {}).get("pdb_residue", int(pos)))

    epitope_cols = df.apply(score_epitope, axis=1)
    df = pd.concat([df, epitope_cols], axis=1)
    Path(OUTPUT_CSV).parent.mkdir(exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Epitope preservation scoring complete: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
