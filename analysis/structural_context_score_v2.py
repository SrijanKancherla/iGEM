"""
Structural context scoring for the v2 construct-based pipeline.

Identical logic to v1 but uses construct_v2/results/ paths.
Critical change: wt=="C" now correctly triggers disulfide-loss flags for
mutations at s_protein positions 76 and 90, which were incorrectly labelled
as alanine in v1.
"""

from pathlib import Path
import pandas as pd

INPUT_CSV  = Path("construct_v2/results/4_epitope_scores.csv")
OUTPUT_CSV = Path("construct_v2/results/5_structural_context.csv")

HYDROPHOBIC = {"A","V","L","I","M","F","W","Y"}
CHARGED     = {"D","E","K","R","H"}
POLAR       = {"S","T","N","Q","C"}


def score_structural_context(row):
    wt           = row["wt"]
    mut          = row["mut"]
    surface_class = row.get("surface_class", "unknown")
    in_contact   = bool(row.get("in_antibody_contact", False))

    bonding_penalty = interface_penalty = structural_penalty = 0
    flags = []

    if wt == "C" and mut != "C":
        bonding_penalty += 25
        flags.append("possible_disulfide_loss")
    if mut == "C" and wt != "C":
        bonding_penalty += 10
        flags.append("new_cysteine")

    if surface_class == "buried" and mut in CHARGED:
        structural_penalty += 16
        flags.append("buried_charge")
    if surface_class == "buried" and wt in HYDROPHOBIC and mut in POLAR | CHARGED:
        structural_penalty += 10
        flags.append("buried_hydrophobic_loss")
    if surface_class == "surface" and mut in HYDROPHOBIC:
        structural_penalty += 4
        flags.append("surface_hydrophobic_gain")

    if mut == "P":
        structural_penalty += 8
        flags.append("new_proline")
    if mut == "G":
        structural_penalty += 4
        flags.append("new_glycine")

    if in_contact:
        interface_penalty += 15
        flags.append("antibody_interface_mutation")

    return pd.Series({
        "bonding_penalty":            bonding_penalty,
        "interface_penalty":          interface_penalty,
        "structural_context_penalty": bonding_penalty + interface_penalty + structural_penalty,
        "structural_flags":           ";".join(flags),
    })


def main():
    df = pd.read_csv(INPUT_CSV)
    scored = df.apply(score_structural_context, axis=1)
    df = pd.concat([df, scored], axis=1)
    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Structural context scoring complete → {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
