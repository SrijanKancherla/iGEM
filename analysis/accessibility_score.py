from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley

from residue_mapping import residue_map_by_pos


INPUT_CSV = "results/2_conservation.csv"
OUTPUT_CSV = "results/3_accessibility.csv"
PDB_PATH = "data/cleaned/hbsag.pdb"

MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}


def surface_class(rsa):
    if rsa < 0.15:
        return "buried"
    if rsa < 0.35:
        return "partial"
    return "surface"


def mutation_accessibility_penalty(row):
    cls = row["surface_class"]
    wt = row["wt"]
    mut = row["mut"]

    penalty = 0
    if cls == "buried" and mut in {"D", "E", "K", "R", "H"}:
        penalty += 12
    if cls == "buried" and mut == "P":
        penalty += 10
    if cls == "buried" and wt in {"F", "I", "L", "M", "V", "W", "Y"} and mut in {"D", "E", "K", "R", "N", "Q"}:
        penalty += 8
    if cls == "surface" and mut in {"F", "I", "L", "M", "V", "W", "Y"}:
        penalty += 4

    return penalty


def calculate_accessibility(pdb_path=PDB_PATH):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    ShrakeRupley().compute(structure, level="R")

    mapping = residue_map_by_pos(pdb_path)
    out = {}

    residues = []
    for model in structure:
        if "A" not in model:
            break
        residues = [r for r in model["A"] if r.id[0] == " " and r.has_id("CA")]
        break

    for pos, residue in enumerate(residues):
        wt = None
        asa = float(getattr(residue, "sasa", 0.0) or 0.0)
        if pos in mapping:
            wt = mapping[pos]
        out[pos] = {
            "asa": asa,
            "rsa": None,
            "surface_class": "unknown",
        }

    return out


def main():
    df = pd.read_csv(INPUT_CSV)
    access = calculate_accessibility()

    df["asa"] = df["pos"].map(lambda pos: access.get(int(pos), {}).get("asa", 0.0))
    df["rsa"] = df.apply(
        lambda row: min(float(row["asa"]) / MAX_ASA.get(row["wt"], 200.0), 1.0),
        axis=1,
    )
    df["surface_class"] = df["rsa"].apply(surface_class)
    df["accessibility_penalty"] = df.apply(mutation_accessibility_penalty, axis=1)

    Path(OUTPUT_CSV).parent.mkdir(exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Accessibility scoring complete: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
