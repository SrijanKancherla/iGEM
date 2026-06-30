"""
Solvent-accessibility scoring for the v2 construct-based pipeline.

Requires: construct_v2/data/cleaned/alphafold_construct.pdb
  AlphaFold residue numbers are 1-indexed and map directly to s_protein_pos (pos + 1).

Falls back to per-class penalties if the PDB is absent, so the pipeline can
continue without a structure (with reduced accuracy for buried/surface calls).
"""

from pathlib import Path
import sys

import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley

INPUT_CSV  = Path("construct_v2/results/2_conservation.csv")
OUTPUT_CSV = Path("construct_v2/results/3_accessibility.csv")
PDB_V2     = Path("construct_v2/data/cleaned/alphafold_construct.pdb")

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
    wt  = row["wt"]
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


def calculate_accessibility_from_pdb(pdb_path):
    """Return {0-indexed pos: {asa, rsa, surface_class}} from AlphaFold PDB."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    ShrakeRupley().compute(structure, level="R")

    out = {}
    for model in structure:
        chain = model["A"] if "A" in model else list(model.get_chains())[0]
        residues = [r for r in chain if r.id[0] == " " and r.has_id("CA")]
        for r in residues:
            # AlphaFold PDB: residue sequence number = s_protein_pos = pos + 1
            pos = int(r.id[1]) - 1   # convert back to 0-indexed
            wt_aa = None
            try:
                from Bio.PDB.Polypeptide import protein_letters_3to1
                wt_aa = protein_letters_3to1.get(r.resname, "X")
            except ImportError:
                pass
            asa = float(getattr(r, "sasa", 0.0) or 0.0)
            out[pos] = {"asa": asa}
        break
    return out


def main():
    df = pd.read_csv(INPUT_CSV)

    if PDB_V2.exists():
        print(f"Computing accessibility from {PDB_V2}")
        access = calculate_accessibility_from_pdb(PDB_V2)
        df["asa"] = df["pos"].map(lambda p: access.get(int(p), {}).get("asa", 0.0))
        df["rsa"] = df.apply(
            lambda r: min(float(r["asa"]) / MAX_ASA.get(r["wt"], 200.0), 1.0),
            axis=1,
        )
        df["surface_class"] = df["rsa"].apply(surface_class)
    else:
        print(f"WARNING: {PDB_V2} not found — using 'unknown' surface class.")
        print("  Run AlphaFold on construct_v2/alphafold_input/hbsag_for_alphafold.fasta")
        print("  and save result to construct_v2/data/cleaned/alphafold_construct.pdb")
        df["asa"] = 0.0
        df["rsa"] = 0.5        # neutral: won't over-penalise any mutation
        df["surface_class"] = "unknown"

    df["accessibility_penalty"] = df.apply(mutation_accessibility_penalty, axis=1)
    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Accessibility scoring complete → {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
