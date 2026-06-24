from pathlib import Path

from Bio.PDB import PDBParser


PDB_PATH = Path("data/cleaned/hbsag.pdb")


def build_residue_mapping(pdb_path=PDB_PATH, chain_id="A"):
    """Map zero-indexed pipeline positions to PDB residue ids."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))

    residues = []
    for model in structure:
        if chain_id not in model:
            continue

        for residue in model[chain_id]:
            hetflag, resseq, icode = residue.id
            if hetflag == " " and residue.has_id("CA"):
                residues.append(
                    {
                        "pos": len(residues),
                        "chain": chain_id,
                        "pdb_residue": int(resseq),
                        "insertion_code": icode.strip(),
                    }
                )
        break

    return residues


def residue_map_by_pos(pdb_path=PDB_PATH, chain_id="A"):
    return {row["pos"]: row for row in build_residue_mapping(pdb_path, chain_id)}


if __name__ == "__main__":
    import pandas as pd

    mapping = build_residue_mapping()
    df = pd.DataFrame(mapping)
    Path("results").mkdir(exist_ok=True)
    df.to_csv("results/residue_mapping.csv", index=False)
    print(f"Saved residue mapping for {len(df)} residues to results/residue_mapping.csv")
