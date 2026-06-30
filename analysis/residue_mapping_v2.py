"""
Residue mapping for the v2 construct-based pipeline.

No PDB required. The mapping is a fixed arithmetic offset:
  s_protein_pos  = pos + 1      (pos is 0-indexed in the 226 aa HBsAg sequence)
  construct_pos  = pos + 13     (1-indexed in the full 238 aa construct including His-tag)

This replaces analysis/residue_mapping.py, which derived positions from the
9UBQ ATOM records and had a non-trivial gap (PDB residues 42-71 unresolved).
"""

from pathlib import Path


HBSAG_FASTA_V2 = Path("construct_v2/data/cleaned/hbsag_construct.fasta")
N_RESIDUES = 226
HIS_TAG_LEN = 12  # MGSSHHHHHHSS


def build_residue_mapping_v2(n_residues: int = N_RESIDUES):
    """Return list of dicts mapping 0-indexed pos → s_protein_pos / construct_pos."""
    return [
        {
            "pos": i,
            "chain": "A",
            "s_protein_pos": i + 1,
            "pdb_residue": i + 1,          # alias kept for backward compat with scoring scripts
            "construct_pos": i + HIS_TAG_LEN + 1,
        }
        for i in range(n_residues)
    ]


def residue_map_by_pos_v2(n_residues: int = N_RESIDUES) -> dict:
    """Return {pos: mapping_dict} keyed by 0-indexed position."""
    return {row["pos"]: row for row in build_residue_mapping_v2(n_residues)}


# Backward-compatible aliases so combine_scores_v2 can call the same interface
def build_residue_mapping(pdb_path=None, chain_id="A"):
    return build_residue_mapping_v2()


def residue_map_by_pos(pdb_path=None, chain_id="A"):
    return residue_map_by_pos_v2()


if __name__ == "__main__":
    import pandas as pd

    mapping = build_residue_mapping_v2()
    df = pd.DataFrame(mapping)
    out = Path("construct_v2/results/residue_mapping_v2.csv")
    out.parent.mkdir(exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Saved residue mapping for {len(df)} residues → {out}")
    print(df.head(10).to_string(index=False))
