"""
Prepare ProteinMPNN inputs for the v2 construct-based pipeline.

Shortlisting strategy (v2):
  PRIMARY  — strong/candidate mutations in MHR (s_protein_pos 99–169) or
             the "a" determinant loop (41–70), excluding disulfide/buried-charge
             mutations.  These are the immunologically relevant positions.
  SECONDARY — top-ranked strong_candidate mutations outside those regions
             (TM domain / N-terminal), for structural context.

Outputs (in construct_v2/phase2/proteinmpnn/inputs/):
  proteinmpnn_candidates.csv    — shortlisted mutations
  fixed_positions.jsonl         — positions fixed during ProteinMPNN redesign
  shortlist_candidates.fasta    — WT + each mutant sequence (for --score_only)
  run_notes.txt                 — ready-to-paste run command

Run from project root:
    python phase2/proteinmpnn/prepare_inputs_v2.py
"""

import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).parent.parent.parent

RANKED_CSV  = ROOT / "construct_v2/results/final_ranked.csv"
WT_FASTA    = ROOT / "construct_v2/data/cleaned/hbsag_construct.fasta"
PDB_PATH    = ROOT / "construct_v2/data/cleaned/alphafold_construct.pdb"
OUTPUT_DIR  = ROOT / "construct_v2/phase2/proteinmpnn/inputs"

# Immunologically important regions (s_protein_pos, 1-indexed)
MHR_START, MHR_END       = 99, 169   # MHR / a-determinant
ADET_LOOP_START, ADET_LOOP_END = 41, 70    # disordered loop (present in construct)

BAD_FLAGS   = {"possible_disulfide_loss", "buried_charge", "low_solubility"}

N_PRIMARY   = 15   # MHR + a-det loop candidates
N_SECONDARY = 8    # top global hits outside epitope


def read_wt(fasta_path: Path) -> str:
    for rec in SeqIO.parse(fasta_path, "fasta"):
        return str(rec.seq).upper()
    raise ValueError(f"No sequence in {fasta_path}")


def has_bad_flag(flags_str) -> bool:
    if pd.isna(flags_str) or str(flags_str).strip() in ("", "nan"):
        return False
    flags = {f.strip() for f in str(flags_str).split(";")}
    return bool(flags & BAD_FLAGS)


def choose_candidates(df: pd.DataFrame) -> pd.DataFrame:
    clean = df[
        df["recommendation"].isin(["strong_candidate", "candidate"])
        & ~df["flags"].apply(has_bad_flag)
    ].copy()

    in_primary_region = (
        clean["s_protein_pos"].between(MHR_START, MHR_END)
        | clean["s_protein_pos"].between(ADET_LOOP_START, ADET_LOOP_END)
    )
    primary   = clean[in_primary_region].sort_values("final_score", ascending=False).head(N_PRIMARY)
    secondary = clean[~in_primary_region].sort_values("final_score", ascending=False).head(N_SECONDARY)

    combined = pd.concat([primary, secondary]).drop_duplicates(subset=["pos", "mut"])
    combined = combined.sort_values("final_score", ascending=False).reset_index(drop=True)
    return combined


def build_fixed_positions(candidates: pd.DataFrame, n_residues: int = 226) -> dict:
    """All residues NOT in the candidate mutable set are fixed."""
    mutable = set(int(r) for r in candidates["pdb_residue"])
    fixed   = sorted(i for i in range(1, n_residues + 1) if i not in mutable)
    pdb_stem = PDB_PATH.stem  # "alphafold_construct"
    return {pdb_stem: {"A": fixed}}


def write_fasta(candidates: pd.DataFrame, wt_seq: str, out_path: Path):
    lines = []
    for _, row in candidates.iterrows():
        mut_seq = list(wt_seq)
        mut_seq[int(row["pos"])] = row["mut"]
        lines.append(f">{row['mutation']}")
        lines.append("".join(mut_seq))
    out_path.write_text("\n".join(lines) + "\n")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df         = pd.read_csv(RANKED_CSV)
    wt_seq     = read_wt(WT_FASTA)
    candidates = choose_candidates(df)

    candidates.to_csv(OUTPUT_DIR / "proteinmpnn_candidates.csv", index=False)

    fixed = build_fixed_positions(candidates, n_residues=len(wt_seq))
    with open(OUTPUT_DIR / "fixed_positions.jsonl", "w") as fh:
        fh.write(json.dumps(fixed) + "\n")

    write_fasta(candidates, wt_seq, OUTPUT_DIR / "shortlist_candidates.fasta")

    mutable_pos = sorted(int(r) for r in candidates["pdb_residue"].unique())
    mutable_jsonl = {PDB_PATH.stem: {"A": mutable_pos}}
    with open(OUTPUT_DIR / "mutable_positions.jsonl", "w") as fh:
        fh.write(json.dumps(mutable_jsonl) + "\n")

    run_cmd = (
        "# Score each shortlisted mutant against the AlphaFold backbone\n"
        "python phase2/proteinmpnn/Workflow/protein_mpnn_run.py \\\n"
        f"  --pdb_path {PDB_PATH.relative_to(ROOT)} \\\n"
        "  --pdb_path_chains A \\\n"
        "  --out_folder construct_v2/phase2/proteinmpnn/results/scoring/ \\\n"
        f"  --fixed_positions_jsonl {(OUTPUT_DIR / 'fixed_positions.jsonl').relative_to(ROOT)} \\\n"
        f"  --path_to_fasta {(OUTPUT_DIR / 'shortlist_candidates.fasta').relative_to(ROOT)} \\\n"
        "  --num_seq_per_target 1 \\\n"
        "  --score_only 1\n\n"
        "# Design new sequences (redesign mutable positions)\n"
        "python phase2/proteinmpnn/Workflow/protein_mpnn_run.py \\\n"
        f"  --pdb_path {PDB_PATH.relative_to(ROOT)} \\\n"
        "  --pdb_path_chains A \\\n"
        "  --out_folder construct_v2/phase2/proteinmpnn/results/design/ \\\n"
        f"  --fixed_positions_jsonl {(OUTPUT_DIR / 'fixed_positions.jsonl').relative_to(ROOT)} \\\n"
        "  --num_seq_per_target 20 \\\n"
        "  --sampling_temp 0.2\n"
    )

    (OUTPUT_DIR / "run_notes.txt").write_text(
        "ProteinMPNN v2 preparation complete.\n\n"
        f"PDB:            {PDB_PATH.relative_to(ROOT)}\n"
        f"Candidates:     {len(candidates)} mutations "
        f"({len(candidates[candidates['s_protein_pos'].between(MHR_START, MHR_END)])} in MHR, "
        f"{len(candidates[candidates['s_protein_pos'].between(ADET_LOOP_START, ADET_LOOP_END)])} in a-det loop, "
        f"{len(candidates[~candidates['s_protein_pos'].between(ADET_LOOP_START, MHR_END)])} outside)\n"
        f"Mutable positions: {mutable_pos}\n\n"
        "Fixed positions: all residues NOT in the mutable set.\n"
        "Critical cysteines (C65, C76, C90, C107, C137, C149) and antibody\n"
        "contacts are fixed because they were excluded from shortlisting.\n\n"
        "Commands:\n"
        + run_cmd
    )

    print(f"ProteinMPNN v2 inputs written to {OUTPUT_DIR}")
    print(f"Total candidates: {len(candidates)}")
    print()
    print(candidates[["rank", "mutation", "s_protein_pos", "wt", "mut",
                       "esm_score", "ddg", "final_score", "recommendation"]].to_string(index=False))


if __name__ == "__main__":
    main()
