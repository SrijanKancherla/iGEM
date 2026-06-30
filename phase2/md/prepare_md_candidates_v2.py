"""
Prepare Phase 2 MD candidate inputs for the v2 construct-based pipeline.

Shortlisting strategy (v2):
  PRIMARY  — mutations in MHR (s_protein_pos 99–169) or the "a" determinant
             loop (41–70), no disulfide-loss / buried-charge flags.
  SECONDARY — top-ranked strong_candidate mutations outside those regions.

Outputs (in construct_v2/phase2/md/inputs/):
  md_candidates.csv           — shortlisted mutations with all scores
  mutation_specs.tsv          — concise spec sheet (construct_pos included)
  mutant_fastas/<MUT>.fasta   — 226 aa mutant sequence for each candidate

Run from project root:
    python phase2/md/prepare_md_candidates_v2.py
"""

from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).parent.parent.parent

RANKED_CSV = ROOT / "construct_v2/results/final_ranked.csv"
WT_FASTA   = ROOT / "construct_v2/data/cleaned/hbsag_construct.fasta"
OUTPUT_DIR = ROOT / "construct_v2/phase2/md/inputs"

MHR_START, MHR_END           = 99, 169
ADET_LOOP_START, ADET_LOOP_END = 41, 70

BAD_FLAGS  = {"possible_disulfide_loss", "buried_charge", "low_solubility"}

N_PRIMARY   = 15
N_SECONDARY = 8


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

    in_primary = (
        clean["s_protein_pos"].between(MHR_START, MHR_END)
        | clean["s_protein_pos"].between(ADET_LOOP_START, ADET_LOOP_END)
    )
    primary   = clean[in_primary].sort_values("final_score", ascending=False).head(N_PRIMARY)
    secondary = clean[~in_primary].sort_values("final_score", ascending=False).head(N_SECONDARY)

    combined = pd.concat([primary, secondary]).drop_duplicates(subset=["pos", "mut"])
    return combined.sort_values("final_score", ascending=False).reset_index(drop=True)


def apply_mutation(wt_seq: str, pos: int, mut: str) -> str:
    chars = list(wt_seq)
    chars[pos] = mut
    return "".join(chars)


def main():
    fasta_dir = OUTPUT_DIR / "mutant_fastas"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    wt_seq     = read_wt(WT_FASTA)
    df         = pd.read_csv(RANKED_CSV)
    candidates = choose_candidates(df)

    candidates.to_csv(OUTPUT_DIR / "md_candidates.csv", index=False)

    spec_lines = []
    for _, row in candidates.iterrows():
        mutation    = row["mutation"]
        pos         = int(row["pos"])
        s_pos       = int(row["s_protein_pos"])
        c_pos       = int(row["construct_pos"]) if "construct_pos" in row else s_pos + 12
        mut_seq     = apply_mutation(wt_seq, pos, row["mut"])
        fasta_path  = fasta_dir / f"{mutation}.fasta"
        fasta_path.write_text(f">{mutation}\n{mut_seq}\n")
        spec_lines.append(
            f"{mutation}"
            f"\tpos={pos}"
            f"\ts_protein_pos={s_pos}"
            f"\tconstruct_pos={c_pos}"
            f"\twt={row['wt']}"
            f"\tmut={row['mut']}"
            f"\tfinal_score={row['final_score']:.3f}"
            f"\tconfidence={row['confidence_score']:.3f}"
        )

    (OUTPUT_DIR / "mutation_specs.tsv").write_text("\n".join(spec_lines) + "\n")
    (OUTPUT_DIR / "run_notes.txt").write_text(
        "MD v2 preparation complete.\n"
        "Sequences are based on the 226 aa expression construct (no His-tag).\n"
        "construct_pos = s_protein_pos + 12 (accounts for 12-aa His-tag in the\n"
        "full 238 aa CDS; use s_protein_pos for structure-based numbering).\n\n"
        "Suggested MD readouts:\n"
        "- RMSD of mutant vs WT AlphaFold backbone (construct_v2/data/cleaned/alphafold_construct.pdb)\n"
        "- RMSF around mutated residue ± 10 residues\n"
        "- RMSF of 'a' determinant loop (s_protein_pos 99–169)\n"
        "- Monitoring of C76–C221 and C90–C65 disulfide distances (confirm intact)\n"
        "- Gross unfolding index (fraction of contacts lost vs WT)\n"
    )

    print(f"MD v2 inputs written to {OUTPUT_DIR}")
    print(f"Total candidates: {len(candidates)}")
    print(f"FASTA files written: {len(list(fasta_dir.glob('*.fasta')))}")
    print()

    display_cols = ["mutation", "s_protein_pos", "wt", "mut",
                    "esm_score", "ddg", "final_score", "recommendation"]
    display_cols = [c for c in display_cols if c in candidates.columns]
    print(candidates[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
