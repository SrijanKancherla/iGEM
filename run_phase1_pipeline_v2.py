"""
Phase 1 pipeline v2 — predictions on the actual expression construct.

Reference sequence: construct_v2/data/cleaned/hbsag_construct.fasta
  226 aa HBsAg, no His-tag, numbered so s_protein_pos = pos + 1.

Run from the project root:
    python run_phase1_pipeline_v2.py [--skip-esm] [--skip-rosetta] [--device cuda]

Pipeline steps
--------------
STEP 1  ESM2 mutation scoring          [NO PDB NEEDED — runs immediately]
STEP 2  Rosetta ddG scanning           [NEEDS AlphaFold PDB — see note below]
STEP 3  Stability filter (ddG < 5)     [depends on step 2]
STEP 4  MSA conservation scoring       [runs with v1 MSA fallback if v2 MSA absent]
STEP 5  Solvent accessibility          [NEEDS AlphaFold PDB — gracefully skipped]
STEP 6  Epitope preservation           [NO PDB NEEDED — runs on sequence]
STEP 7  Structural context             [uses surface_class from step 5]
STEP 8  Expression scoring             [NO PDB NEEDED]
STEP 9  Solubility scoring             [NO PDB NEEDED]
STEP 10 Final ranking                  [combines all above]

AlphaFold note
--------------
Steps 2 and 5 require a predicted structure.
1. Submit construct_v2/alphafold_input/hbsag_for_alphafold.fasta to AlphaFold Server
   (https://alphafoldserver.com) or run locally / via ColabFold.
2. Download the relaxed rank-1 model PDB.
3. Save it as: construct_v2/data/cleaned/alphafold_construct.pdb
Then re-run this script; steps 2 and 5 will activate automatically.

MSA note
--------
Run 'python scripts/build_msa_v2.py' once before step 4 for full conservation
data. Without it the script falls back to the v1 MSA (covers ~113 positions).
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent


def run(cmd: list, label: str):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    result = subprocess.run([sys.executable] + cmd, cwd=str(ROOT))
    if result.returncode != 0:
        print(f"\nERROR: {label} failed (exit {result.returncode}). Stopping.")
        sys.exit(result.returncode)


def check_file(path: Path, label: str) -> bool:
    if path.exists():
        return True
    print(f"\n  SKIP: {label}")
    print(f"        {path} not found.")
    return False


def main():
    parser = argparse.ArgumentParser(description="v2 Phase 1 pipeline on expression construct.")
    parser.add_argument("--skip-esm",    action="store_true", help="Skip ESM2 (reuse existing CSV)")
    parser.add_argument("--skip-rosetta", action="store_true", help="Skip Rosetta (reuse existing CSV)")
    parser.add_argument("--device",      default="cpu", choices=["cpu","cuda"])
    parser.add_argument("--top-k",       type=int, default=10)
    args = parser.parse_args()

    alphafold_pdb  = ROOT / "construct_v2/data/cleaned/alphafold_construct.pdb"
    esm_output     = ROOT / "construct_v2/mutations/esm2_mutations.csv"
    rosetta_output = ROOT / "construct_v2/results/rosetta_ddg.csv"

    # ── STEP 1: ESM2 ─────────────────────────────────────────────────────────
    if args.skip_esm and esm_output.exists():
        print("\nSKIP: ESM2 scoring (--skip-esm; reusing existing CSV)")
    else:
        run(
            ["esm/esm2_mutations_v2.py",
             "--fasta",  "construct_v2/data/cleaned/hbsag_construct.fasta",
             "--output", "construct_v2/mutations/esm2_mutations.csv",
             "--top-k",  str(args.top_k),
             "--device", args.device],
            "STEP 1 — ESM2 mutation scoring (construct v2)",
        )

    # ── STEP 2: Rosetta ddG ──────────────────────────────────────────────────
    if args.skip_rosetta and rosetta_output.exists():
        print("\nSKIP: Rosetta ddG (--skip-rosetta; reusing existing CSV)")
    elif not alphafold_pdb.exists():
        print(f"\n  STEP 2 — Rosetta ddG: SKIPPED (AlphaFold PDB not yet available)")
        print(f"    → Save AlphaFold result to: {alphafold_pdb}")
        print(f"    → Then re-run this script (use --skip-esm to save time)")
        if not rosetta_output.exists():
            print("\n  Cannot continue past step 2 without ddG data. Exiting.")
            print("  Provide the AlphaFold PDB and re-run, or supply a rosetta_ddg.csv manually.")
            sys.exit(0)
        print("  Found existing rosetta_ddg.csv — continuing with cached ddG values.")
    else:
        run(["rosetta/ddg_scan_v2.py"], "STEP 2 — Rosetta ddG scanning (construct v2)")

    # ── STEP 3: ddG filter ───────────────────────────────────────────────────
    run(["analysis/igem_filter_v2.py"], "STEP 3 — Stability filter (ddG < 5 REU)")

    # ── STEP 4: Conservation ─────────────────────────────────────────────────
    run(["analysis/conservation_msa_score_v2.py"], "STEP 4 — MSA conservation scoring")

    # ── STEP 5: Accessibility ────────────────────────────────────────────────
    run(["analysis/accessibility_score_v2.py"], "STEP 5 — Solvent accessibility")

    # ── STEP 6: Epitope ──────────────────────────────────────────────────────
    run(["analysis/epitope_prediction_score_v2.py"], "STEP 6 — Epitope preservation scoring")

    # ── STEP 7: Structural context ───────────────────────────────────────────
    run(["analysis/structural_context_score_v2.py"], "STEP 7 — Structural context scoring")

    # ── STEP 8: Expression ───────────────────────────────────────────────────
    run(["analysis/expression_score_v2.py"], "STEP 8 — Expression scoring")

    # ── STEP 9: Solubility ───────────────────────────────────────────────────
    run(["analysis/solubility_score_v2.py"], "STEP 9 — Solubility scoring")

    # ── STEP 10: Final ranking ───────────────────────────────────────────────
    run(["analysis/combine_scores_v2.py"], "STEP 10 — Final ranking")

    print("\n" + "="*60)
    print("  v2 pipeline complete.")
    print(f"  Results: {ROOT}/construct_v2/results/")
    print(f"  Top candidates: {ROOT}/construct_v2/results/final_ranked_simple.csv")
    print("="*60)


if __name__ == "__main__":
    main()
