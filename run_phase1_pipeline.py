#!/usr/bin/env python3
"""
PHASE 1 PIPELINE: Thermostability + Solubility Scoring

Orchestrates the complete Phase 1 workflow:
1. ESM2-based mutation generation (replacing heuristics)
2. Rosetta ddG calculations (existing)
3. Solubility scoring for E. coli (NEW)
4. Integrated ranking with normalized scoring

Usage:
    python run_phase1_pipeline.py

or with custom settings:
    python run_phase1_pipeline.py --device cuda --top-k 5

Requirements:
    - fair-esm (ESM2 models)
    - biopython
    - pyrosetta (for Rosetta calculations)
    - pandas, numpy
"""

import argparse
import sys
from pathlib import Path
import subprocess
import json
from datetime import datetime


def run_command(cmd, description):
    """
    Run a shell command and report results.

    Args:
        cmd: Command to run
        description: Human-readable description of what's running
    """
    print(f"\n{'='*70}")
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {description}")
    print(f"{'='*70}")
    print(f"Running: {' '.join(cmd)}")
    print()

    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"\n⚠️  Warning: Command exited with code {result.returncode}")
        return False

    print(f"\n✓ {description} complete")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Run Phase 1 protein engineering pipeline"
    )
    parser.add_argument(
        "--device",
        default="cpu",
        choices=["cpu", "cuda"],
        help="Device for ESM2 (cpu or cuda if GPU available)"
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=10,
        help="Keep top-K mutations per position (default: 10)"
    )
    parser.add_argument(
        "--skip-esm",
        action="store_true",
        help="Skip ESM2 mutation generation (use existing mutations)"
    )
    parser.add_argument(
        "--skip-rosetta",
        action="store_true",
        help="Skip Rosetta calculations (use existing ddG scores)"
    )
    parser.add_argument(
        "--skip-solubility",
        action="store_true",
        help="Skip solubility scoring"
    )
    args = parser.parse_args()

    # Create necessary directories
    Path("results").mkdir(exist_ok=True)
    Path("data/mutations").mkdir(parents=True, exist_ok=True)

    pipeline_steps = []

    # ========================================================================
    # STEP 1: ESM2 MUTATION GENERATION
    # ========================================================================
    if not args.skip_esm:
        print("\n" + "█"*70)
        print("PHASE 1, STEP 1: ESM2-Based Mutation Generation")
        print("█"*70)
        print("\nGenerating mutations using pre-trained ESM2 language model...")
        print(f"Device: {args.device}")
        print(f"Top-K per position: {args.top_k}")

        cmd = [
            sys.executable,
            "esm/esm2_mutations.py",
        ]

        if not run_command(cmd, "ESM2 mutation generation"):
            print("❌ ESM2 generation failed. Exiting.")
            sys.exit(1)

        pipeline_steps.append("✓ ESM2 mutations")
    else:
        print("\n⊘ Skipping ESM2 (using existing mutations)")
        pipeline_steps.append("⊘ ESM2 mutations (skipped)")

    # ========================================================================
    # STEP 2: ROSETTA DDG CALCULATION
    # ========================================================================
    if not args.skip_rosetta:
        print("\n" + "█"*70)
        print("PHASE 1, STEP 2: Rosetta ddG Calculation")
        print("█"*70)
        print("\nCalculating thermostability (ddG) using Rosetta...")

        cmd = [
            sys.executable,
            "rosetta/ddg_scan.py",
        ]

        if not run_command(cmd, "Rosetta ddG scan"):
            print("⚠️  Rosetta calculation failed (check PyRosetta installation)")
        else:
            pipeline_steps.append("✓ Rosetta ddG")
    else:
        print("\n⊘ Skipping Rosetta (using existing ddG scores)")
        pipeline_steps.append("⊘ Rosetta ddG (skipped)")

    # ========================================================================
    # STEP 3: FILTERING
    # ========================================================================
    print("\n" + "█"*70)
    print("PHASE 1, STEP 3: Initial Filtering")
    print("█"*70)

    cmd = [sys.executable, "analysis/igem_filter.py"]
    if not run_command(cmd, "Filter catastrophic mutations"):
        print("⚠️  Filtering step failed")
    else:
        pipeline_steps.append("✓ Filtering")

    # ========================================================================
    # STEP 4: CONSERVATION SCORING
    # ========================================================================
    print("\n" + "█"*70)
    print("PHASE 1, STEP 4: Conservation Analysis")
    print("█"*70)

    cmd = [sys.executable, "analysis/conservation_score.py"]
    if not run_command(cmd, "Conservation scoring"):
        print("⚠️  Conservation scoring failed")
    else:
        pipeline_steps.append("✓ Conservation")

    # ========================================================================
    # STEP 5: ACCESSIBILITY ANALYSIS
    # ========================================================================
    print("\n" + "█"*70)
    print("PHASE 1, STEP 5: Solvent Accessibility")
    print("█"*70)

    cmd = [sys.executable, "analysis/accesibility_score.py"]
    if not run_command(cmd, "Accessibility scoring"):
        print("⚠️  Accessibility scoring failed")
    else:
        pipeline_steps.append("✓ Accessibility")

    # ========================================================================
    # STEP 6: EXPRESSION SCORING
    # ========================================================================
    print("\n" + "█"*70)
    print("PHASE 1, STEP 6: Basic Expression Scoring")
    print("█"*70)

    cmd = [sys.executable, "analysis/expression_score.py"]
    if not run_command(cmd, "Expression scoring"):
        print("⚠️  Expression scoring failed")
    else:
        pipeline_steps.append("✓ Expression")

    # ========================================================================
    # STEP 7: SOLUBILITY SCORING (NEW - PHASE 1)
    # ========================================================================
    if not args.skip_solubility:
        print("\n" + "█"*70)
        print("PHASE 1, STEP 7: E. coli Solubility Scoring (NEW)")
        print("█"*70)
        print("\nScoring mutations for E. coli expression solubility...")

        cmd = [sys.executable, "analysis/solubility_score.py"]
        if not run_command(cmd, "Solubility scoring"):
            print("⚠️  Solubility scoring failed")
        else:
            pipeline_steps.append("✓ Solubility (NEW)")
    else:
        print("\n⊘ Skipping solubility scoring")
        pipeline_steps.append("⊘ Solubility (skipped)")

    # ========================================================================
    # STEP 8: COMBINED RANKING (UPDATED)
    # ========================================================================
    print("\n" + "█"*70)
    print("PHASE 1, STEP 8: Integrated Ranking")
    print("█"*70)
    print("\nCombining scores with normalized weighting...")
    print("  - Thermostability: 35%")
    print("  - Solubility:      25%")
    print("  - Immune epitopes: 30%")
    print("  - Other factors:   10%")

    cmd = [sys.executable, "analysis/combine_scores.py"]
    if not run_command(cmd, "Final ranking"):
        print("❌ Ranking failed")
        sys.exit(1)
    else:
        pipeline_steps.append("✓ Final ranking")

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "▓"*70)
    print("PHASE 1 PIPELINE COMPLETE")
    print("▓"*70)

    print("\nPipeline steps executed:")
    for step in pipeline_steps:
        print(f"  {step}")

    print("\nOutput files:")
    print("  - results/final_ranked.csv (comprehensive results)")
    print("  - results/final_ranked_simple.csv (simplified for reading)")

    print("\nNext steps:")
    print("  1. Review top 20 mutations in results/final_ranked_simple.csv")
    print("  2. (Optional) Phase 2: Add real DSSP + AlphaFold validation")
    print("  3. (Optional) Phase 3: FoldX cross-validation")
    print("  4. Synthesize and experimentally validate top candidates")

    print()


if __name__ == "__main__":
    main()
