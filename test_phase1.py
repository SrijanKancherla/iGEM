#!/usr/bin/env python3
"""
Quick test script for Phase 1 components.

This validates ESM2 and solubility scoring without running expensive Rosetta calculations.
Useful for checking setup before running the full pipeline.

Usage:
    python test_phase1.py
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np


def test_esm2_import():
    """Test ESM2 availability."""
    print("Test 1: ESM2 Import")
    print("-" * 50)

    try:
        import esm
        print("✓ ESM2 package available")
        print(f"  ESM version: {esm.__version__ if hasattr(esm, '__version__') else 'unknown'}")
        return True
    except ImportError as e:
        print(f"✗ ESM2 import failed: {e}")
        print("\n  Fix: pip install fair-esm[esmfold]")
        return False


def test_biopython():
    """Test Biopython availability."""
    print("\nTest 2: Biopython Import")
    print("-" * 50)

    try:
        from Bio import SeqIO, SeqUtils
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        print("✓ Biopython available")
        return True
    except ImportError as e:
        print(f"✗ Biopython import failed: {e}")
        print("\n  Fix: pip install biopython")
        return False


def test_solubility_functions():
    """Test solubility scoring functions."""
    print("\nTest 3: Solubility Functions")
    print("-" * 50)

    try:
        from analysis.solubility_score import (
            calculate_gravy,
            calculate_cai,
            score_solubility,
        )

        # Test GRAVY on wild-type
        wt_seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"
        gravy = calculate_gravy(wt_seq)
        print(f"✓ GRAVY calculation works: {gravy:.3f}")

        # Test solubility scoring
        mut_seq = wt_seq.replace("M", "L")  # Simple mutation
        score, components = score_solubility(wt_seq, mut_seq, 0, "M", "L")
        print(f"✓ Solubility scoring works: {score:.2f}/10")
        print(f"  Components: GRAVY={components['gravy_score']:.1f}, "
              f"Charge={components['charge_score']:.1f}, "
              f"Cys={components['cys_score']:.1f}")
        return True

    except Exception as e:
        print(f"✗ Solubility test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_input_files():
    """Test that input files exist."""
    print("\nTest 4: Input Files")
    print("-" * 50)

    required_files = [
        "data/cleaned/hbsag.fasta",
        "data/cleaned/hbsag.pdb",
    ]

    all_exist = True
    for file in required_files:
        if Path(file).exists():
            size = Path(file).stat().st_size
            print(f"✓ {file} ({size} bytes)")
        else:
            print(f"✗ {file} NOT FOUND")
            all_exist = False

    return all_exist


def test_output_directories():
    """Test that output directories can be created."""
    print("\nTest 5: Output Directories")
    print("-" * 50)

    dirs = [
        "results",
        "data/mutations",
    ]

    all_ok = True
    for dir in dirs:
        try:
            Path(dir).mkdir(parents=True, exist_ok=True)
            print(f"✓ {dir} is writable")
        except Exception as e:
            print(f"✗ Cannot create {dir}: {e}")
            all_ok = False

    return all_ok


def test_pandas_workflow():
    """Test basic pandas workflow for pipeline."""
    print("\nTest 6: Pandas Workflow")
    print("-" * 50)

    try:
        # Create a mock mutations dataframe
        mutations = [
            {"pos": 0, "wt": "M", "mut": "L", "esm_score": -2.5, "solubility_score": 7.5},
            {"pos": 10, "wt": "K", "mut": "R", "esm_score": -1.2, "solubility_score": 8.0},
            {"pos": 20, "wt": "C", "mut": "S", "esm_score": -3.1, "solubility_score": 6.0},
        ]

        df = pd.DataFrame(mutations)

        # Test normalization
        df["esm_norm"] = (df["esm_score"] - df["esm_score"].min()) / (df["esm_score"].max() - df["esm_score"].min())
        df["final_score"] = 0.5 * df["esm_norm"] + 0.5 * (df["solubility_score"] / 10)

        print(f"✓ Workflow works with {len(df)} mutations")
        print(f"  Final score range: {df['final_score'].min():.3f} - {df['final_score'].max():.3f}")
        print(f"\n  Sample ranked mutation:")
        print(df[["pos", "wt", "mut", "final_score"]].to_string(index=False))

        return True

    except Exception as e:
        print(f"✗ Pandas workflow failed: {e}")
        return False


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("PHASE 1 COMPONENT TEST SUITE")
    print("="*70 + "\n")

    results = {
        "ESM2 Import": test_esm2_import(),
        "Biopython": test_biopython(),
        "Solubility Functions": test_solubility_functions(),
        "Input Files": test_input_files(),
        "Output Directories": test_output_directories(),
        "Pandas Workflow": test_pandas_workflow(),
    }

    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)

    passed = sum(1 for v in results.values() if v)
    total = len(results)

    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status:8} {test_name}")

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n✓ All tests passed! Ready to run Phase 1 pipeline.")
        print("\nNext step:")
        print("  python run_phase1_pipeline.py --device cpu")
        return 0
    else:
        print("\n✗ Some tests failed. See output above for fixes.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
