#!/usr/bin/env bash

set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"
DEVICE="${DEVICE:-cpu}"
TOP_K="${TOP_K:-10}"
SKIP_ESM="${SKIP_ESM:-0}"
SKIP_ROSETTA="${SKIP_ROSETTA:-0}"

echo "========================================"
echo "HBsAg thermostability design pipeline"
echo "========================================"

mkdir -p data/mutations results

echo "[0/10] Residue mapping"
"$PYTHON_BIN" analysis/residue_mapping.py

if [[ "$SKIP_ESM" == "1" ]]; then
  echo "[1/10] ESM2 mutation generation skipped"
else
  echo "[1/10] ESM2 mutation generation"
  "$PYTHON_BIN" esm/esm2_mutations.py \
    --device "$DEVICE" \
    --top-k "$TOP_K" \
    --output data/mutations/esm2_mutations.csv
fi

if [[ "$SKIP_ROSETTA" == "1" ]]; then
  echo "[2/10] Rosetta ddG skipped; using existing results/rosetta_ddg.csv"
else
  echo "[2/10] Rosetta ddG scan"
  "$PYTHON_BIN" rosetta/ddg_scan.py
fi

echo "[3/10] Initial stability filter"
"$PYTHON_BIN" analysis/igem_filter.py

echo "[4/10] Conservation scoring"
"$PYTHON_BIN" analysis/conservation_score.py

echo "[5/10] Solvent accessibility scoring"
"$PYTHON_BIN" analysis/accessibility_score.py

echo "[6/10] Epitope preservation scoring"
"$PYTHON_BIN" analysis/epitope_preservation_score.py

echo "[7/10] Structural context scoring"
"$PYTHON_BIN" analysis/structural_context_score.py

echo "[8/10] Expression scoring"
"$PYTHON_BIN" analysis/expression_score.py

echo "[9/10] Solubility scoring"
"$PYTHON_BIN" analysis/solubility_score.py

echo "[10/10] Final ranking"
"$PYTHON_BIN" analysis/combine_scores.py

echo ""
echo "Pipeline complete."
echo "Main output: results/final_ranked.csv"
echo "Readable output: results/final_ranked_simple.csv"
