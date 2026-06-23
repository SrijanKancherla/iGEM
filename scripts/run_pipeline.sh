#!/bin/bash

set -e

echo "================================="
echo "HBsAg DESIGN PIPELINE"
echo "================================="

echo "[1/5] Rosetta ddG..."
python rosetta/ddg_scan.py

echo "[2/5] Structural filter..."
python analysis/igem_filter.py

echo "[3/5] Conservation..."
python analysis/conservation_score.py

echo "[4/5] Accessibility..."
python analysis/accessibility_score.py

echo "[5/5] Expression + ranking..."
python analysis/expression_score.py
python analysis/final_rank.py

echo ""
echo "DONE"
echo ""

echo "Results:"
echo "results/final_ranked.csv"