#!/usr/bin/env bash

set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"
TOP_N_MPNN="${TOP_N_MPNN:-20}"
TOP_N_MD="${TOP_N_MD:-10}"

echo "========================================"
echo "Phase 2 preparation: ProteinMPNN + MD"
echo "========================================"

if [[ ! -f results/final_ranked.csv ]]; then
  echo "Missing results/final_ranked.csv. Run Phase 1 first:"
  echo "  bash scripts/run_pipeline.sh"
  exit 1
fi

echo "[1/2] Preparing ProteinMPNN candidate inputs"
"$PYTHON_BIN" phase2/proteinmpnn/prepare_inputs.py --top-n "$TOP_N_MPNN"

echo "[2/2] Preparing MD candidate inputs"
"$PYTHON_BIN" phase2/md/prepare_md_candidates.py --top-n "$TOP_N_MD"

echo ""
echo "Phase 2 preparation complete."
echo "ProteinMPNN inputs: phase2/proteinmpnn/inputs/"
echo "MD inputs: phase2/md/inputs/"
