# HBsAg Thermostability Engineering Pipeline

Computational design of thermostable Hepatitis B surface antigen mutations while preserving immune function.

## Quick Start

```bash
# Install
pip install -r requirements.txt

# Run (Phase 1: ESM2 + Rosetta + Solubility)
python run_phase1_pipeline.py --device cpu

# Or with GPU (much faster)
python run_phase1_pipeline.py --device cuda
```

## Results

Top candidates ranked by combined score in `results/final_ranked.csv`:
- **Top 5:** A76K, A76D, G159L, G159A, A90S
- **Predicted Tm improvement:** +1.5 to +2.5°C per mutation
- **E. coli solubility:** Safe (score ≥ 9.6/10)

## What It Does

1. **Mutation generation** (ESM2) — Ranks mutations by evolutionary plausibility
2. **Thermostability** (Rosetta) — Calculates ΔΔG for stability improvement
3. **Solubility** (E. coli) — Ensures expression won't fail from aggregation
4. **Epitope preservation** — Maintains immune recognition
5. **Ranking** — Integrates all metrics into final score

## Validation

```bash
python test_phase1.py  # Quick validation
```

## For Detailed Information

See `archive/results_analysis/` folder for:
- `METHODOLOGY_AND_SELECTION_REPORT.md` — Full technical documentation
- `PROFESSOR_REVIEW_GUIDE.md` — For academic review
- `final_ranked.csv` — All 2,931 ranked mutations

## Requirements

- Python 3.8+
- PyRosetta (for Rosetta calculations)
- ESM2 (for mutation scoring)
- Biopython, pandas, numpy

See `requirements.txt` for full list.

---

**Status:** Phase 1 complete ✓  
**For questions:** See archive/ folder documentation
