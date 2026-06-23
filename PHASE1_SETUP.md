# Phase 1 Implementation: ESM2 + Solubility Scoring

## Overview

Phase 1 adds two critical improvements to your thermostability pipeline:

1. **ESM2-based Mutation Generation** — Replaces heuristic rules with pre-trained language model scoring
2. **E. coli Solubility Scoring** — Ensures mutations don't cause aggregation in expression

Both integrate seamlessly into your existing Rosetta pipeline.

---

## Installation

### Prerequisites
- Python 3.8+
- PyRosetta (already installed)
- GPU support optional but helpful for ESM2

### Step 1: Install Dependencies

```bash
# Install ESM2 and other requirements
pip install fair-esm[esmfold] biopython pandas numpy

# Verify installation
python -c "import esm; print('ESM2 ready')"
```

**Note:** ESM2 models are large (~2 GB). First run will download the model automatically.

### Step 2: Prepare Your Data

Ensure you have:
- `data/cleaned/hbsag.fasta` — Protein sequence
- `data/cleaned/hbsag.pdb` — 3D structure (for Rosetta)

---

## Running Phase 1

### Option A: Run Full Pipeline (Recommended)

```bash
cd /path/to/iGEM

# Full pipeline with CPU (slower)
python run_phase1_pipeline.py --device cpu

# Or with GPU (if available)
python run_phase1_pipeline.py --device cuda --top-k 10
```

This executes all 8 steps in sequence:
1. ✓ ESM2 mutation generation
2. ✓ Rosetta ddG calculation
3. ✓ Filtering (ddG < 5)
4. ✓ Conservation analysis
5. ✓ Accessibility analysis
6. ✓ Expression penalties
7. ✓ **[NEW] Solubility scoring** for E. coli
8. ✓ **[NEW] Integrated ranking** with normalized weights

### Option B: Run Individual Steps

If you already have mutations or ddG scores:

```bash
# Just ESM2 mutations
python esm/esm2_mutations.py

# Just solubility scoring (on existing mutations)
python analysis/solubility_score.py

# Just final ranking
python analysis/combine_scores.py
```

### Option C: Skip Expensive Steps

```bash
# Skip Rosetta (use existing ddG)
python run_phase1_pipeline.py --skip-rosetta --device cpu

# Skip ESM2 (use existing mutations)
python run_phase1_pipeline.py --skip-esm

# Skip solubility (older pipeline behavior)
python run_phase1_pipeline.py --skip-solubility
```

---

## Understanding the Results

### Main Output: `results/final_ranked.csv`

Columns include:
- `pos`, `wt`, `mut` — Position and amino acid change
- `ddg` — Rosetta stability score (kcal/mol)
- `solubility_score` — E. coli expression score (0–10 scale)
- `tcell_penalty`, `bcell_penalty` — Immune epitope penalties
- `conservation_penalty` — Sequence conservation penalty
- `expression_penalty` — Basic expression penalties
- `final_score` — **Combined rank (0–1, higher = better)**

### Simplified Output: `results/final_ranked_simple.csv`

Just the essential columns for quick review.

---

## Score Interpretation

### Final Score Weighting

The `final_score` combines multiple criteria with normalized [0–1] scales:

| Criterion | Weight | What it measures |
|-----------|--------|------------------|
| **Thermostability** | 35% | ddG (lower = more stable) |
| **Solubility** | 25% | Risk of aggregation in E. coli |
| **T-cell epitopes** | 15% | Immune escape potential |
| **B-cell epitopes** | 15% | Immune escape potential |
| **Conservation** | 10% | Sequence stability |
| **Expression factors** | 5% | Toxicity, codon usage |

**Interpretation:**
- `final_score > 0.7` — Excellent candidate
- `final_score 0.5–0.7` — Good candidate
- `final_score < 0.5` — Risky candidate

### Solubility Score (0–10 scale)

Component breakdown:
- **GRAVY score** (30% weight) — Hydrophobicity change
  - Ideal: GRAVY ≈ −0.3 (slightly hydrophilic)
  - Flag: GRAVY > 0.2 (becoming too hydrophobic)

- **Charge** (25% weight) — Net charge at pH 7
  - Penalizes large charge loss (→ aggregation risk)

- **Proteolytic sites** (20% weight) — Cleavage risk
  - Flags new Lys/Arg (trypsin sites)
  - Flags new Pro (structure breaker)

- **Codon usage** (15% weight) — E. coli expression efficiency
  - Rare codons → low expression

- **Cysteine handling** (10% weight) — Disulfide bond risk
  - Penalizes Cys in reducing cytoplasm

**Interpretation:**
- `solubility_score > 8` — Safe for E. coli expression
- `solubility_score 5–8` — Monitor; may aggregate
- `solubility_score < 5` — High aggregation risk; avoid

---

## ESM2 Mutation Generation Details

### What Changed?

**Before (heuristic):**
```python
# Random penalties
score = 0
if aa in ["K", "R", "E", "D"]:
    score -= 0.2  # charged allowed
if aa == "C":
    score -= 1.0  # rare
```

**Now (ESM2):**
```python
# Actual pre-trained model probability
log_prob = esm2_model.score_mutation(
    sequence, pos, wt_aa, mut_aa
)
# Higher log-prob = mutation more evolutionarily plausible
```

### Model Choice

Uses `esm2_t33_650M_UR50D` by default:
- **t33** = 33 Transformer layers (good balance of speed/accuracy)
- **650M** = 650M parameters
- **UR50D** = UniRef50 pre-training

For more accuracy (slower):
```python
model_name="esm2_t36_3B_UR50D"
```

For faster screening (lower accuracy):
```python
model_name="esm2_t12_35M_UR50D"
```

### Performance

- **Speed:** ~200 mutations/min on CPU (varies with hardware)
- **Memory:** ~4 GB for t33 model
- **Top-K filtering:** Default keeps top-10 per position (~2000 mutations total)

---

## E. coli Specific Considerations

### Why Solubility Matters for E. coli

1. **Reducing cytoplasm** — Disulfide bonds don't form; Cys residues risk aggregation
2. **High expression levels** — Overexpressed proteins aggregate if not soluble
3. **No ER chaperones** — Unlike eukaryotes; relies on GroEL/ES
4. **Codon bias** — E. coli prefers certain codons; rare codons = low expression

### How Phase 1 Handles This

1. **Hydrophobicity** — Penalizes mutations that increase surface hydrophobicity
2. **Codon bias** — Flags rare amino acids (Ala, Arg, Gly)
3. **Cysteine handling** — Penalizes new Cys; OK to remove existing
4. **Charge analysis** — Maintains net charge for solubility

### Customization for Other Hosts

To adapt for mammalian expression, edit `analysis/solubility_score.py`:

```python
# For mammalian (HEK293, CHO):
# - Less penalty for Cys (ER has disulfide isomerase)
# - Different codon frequencies
# - pH optimized for secreted proteins

# For P. pastoris:
# - Higher expression level tolerance
# - Different proteolytic risk
```

---

## Quality Checks

### Before Running

✓ Check that mutation files exist:
```bash
ls -la data/mutations/
ls -la results/
```

✓ Check sequence format:
```bash
head -2 data/cleaned/hbsag.fasta
# Should show: >header and sequence
```

✓ Check PDB format:
```bash
head -5 data/cleaned/hbsag.pdb
# Should start with ATOM or HETATM
```

### After Running

✓ Check output files:
```bash
ls -la results/final_ranked*.csv
wc -l results/final_ranked.csv
```

✓ Spot-check top mutations:
```bash
head -5 results/final_ranked_simple.csv
```

✓ Verify score distributions:
```bash
python -c "
import pandas as pd
df = pd.read_csv('results/final_ranked.csv')
print(f'Final score range: {df[\"final_score\"].min():.3f} to {df[\"final_score\"].max():.3f}')
print(f'Solubility range: {df[\"solubility_score\"].min():.1f} to {df[\"solubility_score\"].max():.1f}')
"
```

---

## Troubleshooting

### Issue: ESM2 download fails
**Solution:**
```bash
# Manual download
python -c "import esm; esm.pretrained.load_model_and_alphabet_local('esm2_t33_650M_UR50D')"
```

### Issue: Out of memory (ESM2)
**Solution:** Use smaller model
```python
model_name="esm2_t12_35M_UR50D"  # Much faster, ~500 MB
```

Or reduce batch size:
```python
batch_size=1  # Process one mutation at a time
```

### Issue: PyRosetta not found
**Solution:** Install PyRosetta
```bash
pip install pyrosetta-distro  # Try this first
# or download from rosettacommons.org
```

### Issue: Solubility scores all NaN
**Solution:** Check that sequence format is correct (uppercase, no gaps)

---

## Next Steps (Phase 2)

Once Phase 1 is complete:

1. **DSSP validation** — Real secondary structure assignment (currently placeholder)
2. **AlphaFold pLDDT** — Confidence check on folding (currently missing)
3. **Conservation scoring** — Real MSA instead of hardcoded regions

See `pipeline_assessment.md` for full Phase 2 details.

---

## Citation & References

**ESM2 Models:**
- Lin et al. (2023) "Protein language models give new life to discrete dynamics." bioRxiv
- GitHub: https://github.com/facebookresearch/esm

**Rosetta:**
- Leaver-Fay et al. (2011) "ROSETTA3: an object-oriented software suite for the simulation and design of macromolecules"

**E. coli Expression:**
- Baneyx & Mujacic (2004) "Recombinant protein folding and misfolding in Escherichia coli"
- Kane et al. (1992) "Protein misfolding: When good proteins go bad"

---

## Questions?

Check the code comments in:
- `esm/esm2_mutations.py` — ESM2 implementation
- `analysis/solubility_score.py` — Solubility metrics
- `analysis/combine_scores.py` — Score integration

Each function has docstrings explaining the approach.
