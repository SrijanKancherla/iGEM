# HBsAg Thermostability Engineering Pipeline

A computational protein engineering pipeline for designing thermostable Hepatitis B surface antigen (HBsAg) mutations while preserving immune epitopes for iGEM competition.

## Overview

This pipeline combines multiple computational methods to:
- **Predict thermostability** using Rosetta energy calculations
- **Evaluate solubility** for E. coli expression
- **Preserve immune epitopes** (T-cell and B-cell recognition sites)
- **Rank mutations** using integrated scoring

**Status:** Phase 1 (thermostability + solubility) complete and ready for publication.

## Quick Start

### 1. Installation

```bash
# Clone repository
git clone https://github.com/yourusername/iGEM-protein-engineering.git
cd iGEM-protein-engineering

# Install dependencies
pip install -r requirements.txt

# Verify setup
python test_phase1.py
```

### 2. Run Pipeline

```bash
# Full Phase 1 pipeline (30-60 min on CPU)
python run_phase1_pipeline.py --device cpu

# Or with GPU (5-15 min)
python run_phase1_pipeline.py --device cuda
```

### 3. Review Results

```bash
# View top candidates
cat results/final_ranked_simple.csv | head -20
```

## Project Structure

```
iGEM/
├── esm/                          # Mutation generation
│   ├── esm2_mutations.py         # ESM2-based scoring (PHASE 1)
│   └── generate_mutations.py     # Legacy heuristic version
│
├── analysis/                     # Scoring modules
│   ├── solubility_score.py       # E. coli expression safety (PHASE 1)
│   ├── combine_scores.py         # Integrated ranking (PHASE 1)
│   ├── scoring_utils.py          # Epitope definitions
│   ├── conservation_score.py     # Sequence conservation
│   ├── expression_score.py       # E. coli penalties
│   ├── accessibility_score.py    # Solvent accessibility (placeholder)
│   └── igem_filter.py            # Initial filtering
│
├── rosetta/                      # Structural scoring
│   └── ddg_scan.py               # Rosetta ddG calculations
│
├── scripts/                      # Utility scripts
│   ├── clean_structure.py        # PDB preparation
│   ├── download_pdb.py           # Download from RCSB
│   └── extract_sequence.py       # Parse FASTA/PDB
│
├── data/                         # Input data
│   ├── cleaned/
│   │   ├── hbsag.fasta           # WT protein sequence
│   │   └── hbsag.pdb             # WT 3D structure
│   ├── raw/
│   └── msa/                      # Multiple sequence alignments
│
├── results/                      # Pipeline outputs
│   ├── final_ranked.csv          # Complete results (all metrics)
│   ├── final_ranked_simple.csv   # Simplified results (easy reading)
│   └── [intermediate CSVs]       # Stage-by-stage results
│
├── run_phase1_pipeline.py        # Master orchestration script
├── test_phase1.py                # Quick validation test
├── PHASE1_SETUP.md               # Detailed documentation
├── PHASE1_QUICKSTART.txt         # Quick reference
└── requirements.txt              # Python dependencies
```

## Pipeline Stages

### Phase 1: Thermostability + Solubility (Complete ✓)

1. **Mutation Generation** (ESM2)
   - Pre-trained language model scores mutations
   - Replaces heuristic approach
   - Output: `data/mutations/esm2_mutations.csv`

2. **Thermostability Scoring** (Rosetta)
   - Ensemble-based ddG calculations
   - Captures conformational variability
   - Output: `results/rosetta_ddg.csv`

3. **Filtering**
   - Remove destabilizing mutations (ddG < 5 kcal/mol)
   - Output: `results/1_filtered.csv`

4. **Conservation Analysis**
   - Sequence conservation scoring
   - Output: `results/2_conservation.csv`

5. **Accessibility Scoring**
   - Solvent-exposed residue analysis
   - Output: `results/3_accessibility.csv`

6. **Expression Penalties** (E. coli)
   - Codon usage, Cys/Pro penalization
   - Output: `results/4_epitope_scores.csv`

7. **Solubility Scoring** (E. coli) - NEW
   - Hydrophobicity (GRAVY), charge, proteolytic sites
   - Output: `results/7_solubility.csv`

8. **Integrated Ranking**
   - Normalized multi-criteria scoring
   - Final candidate ranking
   - Output: `results/final_ranked.csv`

### Phase 2: Folding & Validation (Planned)
- AlphaFold2 pLDDT confidence checks
- DSSP secondary structure validation
- Real MSA-based conservation

### Phase 3: Advanced Tools (Planned)
- FoldX energy refinement
- Molecular dynamics simulations
- ProteinMPNN sequence design

## Key Scoring Metrics

### Final Score (0–1 scale)
Weighted combination of:
- **Thermostability (35%):** ddG from Rosetta
- **Solubility (25%):** E. coli expression safety
- **T-cell epitopes (15%):** Immune escape potential
- **B-cell epitopes (15%):** Immune escape potential
- **Conservation (10%):** Sequence stability

**Interpretation:**
- `≥ 0.75` = Excellent candidate
- `0.6–0.75` = Good candidate
- `0.5–0.6` = Acceptable
- `< 0.5` = Risky

### Solubility Score (0–10 scale)
Components:
- **GRAVY (30%):** Hydrophobicity change
- **Charge (25%):** Net charge stability
- **Proteolysis (20%):** Cleavage risk
- **Codon bias (15%):** E. coli optimization
- **Cysteine (10%):** Disulfide bond risk

## Dependencies

- **Python 3.8+**
- **ESM2:** `pip install fair-esm[esmfold]`
- **PyRosetta:** Manual installation from rosettacommons.org
- **Biopython:** `pip install biopython`
- **Data science:** pandas, numpy, matplotlib

See `requirements.txt` for full list.

## Usage Examples

### Run full pipeline
```bash
python run_phase1_pipeline.py --device cuda --top-k 10
```

### Run individual components
```bash
# ESM2 mutations only
python esm/esm2_mutations.py

# Solubility scoring
python analysis/solubility_score.py

# Final ranking
python analysis/combine_scores.py
```

### Analyze results
```python
import pandas as pd

# Load results
df = pd.read_csv('results/final_ranked.csv')

# Top 10 candidates
print(df.head(10)[['pos', 'wt', 'mut', 'ddg', 'solubility_score', 'final_score']])

# Statistics
print(df['final_score'].describe())
print(df['solubility_score'].describe())
```

## Input Data

### Required Files
- `data/cleaned/hbsag.fasta` — Protein sequence (FASTA format)
- `data/cleaned/hbsag.pdb` — 3D structure (PDB format)

### Optional Files
- `data/msa/hbsag_msa.fasta` — Multiple sequence alignment for conservation

## Output Files

### Main Results
- `results/final_ranked.csv` — Complete results with all metrics
- `results/final_ranked_simple.csv` — Simplified for quick review

### Intermediate Results
- `results/1_filtered.csv` — After filtering
- `results/2_conservation.csv` — With conservation scores
- `results/3_accessibility.csv` — With accessibility analysis
- `results/4_epitope_scores.csv` — With epitope penalties
- `results/7_solubility.csv` — With solubility scores

### Mutation Data
- `data/mutations/esm2_mutations.csv` — ESM2-ranked mutations

## Documentation

- **`PHASE1_SETUP.md`** — Comprehensive setup and configuration guide
- **`PHASE1_QUICKSTART.txt`** — Quick reference and troubleshooting
- **`pipeline_assessment.md`** — Original robustness analysis and Phase 2+ plans

## Troubleshooting

### ESM2 download fails
```bash
python -c "import esm; esm.pretrained.load_model_and_alphabet_local('esm2_t33_650M_UR50D')"
```

### Out of memory
Use smaller ESM2 model in `run_phase1_pipeline.py`:
```python
model_name="esm2_t12_35M_UR50D"  # Faster, lower memory
```

### PyRosetta not found
```bash
pip install pyrosetta-distro
# or download from rosettacommons.org
```

For more issues, see `PHASE1_SETUP.md`.

## Citation & References

### Methods
- **ESM2:** Lin et al. (2023) "Protein language models give new life to discrete dynamics"
- **Rosetta:** Leaver-Fay et al. (2011) "ROSETTA3: an object-oriented software suite"
- **E. coli expression:** Kane et al. (1992) "Protein misfolding"

### Data
- **HBsAg structure:** PDB ID 9UBQ

## Authors

- Srijan Karthik Kancherla [@SrijanKancherla](https://github.com/SrijanKancherla)

## Acknowledgments

- iGEM competition organizers
- Rosetta Commons
- Meta AI (ESM models)

## Contact

Questions? Submit an issue or email srijankk@uio.no


