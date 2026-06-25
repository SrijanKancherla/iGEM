# HBsAg Thermostability Design Pipeline

This repository contains a theoretical protein-modelling pipeline for ranking single amino-acid mutations in HBsAg. The goal is to identify candidate mutations that may improve thermostability while preserving immune epitopes and avoiding obvious risks to folding, solubility, bonding, solvent accessibility, and expression.

The pipeline is designed for in-silico exploration only. It is not a clinical, diagnostic, or experimental validation workflow.

## Method Summary

The pipeline scores candidate mutations through a staged workflow:

1. Generate candidate substitutions with ESM2 sequence-likelihood scoring.
2. Estimate structural stability with Rosetta ddG.
3. Remove catastrophic mutations with an initial ddG filter.
4. Add MSA-based conservation penalties, or a clearly marked fallback score if no MSA is provided.
5. Estimate solvent accessibility from the cleaned PDB structure.
6. Predict and score epitope preservation using HBsAg antigenic regions, antibody-contact-sensitive residues, internal MHC-I/MHC-II motif scoring, and B-cell antigenicity changes.
7. Score structural context risks such as cysteine loss, buried charge, proline/glycine introduction, and interface mutation.
8. Score expression penalties.
9. Score E. coli solubility risk.
10. Combine all scores into a final rank, confidence score, flags, and recommendation.

The final ranking uses hard-risk flags plus a weighted cumulative score. This avoids treating a mutation that breaks a protected epitope or creates a buried charge as good only because it scores well in another category.

## Repository Layout

```text
data/raw/9UBQ.pdb                 Input HBsAg-antibody complex structure
data/cleaned/hbsag.pdb            Cleaned chain A structure
data/cleaned/hbsag.fasta          Extracted HBsAg sequence
data/mutations/                   Candidate mutation CSVs
analysis/                         Scoring and ranking modules
esm/                              ESM2 candidate generation
rosetta/                          PyRosetta ddG scan
scripts/run_pipeline.sh           Main runnable pipeline script
scripts/run_phase2_analysis.sh    Phase 2 input-preparation script
phase2/                           ProteinMPNN and MD preparation modules
results/                          Generated intermediate and final outputs
```

## Main Modules

### `esm/esm2_mutations.py`

Generates candidate mutations using masked-token ESM2 scoring. For each sequence position, it scores the 19 possible substitutions and keeps the top candidates per position.

Output:

```text
data/mutations/esm2_mutations.csv
```

Important columns:

```text
pos, wt, mut, esm_score
```

### `rosetta/ddg_scan.py`

Uses PyRosetta to estimate mutation ddG. Lower ddG is treated as better for thermostability. The module relaxes a WT ensemble, mutates candidate residues, repacks/minimizes, and calculates mutant-minus-WT score.

Output:

```text
results/rosetta_ddg.csv
```

Important columns:

```text
pos, wt, mut, esm_score, ddg
```

### `analysis/igem_filter.py`

Applies the first stability filter. Currently removes mutations with:

```text
ddg >= 5
```

Output:

```text
results/1_filtered.csv
```

### `analysis/conservation_msa_score.py`

Adds a conservation penalty from an aligned FASTA MSA when available.

Expected optional MSA input:

```text
data/msa/hbsag_msa.fasta
```

The first MSA sequence should match `data/cleaned/hbsag.fasta` after removing gaps. The module maps ungapped sequence positions to MSA columns, computes residue frequencies, consensus residue, normalized entropy, and a mutation-specific conservation penalty.

If no MSA is present, the pipeline still runs and marks:

```text
conservation_mode = fallback_no_msa
```

Output:

```text
results/2_conservation.csv
```

Important columns:

```text
alignment_column
msa_depth
msa_observed
msa_consensus
wt_frequency
mut_frequency
position_entropy
conservation_score_raw
conservation_penalty
conservation_mode
```

### `analysis/residue_mapping.py`

Maps zero-indexed pipeline positions to PDB residue numbers. This matters because mutation CSVs use positions like `0, 1, 2`, while the PDB chain can start at a different residue number.

Output:

```text
results/residue_mapping.csv
```

Important columns:

```text
pos, chain, pdb_residue, insertion_code
```

### `analysis/accessibility_score.py`

Calculates approximate solvent accessibility with Biopython's Shrake-Rupley implementation. It classifies residues as:

```text
buried, partial, surface
```

It then penalizes risky substitutions based on structural exposure, such as buried charges or surface hydrophobic gains.

Output:

```text
results/3_accessibility.csv
```

Important columns:

```text
asa, rsa, surface_class, accessibility_penalty
```

### `analysis/epitope_prediction_score.py`

Predicts and scores whether a mutation may disrupt immune-sensitive regions:

```text
MHR: 99-169
Major epitope: 124-147
Antibody-contact-sensitive residues from 9UBQ context
Internal MHC-I anchor-motif model
Internal MHC-II core-motif model
B-cell antigenicity delta
```

This is an internal theoretical predictor, not a replacement for NetMHCpan, IEDB tools, or experimental epitope mapping. Conservative substitutions are penalized less than disruptive substitutions.

Output:

```text
results/4_epitope_scores.csv
```

Important columns:

```text
in_mhr, in_major_epitope, in_antibody_contact,
wt_mhci_score, mut_mhci_score, mhci_delta,
wt_mhcii_score, mut_mhcii_score, mhcii_delta,
wt_tcell_epitope_score, mut_tcell_epitope_score,
tcell_epitope_delta, bcell_antigenicity_delta,
tcell_penalty, bcell_penalty, epitope_penalty,
epitope_prediction_mode
```

### `analysis/structural_context_score.py`

Adds local structural-risk penalties for:

```text
possible disulfide loss
new cysteine
buried charge
buried hydrophobic loss
surface hydrophobic gain
new proline
new glycine
antibody-interface mutation
```

Output:

```text
results/5_structural_context.csv
```

Important columns:

```text
bonding_penalty, interface_penalty,
structural_context_penalty, structural_flags
```

### `analysis/expression_score.py`

Adds basic expression penalties, especially cysteine/proline-related expression or folding risks.

Output:

```text
results/6_expression.csv
```

### `analysis/solubility_score.py`

Scores E. coli solubility using sequence-level heuristics:

```text
GRAVY hydrophobicity change
net charge change
proteolytic-site risk
codon-use proxy
cysteine handling
```

Output:

```text
results/7_solubility.csv
```

Important columns:

```text
solubility_score, gravy_score, gravy_delta,
charge_score, charge_delta, proteolytic_score,
codon_score, cys_score
```

### `analysis/combine_scores.py`

Combines all module outputs into a final ranking.

Weighted score:

```text
0.25 stability_score
0.20 epitope_preservation_score
0.15 solubility_score
0.15 conservation_score
0.10 accessibility_score
0.10 structural_context_score
0.05 expression_score
```

It also adds:

```text
confidence_score
flags
recommendation
```

Recommendations are:

```text
strong_candidate
candidate
review_manually
reject
```

Outputs:

```text
results/final_ranked.csv
results/final_ranked_simple.csv
```

## Final Output Columns

The simplified final table contains the most useful columns:

```text
rank
pos
pdb_residue
wt
mut
mutation
esm_score
ddg
stability_score
rsa
surface_class
conservation_score
epitope_penalty
tcell_penalty
bcell_penalty
wt_tcell_epitope_score
mut_tcell_epitope_score
tcell_epitope_delta
bcell_antigenicity_delta
solubility_score
expression_score
structural_context_penalty
final_score
confidence_score
flags
recommendation
```

## Running the Pipeline

Install dependencies:

```bash
pip install -r requirements.txt
```

PyRosetta must be installed separately because it is not distributed through normal public PyPI.

Run the full pipeline:

```bash
bash scripts/run_pipeline.sh
```

Run with GPU ESM2:

```bash
DEVICE=cuda bash scripts/run_pipeline.sh
```

Keep fewer ESM2 candidates per position:

```bash
TOP_K=5 bash scripts/run_pipeline.sh
```

Reuse existing ESM2 candidates:

```bash
SKIP_ESM=1 bash scripts/run_pipeline.sh
```

Reuse existing Rosetta ddG results:

```bash
SKIP_ROSETTA=1 bash scripts/run_pipeline.sh
```

Use a specific Python executable:

```bash
PYTHON_BIN=.venv/bin/python bash scripts/run_pipeline.sh
```

## Phase 2 Analysis

After Phase 1 generates `results/final_ranked.csv`, prepare focused inputs for ProteinMPNN and MD:

```bash
bash scripts/run_phase2_analysis.sh
```

Useful variants:

```bash
TOP_N_MPNN=30 bash scripts/run_phase2_analysis.sh
TOP_N_MD=5 bash scripts/run_phase2_analysis.sh
PYTHON_BIN=.venv/bin/python bash scripts/run_phase2_analysis.sh
```

Phase 2 outputs:

```text
phase2/proteinmpnn/inputs/proteinmpnn_candidates.csv
phase2/proteinmpnn/inputs/mutable_positions.jsonl
phase2/proteinmpnn/inputs/run_notes.txt
phase2/md/inputs/md_candidates.csv
phase2/md/inputs/mutation_specs.tsv
phase2/md/inputs/mutant_fastas/
phase2/md/inputs/run_notes.txt
```

Phase 2 prepares inputs only. It does not run ProteinMPNN or MD engines because those require separate external installation and project-specific settings.

## Optional External Modules

These are useful additions, but they require external tools or heavier setup:

### FoldX

Add as a second stability estimate:

```text
foldx/foldx_scan.py -> results/foldx_ddg.csv
```

Use it to compare with Rosetta:

```text
Rosetta stable + FoldX stable = higher confidence
Rosetta/FoldX disagreement = manual review
```

### ProteinMPNN

Use as a structure-conditioned candidate generator. The Phase 2 script prepares candidate and mutable-position inputs:

```text
phase2/proteinmpnn/prepare_inputs.py
```

Recommended constraints:

```text
freeze epitope residues
freeze antibody-contact residues
freeze critical cysteines
freeze highly conserved residues
```

### Molecular Dynamics

Use only as final validation for top candidates, not as a bulk-screening step:

```text
phase2/md/prepare_md_candidates.py
```

Good MD checks:

```text
RMSD
RMSF
epitope-loop flexibility
gross unfolding
interface/contact stability
```

## Notes and Limitations

The pipeline is intentionally conservative. It is meant to rank theoretical candidates, not prove thermostability.

Current limitations:

```text
Conservation is MSA-based when `data/msa/hbsag_msa.fasta` is provided; otherwise it falls back to placeholder rules.
T-cell scoring uses internal motif-style prediction and should be replaced with NetMHCpan/IEDB-style predictors for serious work.
Rosetta ddG should be calibrated before being treated quantitatively.
Solubility is heuristic and should be treated as risk scoring, not prediction.
FoldX, ProteinMPNN, and MD are external/deeper analysis steps, not default Phase 1 stages.
```

The most important next scientific upgrade is replacing the internal epitope motif model with a real external MHC-binding predictor and adding a curated HBsAg MSA. The most important engineering upgrade is adding regression tests around the expected intermediate CSV schema.
