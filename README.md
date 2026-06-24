# HBsAg Thermostability Design Pipeline

This repository contains a theoretical protein-modelling pipeline for ranking single amino-acid mutations in HBsAg. The goal is to identify candidate mutations that may improve thermostability while preserving immune epitopes and avoiding obvious risks to folding, solubility, bonding, solvent accessibility, and expression.

The pipeline is designed for in-silico exploration only. It is not a clinical, diagnostic, or experimental validation workflow.

## Method Summary

The pipeline scores candidate mutations through a staged workflow:

1. Generate candidate substitutions with ESM2 sequence-likelihood scoring.
2. Estimate structural stability with Rosetta ddG.
3. Remove catastrophic mutations with an initial ddG filter.
4. Add conservation penalties.
5. Estimate solvent accessibility from the cleaned PDB structure.
6. Score epitope preservation using HBsAg antigenic regions and antibody-contact-sensitive residues.
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

### `analysis/conservation_score.py`

Adds a conservation penalty. The current implementation uses known/high-risk HBsAg positions as a placeholder until an MSA-based conservation model is added.

Output:

```text
results/2_conservation.csv
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

### `analysis/epitope_preservation_score.py`

Scores whether a mutation falls in HBsAg immune-sensitive regions:

```text
MHR: 99-169
Major epitope: 124-147
Antibody-contact-sensitive residues from 9UBQ context
```

Conservative substitutions are penalized less than disruptive substitutions.

Output:

```text
results/4_epitope_scores.csv
```

Important columns:

```text
in_mhr, in_major_epitope, in_antibody_contact,
tcell_penalty, bcell_penalty, epitope_penalty
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

## Optional Future Modules

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

Use as a structure-conditioned candidate generator:

```text
proteinmpnn/generate_candidates.py -> data/mutations/proteinmpnn_candidates.csv
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
md/prepare_top_candidates.py
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
Conservation is still placeholder-based until an MSA is added.
T-cell scoring does not yet use real MHC-binding prediction.
Rosetta ddG should be calibrated before being treated quantitatively.
Solubility is heuristic and should be treated as risk scoring, not prediction.
FoldX, ProteinMPNN, and MD are documented extension points, not default stages.
```

The most important next scientific upgrade is MSA-based conservation plus real epitope prediction. The most important engineering upgrade is adding regression tests around the expected intermediate CSV schema.
