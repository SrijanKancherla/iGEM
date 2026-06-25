# Phase 2 Analysis

Phase 2 is for deeper analysis of the best Phase 1 candidates.

It does not replace the Phase 1 ranking. Instead, it prepares focused candidate sets for external modelling tools:

```text
ProteinMPNN: structure-conditioned sequence redesign around selected candidates
Molecular dynamics: final stability/flexibility checks on top candidates
```

Run the preparation script after Phase 1 has produced `results/final_ranked.csv`:

```bash
bash scripts/run_phase2_analysis.sh
```

Outputs are written under:

```text
phase2/proteinmpnn/inputs/
phase2/md/inputs/
```

These scripts prepare inputs and candidate manifests. They do not run ProteinMPNN or MD engines directly because those require separate external installations and environment-specific setup.

