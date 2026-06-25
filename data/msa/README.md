# MSA Input

Place an aligned HBsAg FASTA file here:

```text
data/msa/hbsag_msa.fasta
```

The first sequence should match `data/cleaned/hbsag.fasta` after removing gaps.

If this file is missing, Phase 1 still runs, but conservation scoring falls back to the built-in placeholder high-risk residue rules and marks rows with:

```text
conservation_mode = fallback_no_msa
```

