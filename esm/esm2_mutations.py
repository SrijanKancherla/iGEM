"""
ESM2-based mutation generation with likelihood scoring.

Replaces heuristic mutation generation with actual pre-trained language model scores.
Uses ESM2-t33 for computational efficiency while maintaining high quality predictions.

Installation requirement:
    pip install fair-esm
"""

import pandas as pd
import numpy as np
import torch
from pathlib import Path

# Try to import ESM; provide helpful error if not available
try:
    import esm
except ImportError:
    raise ImportError(
        "ESM package not found. Install with:\n"
        "  pip install fair-esm[esmfold]\n"
        "or for CPU-only:\n"
        "  pip install fair-esm"
    )


def load_esm2_model(model_name="esm2_t33_650M_UR50D", device="cpu"):
    """
    Load ESM2 model for mutation scoring.

    Args:
        model_name: ESM2 model variant (trade-off between speed and accuracy)
            - esm2_t6_8M_UR50D (fast, lower accuracy)
            - esm2_t12_35M_UR50D (medium)
            - esm2_t33_650M_UR50D (recommended, good balance)
            - esm2_t36_3B_UR50D (slowest, highest accuracy)
        device: "cpu" or "cuda" (GPU required for large models)

    Returns:
        model, alphabet
    """
    print(f"Loading {model_name}...")
    if hasattr(esm.pretrained, model_name):
        model, alphabet = getattr(esm.pretrained, model_name)()
    else:
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(model_name)
    model = model.to(device)
    model.eval()
    print(f"Model loaded on {device}")
    return model, alphabet


def score_mutation(
    model,
    alphabet,
    sequence,
    position,
    wt_aa,
    mut_aa,
    device="cpu"
):
    """
    Score a single mutation using ESM2 masked prediction.

    Uses the masked token prediction approach: mask the position and compute
    log-probability of the mutant amino acid.

    Args:
        model: ESM2 model
        alphabet: ESM2 alphabet
        sequence: Full protein sequence
        position: 0-indexed position in sequence
        wt_aa: Wild-type amino acid
        mut_aa: Mutant amino acid
        device: Computation device

    Returns:
        log_prob: Log-probability of the mutation (higher = more likely)
    """
    # Tokenize with batch_converter so BOS/EOS are included
    batch_converter = alphabet.get_batch_converter()
    _, _, tokens = batch_converter([("seq", sequence)])
    tokens = tokens.to(device)

    # Mask the target position in token space (+1 for BOS at index 0)
    tokens_masked = tokens.clone()
    tokens_masked[0, position + 1] = alphabet.mask_idx

    # Get logits
    with torch.no_grad():
        results = model(tokens_masked, repr_layers=[])
        logits = results["logits"][0]

    # Extract logits at the masked position (+1 for BOS token)
    pos_logits = logits[position + 1]

    # Get log probabilities
    log_probs = torch.log_softmax(pos_logits, dim=0)

    # Look up the mutation amino acid index
    mut_idx = alphabet.get_idx(mut_aa)

    if mut_idx is None:
        return float("-inf")

    log_prob = float(log_probs[mut_idx].cpu())

    return log_prob


def generate_mutations_esm2(
    fasta_file,
    output_csv,
    top_k_per_position=10,
    ddg_cutoff=5.0,
    model_name="esm2_t33_650M_UR50D",
    device="cpu",
    batch_size=1,
):
    """
    Generate candidate mutations ranked by ESM2 likelihood.

    For each position, scores all 19 amino acid substitutions and ranks them.
    Can optionally filter to top-k per position to reduce computational load.

    Args:
        fasta_file: Path to input FASTA file
        output_csv: Path to output CSV with ESM scores
        top_k_per_position: Only keep top-k mutations per position (None = all)
        ddg_cutoff: Rough filter (keep mutations predicted to be < ddg_cutoff kcal/mol)
        model_name: ESM2 model variant
        device: "cpu" or "cuda"
        batch_size: Batch size for processing (1 recommended for memory efficiency)
    """
    # Load model
    model, alphabet = load_esm2_model(model_name, device)

    # Read sequence
    with open(fasta_file, "r") as f:
        lines = f.readlines()

    # Simple FASTA parsing (assumes single sequence)
    sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
    sequence = sequence.upper()

    print(f"Sequence length: {len(sequence)}")
    print(f"Scoring mutations with ESM2...")

    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    mutations = []

    # Score all possible mutations
    for pos in range(len(sequence)):
        wt = sequence[pos]
        pos_scores = []

        for mut in amino_acids:
            if mut == wt:
                continue

            # Score mutation
            log_prob = score_mutation(
                model, alphabet, sequence, pos, wt, mut, device
            )

            pos_scores.append({
                "pos": pos,
                "wt": wt,
                "mut": mut,
                "esm_score": log_prob
            })

        # Sort by ESM score and optionally filter
        pos_scores.sort(key=lambda x: x["esm_score"], reverse=True)

        if top_k_per_position is not None:
            pos_scores = pos_scores[:top_k_per_position]

        mutations.extend(pos_scores)

        if (pos + 1) % 50 == 0:
            print(f"  Processed {pos + 1}/{len(sequence)} positions")

    # Convert to DataFrame
    df = pd.DataFrame(mutations)
    df = df.sort_values("esm_score", ascending=False)

    # Save
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} mutations to {output_csv}")
    print(f"ESM score range: [{df['esm_score'].min():.2f}, {df['esm_score'].max():.2f}]")

    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate ESM2-ranked mutation candidates.")
    parser.add_argument("--fasta", default="data/cleaned/hbsag.fasta")
    parser.add_argument("--output", default="data/mutations/esm2_mutations.csv")
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    parser.add_argument("--model-name", default="esm2_t33_650M_UR50D")
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    df = generate_mutations_esm2(
        args.fasta,
        args.output,
        top_k_per_position=args.top_k,
        device=args.device,
        model_name=args.model_name,
    )

    print("\nTop 20 mutations by ESM2 score:")
    print(df.head(20))
