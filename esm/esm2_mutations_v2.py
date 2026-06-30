"""
ESM2 mutation scoring for the v2 construct-based pipeline.

Reference sequence: 226 aa HBsAg (construct without His-tag).
  construct_v2/data/cleaned/hbsag_construct.fasta

Numbering:
  pos          - 0-indexed position in the 226 aa sequence
  s_protein_pos = pos + 1   (standard HBsAg S-protein position)
  construct_pos = pos + 13  (1-indexed in the full 238 aa His-tagged construct)

Mutations at the old A76/A90 positions will now correctly show WT=C.
The His-tag (construct positions 1-12) is excluded from scoring.
"""

import pandas as pd
import numpy as np
import torch
from pathlib import Path

try:
    import esm
except ImportError:
    raise ImportError("Install ESM: pip install fair-esm")


FASTA_V2    = Path("construct_v2/data/cleaned/hbsag_construct.fasta")
OUTPUT_V2   = Path("construct_v2/mutations/esm2_mutations.csv")
HIS_TAG_LEN = 12


def load_esm2_model(model_name="esm2_t33_650M_UR50D", device="cpu"):
    print(f"Loading {model_name}...")
    if hasattr(esm.pretrained, model_name):
        model, alphabet = getattr(esm.pretrained, model_name)()
    else:
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(model_name)
    model = model.to(device)
    model.eval()
    print(f"Model loaded on {device}")
    return model, alphabet


def score_mutation(model, alphabet, sequence, position, mut_aa, device="cpu"):
    batch_converter = alphabet.get_batch_converter()
    _, _, tokens = batch_converter([("seq", sequence)])
    tokens = tokens.to(device)
    tokens_masked = tokens.clone()
    tokens_masked[0, position + 1] = alphabet.mask_idx  # +1 for BOS token
    with torch.no_grad():
        logits = model(tokens_masked, repr_layers=[])["logits"][0]
    log_probs = torch.log_softmax(logits[position + 1], dim=0)
    mut_idx = alphabet.get_idx(mut_aa)
    if mut_idx is None:
        return float("-inf")
    return float(log_probs[mut_idx].cpu())


def generate_mutations_v2(
    fasta_file=FASTA_V2,
    output_csv=OUTPUT_V2,
    top_k_per_position=10,
    model_name="esm2_t33_650M_UR50D",
    device="cpu",
):
    model, alphabet = load_esm2_model(model_name, device)

    with open(fasta_file) as f:
        lines = f.readlines()
    sequence = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
    n = len(sequence)
    print(f"Sequence length: {n} aa")
    print("Scoring all single-residue substitutions with ESM2...")

    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    mutations = []

    for pos in range(n):
        wt = sequence[pos]
        pos_scores = []
        for mut in amino_acids:
            if mut == wt:
                continue
            log_prob = score_mutation(model, alphabet, sequence, pos, mut, device)
            pos_scores.append({
                "pos":           pos,
                "s_protein_pos": pos + 1,
                "construct_pos": pos + HIS_TAG_LEN + 1,
                "wt":            wt,
                "mut":           mut,
                "esm_score":     log_prob,
            })
        pos_scores.sort(key=lambda x: x["esm_score"], reverse=True)
        if top_k_per_position is not None:
            pos_scores = pos_scores[:top_k_per_position]
        mutations.extend(pos_scores)
        if (pos + 1) % 50 == 0:
            print(f"  {pos + 1}/{n} positions processed")

    df = pd.DataFrame(mutations).sort_values("esm_score", ascending=False)
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"\nSaved {len(df)} scored mutations → {output_csv}")
    print(f"ESM score range: [{df['esm_score'].min():.2f}, {df['esm_score'].max():.2f}]")
    print("\nTop 20 mutations:")
    print(df[["s_protein_pos","wt","mut","esm_score"]].head(20).to_string(index=False))
    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="ESM2 mutation scoring for construct v2.")
    parser.add_argument("--fasta",      default=str(FASTA_V2))
    parser.add_argument("--output",     default=str(OUTPUT_V2))
    parser.add_argument("--top-k",      type=int, default=10)
    parser.add_argument("--device",     default="cpu", choices=["cpu", "cuda"])
    parser.add_argument("--model-name", default="esm2_t33_650M_UR50D")
    args = parser.parse_args()

    generate_mutations_v2(
        fasta_file=args.fasta,
        output_csv=args.output,
        top_k_per_position=args.top_k,
        model_name=args.model_name,
        device=args.device,
    )
