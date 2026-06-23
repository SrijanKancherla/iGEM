import pandas as pd
import numpy as np

from analysis.scoring_utils import (
    tcell_penalty,
    bcell_penalty,
    conservation_penalty,
    expression_penalty
)

# ============================================================================
# LOAD DATA
# ============================================================================

df = pd.read_csv("results/3_epitope_scores.csv")

# PHASE 1: Add solubility scores (NEW)
solubility_df = pd.read_csv("results/solubility_scores.csv")
df = df.merge(
    solubility_df[["pos", "wt", "mut", "solubility_score"]],
    on=["pos", "wt", "mut"],
    how="left"
)

# ============================================================================
# COMPUTE PENALTIES
# ============================================================================

df["tcell_penalty"] = df["pos"].apply(tcell_penalty)
df["bcell_penalty"] = df["pos"].apply(bcell_penalty)
df["conservation_penalty"] = df["pos"].apply(conservation_penalty)
df["expression_penalty"] = df.apply(lambda r: expression_penalty(r["wt"], r["mut"]), axis=1)

# ============================================================================
# NORMALIZE SCORES TO [0, 1] SCALE
# ============================================================================
# This makes combination more interpretable and principled.

# ddG: lower (more stable) is better. Normalize to [0, 1] where 1 = best stability
# Assume ddG < -5 is excellent, ddG > 5 is poor
ddg_min = df["ddg"].quantile(0.1)  # 10th percentile is "good"
ddg_max = df["ddg"].quantile(0.9)  # 90th percentile is "bad"
df["ddg_norm"] = 1.0 - np.clip((df["ddg"] - ddg_min) / (ddg_max - ddg_min), 0, 1)

# Solubility: already on 0-10 scale, normalize to [0, 1]
df["solubility_norm"] = df["solubility_score"] / 10.0

# Penalties: convert to [0, 1] where 1 = no penalty, 0 = worst penalty
# Each penalty has different scale; normalize them
tcell_max = df["tcell_penalty"].max()
bcell_max = df["bcell_penalty"].max()
cons_max = df["conservation_penalty"].max()
expr_max = df["expression_penalty"].max()

df["tcell_norm"] = 1.0 - (df["tcell_penalty"] / max(tcell_max, 1))
df["bcell_norm"] = 1.0 - (df["bcell_penalty"] / max(bcell_max, 1))
df["conservation_norm"] = 1.0 - (df["conservation_penalty"] / max(cons_max, 1))
df["expression_norm"] = 1.0 - (df["expression_penalty"] / max(expr_max, 1))

# ============================================================================
# FINAL SCORE WITH PRINCIPLED WEIGHTING
# ============================================================================
# Weights based on design priorities for E. coli thermostable protein:
# - Thermostability (main goal): 35%
# - Solubility (critical for E. coli): 25%
# - Immune epitope preservation: 25%
# - Other expression factors: 15%

df["final_score"] = (
    0.35 * df["ddg_norm"]           # Thermostability
    + 0.25 * df["solubility_norm"]  # Solubility (NEW)
    + 0.15 * df["tcell_norm"]       # T-cell epitope protection
    + 0.15 * df["bcell_norm"]       # B-cell epitope protection
    + 0.10 * df["conservation_norm"] # Sequence conservation
    + 0.05 * df["expression_norm"]  # Expression penalties
)

# Ensure score is in [0, 1]
df["final_score"] = np.clip(df["final_score"], 0, 1)

# Sort by final score (descending = best first)
df = df.sort_values("final_score", ascending=False)

# Save comprehensive results
df.to_csv("results/final_ranked.csv", index=False)

# Save simplified results (for easier reading)
output_cols = [
    "pos", "wt", "mut", "ddg", "solubility_score",
    "tcell_penalty", "bcell_penalty", "conservation_penalty",
    "expression_penalty", "final_score"
]
df[output_cols].to_csv("results/final_ranked_simple.csv", index=False)

print("\n" + "="*70)
print("PHASE 1 RESULTS: Thermostability + Solubility Pipeline")
print("="*70)
print(f"\nTotal mutations analyzed: {len(df)}")
print(f"\nFinal score range: [{df['final_score'].min():.3f}, {df['final_score'].max():.3f}]")
print(f"Mean final score: {df['final_score'].mean():.3f}")

print("\n" + "-"*70)
print("TOP 20 MUTATIONS (Ranked by combined score)")
print("-"*70)
print(df[output_cols].head(20).to_string(index=False))

print("\n" + "-"*70)
print("WEIGHT CONTRIBUTION ANALYSIS")
print("-"*70)
print(f"Thermostability:  35% (ddG-based)")
print(f"Solubility:       25% (E. coli expression)")
print(f"T-cell epitopes:  15% (immune escape)")
print(f"B-cell epitopes:  15% (immune escape)")
print(f"Conservation:     10% (sequence stability)")
print(f"Expression:        5% (E. coli toxicity)")
print()

print("Saved results to:")
print("  - results/final_ranked.csv (comprehensive)")
print("  - results/final_ranked_simple.csv (simplified)")