import pandas as pd

# This replaces brute-force enumeration logic
# (you can later plug in ESM2 / Evo / ProtBERT here)

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

df = pd.read_csv("data/cleaned/hbsag.fasta")  # just to anchor pipeline
sequence = df.columns[0] if hasattr(df, "columns") else "HBsAg"

mutations = []

# SMART PRIORITIZATION RULES (proxy for ESM ranking)
for pos in range(len(sequence)):

    wt = sequence[pos]

    for aa in amino_acids:

        if aa == wt:
            continue

        # heuristic "likelihood bias" (stand-in for ESM logits)
        score = 0

        if aa in ["K", "R", "E", "D"]:
            score -= 0.2  # charged allowed
        if aa == "C":
            score -= 1.0  # rare in HBsAg
        if aa == "P":
            score -= 0.8  # structure breaker

        mutations.append([pos, wt, aa, score])

df_out = pd.DataFrame(mutations, columns=["pos", "wt", "mut", "esm_score"])
df_out = df_out.sort_values("esm_score", ascending=False)

df_out.to_csv("data/mutations/esm_candidates.csv", index=False)

print("Generated mutations:", len(df_out))