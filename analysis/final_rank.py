import pandas as pd

from scoring_utils import (
    tcell_penalty,
    bcell_penalty
)

df = pd.read_csv(
    "results/4_expression.csv"
)

df["tcell_penalty"] = (
    df["pos"]
    .apply(tcell_penalty)
)

df["bcell_penalty"] = (
    df["pos"]
    .apply(bcell_penalty)
)

df["final_score"] = (
    df["ddg"]
    + df["conservation_penalty"]
    + df["expression_penalty"]
    + df["tcell_penalty"]
    + df["bcell_penalty"]
)

df = df.sort_values(
    "final_score"
)

df.to_csv(
    "results/final_ranked.csv",
    index=False
)

print(
    df.head(20)
)