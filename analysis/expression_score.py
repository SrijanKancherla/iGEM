import pandas as pd

from scoring_utils import (
    expression_penalty
)

df = pd.read_csv(
    "results/5_structural_context.csv"
)

df["expression_penalty"] = df.apply(
    lambda r:
    expression_penalty(
        r["wt"],
        r["mut"]
    ),
    axis=1
)

df.to_csv(
    "results/6_expression.csv",
    index=False
)

print("Expression scoring complete")
