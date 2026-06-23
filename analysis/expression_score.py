import pandas as pd

from scoring_utils import (
    expression_penalty
)

df = pd.read_csv(
    "results/3_accessibility.csv"
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
    "results/4_expression.csv",
    index=False
)

print("Expression scoring complete")