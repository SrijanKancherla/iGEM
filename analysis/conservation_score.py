import pandas as pd

from scoring_utils import conservation_penalty

df = pd.read_csv(
    "results/1_filtered.csv"
)

df["conservation_penalty"] = (
    df["pos"]
    .apply(conservation_penalty)
)

df.to_csv(
    "results/2_conservation.csv",
    index=False
)

print("Conservation scoring complete")