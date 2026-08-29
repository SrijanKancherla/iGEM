import pandas as pd

df = pd.read_csv("results/rosetta_ddg.csv")

# remove catastrophic mutations

df = df[df["ddg"] < 5]

df.to_csv(
    "results/1_filtered.csv",
    index=False
)

print("Remaining mutations:", len(df))