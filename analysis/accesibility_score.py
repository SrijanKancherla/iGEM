import pandas as pd

def accessibility_score(pos):

    # placeholder until DSSP

    surface_region = range(99,170)

    return 1 if pos in surface_region else 0


df = pd.read_csv(
    "results/2_conservation.csv"
)

df["surface"] = (
    df["pos"]
    .apply(accessibility_score)
)

df.to_csv(
    "results/3_accessibility.csv",
    index=False
)

print("Accessibility scoring complete")