import pandas as pd
import os

os.makedirs("results", exist_ok=True)

df = pd.read_csv("results/rosetta_ddg.csv")

# lower ddg = more stable
top5 = df.sort_values("ddg").head(5)

top5.to_csv("results/top5_mutations.csv", index=False)

print("\nTOP 5 MUTATIONS:")
print(top5)