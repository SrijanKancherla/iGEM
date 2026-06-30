"""
Phase 1 filter for v2 pipeline: drop mutations with ddG >= 5 REU.
Reads from construct_v2/results/rosetta_ddg.csv.
"""

import pandas as pd
from pathlib import Path

INPUT_CSV  = Path("construct_v2/results/rosetta_ddg.csv")
OUTPUT_CSV = Path("construct_v2/results/1_filtered.csv")


def main():
    df = pd.read_csv(INPUT_CSV)
    before = len(df)
    df = df[df["ddg"] < 5]
    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Filtered {before} → {len(df)} mutations (ddG < 5 REU) → {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
