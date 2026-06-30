"""Expression scoring for the v2 construct-based pipeline."""

import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from scoring_utils import expression_penalty

INPUT_CSV  = Path("construct_v2/results/5_structural_context.csv")
OUTPUT_CSV = Path("construct_v2/results/6_expression.csv")


def main():
    df = pd.read_csv(INPUT_CSV)
    df["expression_penalty"] = df.apply(
        lambda r: expression_penalty(r["wt"], r["mut"]), axis=1
    )
    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Expression scoring complete → {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
