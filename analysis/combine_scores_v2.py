"""
Final ranking for the v2 construct-based pipeline.

Mutation naming now uses the correct WT residues from the construct:
  C76K / C76D / C76N  (not A76K/D/N)
  C90S                (not A90S)
  A159X               (construct already has A at this position; these will score low)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from residue_mapping_v2 import residue_map_by_pos_v2

INPUT_CSV         = Path("construct_v2/results/7_solubility.csv")
OUTPUT_CSV        = Path("construct_v2/results/final_ranked.csv")
SIMPLE_OUTPUT_CSV = Path("construct_v2/results/final_ranked_simple.csv")


def norm_higher_better(series):
    series = series.astype(float)
    lo, hi = series.quantile(0.05), series.quantile(0.95)
    if hi == lo:
        return pd.Series(0.5, index=series.index)
    return ((series - lo) / (hi - lo)).clip(0, 1)


def norm_lower_better(series):
    return 1.0 - norm_higher_better(series)


def penalty_to_score(series):
    series = series.fillna(0).astype(float)
    max_val = max(float(series.max()), 1.0)
    return (1.0 - (series / max_val)).clip(0, 1)


def _nonempty(val):
    if val is None:
        return False
    try:
        import math
        if math.isnan(float(val)):
            return False
    except (TypeError, ValueError):
        pass
    return str(val).strip() not in ("", "nan", "NaN")


def build_flags(row):
    flags = []
    if bool(row.get("in_antibody_contact", False)):
        flags.append("protected_epitope_or_contact")
    if row.get("bonding_penalty", 0) >= 25:
        flags.append("possible_disulfide_loss")
    sf = row.get("structural_flags", "")
    if _nonempty(sf):
        flags.extend(str(sf).split(";"))
    if row.get("ddg", 0) >= 5:
        flags.append("bad_rosetta_ddg")
    if row.get("solubility_score", 10) < 5:
        flags.append("low_solubility")
    if row.get("accessibility_penalty", 0) >= 12:
        flags.append("accessibility_risk")
    if row.get("tcell_epitope_delta", 0) <= -0.2 and row.get("wt_tcell_epitope_score", 0) >= 0.55:
        flags.append("predicted_tcell_epitope_loss")
    if row.get("bcell_antigenicity_delta", 0) <= -0.15 and bool(row.get("in_mhr", False)):
        flags.append("predicted_bcell_antigenicity_loss")
    return ";".join(sorted(set(f for f in flags if f)))


def recommendation(row):
    flags = (set(str(row.get("flags","")).split(";")) if _nonempty(row.get("flags",""))
             else set())
    severe = {"possible_disulfide_loss", "buried_charge", "bad_rosetta_ddg", "low_solubility"}
    if flags & severe:
        return "reject"
    if row["final_score"] >= 0.75 and row["confidence_score"] >= 0.65 and not flags:
        return "strong_candidate"
    if row["final_score"] >= 0.60:
        return "candidate"
    return "review_manually"


def main():
    df = pd.read_csv(INPUT_CSV)
    mapping = residue_map_by_pos_v2()

    # Ensure s_protein_pos and pdb_residue columns exist
    if "s_protein_pos" not in df.columns:
        df["s_protein_pos"] = df["pos"].astype(int) + 1
    if "pdb_residue" not in df.columns:
        df["pdb_residue"] = df["s_protein_pos"]

    # Mutation label: e.g. C76K, C90S, A159V
    if "mutation" not in df.columns:
        df["mutation"] = (df["wt"]
                          + df["s_protein_pos"].astype(int).astype(str)
                          + df["mut"])

    df["stability_score"]   = norm_lower_better(df["ddg"])
    df["solubility_norm"]   = (df["solubility_score"].fillna(0) / 10.0).clip(0, 1)

    def get_penalty(col):
        return df[col] if col in df.columns else pd.Series(0.0, index=df.index)

    df["conservation_score"]        = penalty_to_score(get_penalty("conservation_penalty"))
    df["epitope_preservation_score"] = penalty_to_score(get_penalty("epitope_penalty"))
    df["accessibility_score"]        = penalty_to_score(get_penalty("accessibility_penalty"))
    df["structural_context_score"]   = penalty_to_score(get_penalty("structural_context_penalty"))
    df["expression_score"]           = penalty_to_score(get_penalty("expression_penalty"))

    df["final_score"] = (
        0.25 * df["stability_score"]
        + 0.20 * df["epitope_preservation_score"]
        + 0.15 * df["solubility_norm"]
        + 0.15 * df["conservation_score"]
        + 0.10 * df["accessibility_score"]
        + 0.10 * df["structural_context_score"]
        + 0.05 * df["expression_score"]
    ).clip(0, 1)

    method_votes = []
    for _, row in df.iterrows():
        votes = [
            row.get("esm_score", -999) > -8,
            row.get("ddg", 999) < 0,
            row.get("epitope_penalty", 0) == 0,
            row.get("tcell_epitope_delta", 0) > -0.2,
            row.get("solubility_score", 0) >= 7,
            row.get("structural_context_penalty", 0) <= 5,
        ]
        method_votes.append(sum(votes) / len(votes))
    df["confidence_score"] = method_votes

    df["flags"]          = df.apply(build_flags, axis=1)
    df["recommendation"] = df.apply(recommendation, axis=1)

    order = {"strong_candidate": 0, "candidate": 1, "review_manually": 2, "reject": 3}
    df["_rec_order"] = df["recommendation"].map(order)
    df = df.sort_values(["_rec_order", "final_score"], ascending=[True, False])
    df.insert(0, "rank", range(1, len(df) + 1))
    df = df.drop(columns=["_rec_order"])

    output_cols = [
        "rank", "pos", "s_protein_pos", "construct_pos", "wt", "mut", "mutation",
        "esm_score", "ddg", "stability_score", "rsa", "surface_class",
        "conservation_score", "epitope_penalty", "tcell_penalty", "bcell_penalty",
        "wt_tcell_epitope_score", "mut_tcell_epitope_score", "tcell_epitope_delta",
        "bcell_antigenicity_delta", "solubility_score", "expression_score",
        "structural_context_penalty", "final_score", "confidence_score",
        "flags", "recommendation",
    ]
    output_cols = [c for c in output_cols if c in df.columns]

    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    df[output_cols].to_csv(SIMPLE_OUTPUT_CSV, index=False)

    print("Final ranking complete (v2 construct-based pipeline)")
    print(f"Mutations analyzed: {len(df)}")
    print(f"Saved: {OUTPUT_CSV}")
    print(f"Saved: {SIMPLE_OUTPUT_CSV}")
    print(df[output_cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
