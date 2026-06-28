"""
ProteinMPNN score analysis for HBsAg WT and mutants.

The score-only run scored each mutant sequence against the WT PDB backbone.
Lower score = -log_prob = higher model confidence = better sequence-structure fit.

FASTA order from shortlist_candidates.fasta:
  fasta_1 = A76K
  fasta_2 = A76D
  fasta_3 = A90S
  fasta_4 = A76N
  fasta_5 = A90T
  fasta_6 = G159A
  pdb     = WT (native PDB sequence)

Produces in phase2/proteinmpnn/results/:
  mpnn_scores.png       -- mean score per sequence (bar chart)
  mpnn_delta_score.png  -- ΔScore vs WT (mutant - WT); positive = worse fit
  mpnn_summary.csv      -- numeric table
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCORE_DIR = Path(__file__).parent / "results" / "scoring" / "score_only"
OUT_DIR = Path(__file__).parent / "results"
OUT_DIR.mkdir(exist_ok=True)

# Mapping: npz filename stem → label
FASTA_MAP = {
    "hbsag_pdb":    "WT",
    "hbsag_fasta_1": "A76K",
    "hbsag_fasta_2": "A76D",
    "hbsag_fasta_3": "A90S",
    "hbsag_fasta_4": "A76N",
    "hbsag_fasta_5": "A90T",
    "hbsag_fasta_6": "G159A",
}

COLORS = {
    "WT":    "#444444",
    "A76K":  "#3a7ebf",
    "A76D":  "#e05555",
    "A90S":  "#4caf50",
    "A76N":  "#9b59b6",
    "A90T":  "#f39c12",
    "G159A": "#1abc9c",
}


# ── load scores ──────────────────────────────────────────────────────────────

rows = []
for stem, label in FASTA_MAP.items():
    path = SCORE_DIR / f"{stem}.npz"
    if not path.exists():
        print(f"WARNING: {path} not found, skipping")
        continue
    d = np.load(str(path))
    # 'score' is the region-masked score (mutable positions only or all)
    # 'global_score' covers the full sequence
    score_mean = float(np.mean(d["score"]))
    global_mean = float(np.mean(d["global_score"]))
    rows.append({
        "label": label,
        "mean_score": round(score_mean, 4),
        "mean_global_score": round(global_mean, 4),
    })

df = pd.DataFrame(rows)

# Compute ΔScore vs WT (positive = worse sequence-structure fit than WT)
wt_score = df.loc[df["label"] == "WT", "mean_score"].values[0]
wt_global = df.loc[df["label"] == "WT", "mean_global_score"].values[0]
df["delta_score"] = (df["mean_score"] - wt_score).round(4)
df["delta_global_score"] = (df["mean_global_score"] - wt_global).round(4)

df.to_csv(OUT_DIR / "mpnn_summary.csv", index=False)
print("Saved: mpnn_summary.csv")
print("\n=== ProteinMPNN Score Summary ===")
print(df.to_string(index=False))


# ── Figure 1: raw mean scores ─────────────────────────────────────────────────

order = ["WT", "A76K", "A76D", "A90S", "A76N", "A90T", "G159A"]
df_plot = df.set_index("label").reindex(order).reset_index()

fig, ax = plt.subplots(figsize=(9, 5))
bars = ax.bar(
    df_plot["label"],
    df_plot["mean_global_score"],
    color=[COLORS.get(l, "#888") for l in df_plot["label"]],
    width=0.6, edgecolor="white",
)
for bar, val in zip(bars, df_plot["mean_global_score"]):
    ax.text(bar.get_x() + bar.get_width() / 2, val + 0.005,
            f"{val:.4f}", ha="center", va="bottom", fontsize=8.5)

ax.axhline(wt_global, ls="--", color="#444444", lw=1.2, alpha=0.7, label=f"WT baseline ({wt_global:.4f})")
ax.set_ylabel("Mean ProteinMPNN score (-log prob, global)", fontsize=11)
ax.set_title("ProteinMPNN sequence–structure compatibility\n(lower = better fit to backbone)", fontsize=12)
ax.set_ylim(0, df_plot["mean_global_score"].max() * 1.2)
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
fig.savefig(OUT_DIR / "mpnn_scores.png", dpi=150)
plt.close(fig)
print("Saved: mpnn_scores.png")


# ── Figure 2: ΔScore vs WT ────────────────────────────────────────────────────

df_mut = df_plot[df_plot["label"] != "WT"].copy()
colors_delta = ["#e05555" if v > 0 else "#4caf50" for v in df_mut["delta_global_score"]]

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(
    df_mut["label"],
    df_mut["delta_global_score"],
    color=colors_delta,
    width=0.6, edgecolor="white",
)
for bar, val in zip(bars, df_mut["delta_global_score"]):
    ypos = val + 0.0005 if val >= 0 else val - 0.002
    ax.text(bar.get_x() + bar.get_width() / 2, ypos,
            f"{val:+.4f}", ha="center", va="bottom", fontsize=9)

ax.axhline(0, color="#444", lw=1.2)
ax.set_ylabel("ΔScore vs WT (positive = worse fit)", fontsize=11)
ax.set_title("ProteinMPNN ΔScore relative to WT backbone\n(green = improved fit, red = reduced fit)", fontsize=12)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
fig.savefig(OUT_DIR / "mpnn_delta_score.png", dpi=150)
plt.close(fig)
print("Saved: mpnn_delta_score.png")
