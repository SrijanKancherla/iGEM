"""
AlphaFold2 structural quality analysis — v2 construct-based pipeline.

Reads from: construct_v2/alphafold_input/<MUTANT>/
Produces in: construct_v2/phase2/alphafold/results/
  plddt_comparison.png   -- per-residue pLDDT for all structures
  pae_heatmaps.png       -- PAE matrices (3 x 3 grid)
  summary_metrics.png    -- pTM + mean pLDDT + epitope pLDDT + RMSD bar charts
  rmsd_vs_wt.png         -- CA-RMSD of each mutant vs WT (full + epitope region)
  alphafold_summary.csv  -- all numeric metrics
  final_shortlist.csv    -- combined Phase1 + MPNN + AlphaFold ranking

Notes:
- ColabFold output is unrelaxed only; unrelaxed rank-001 is used throughout.
- Residue numbering: AlphaFold residues 1-226 = s_protein_pos 1-226 (direct).
- Epitope region: s_protein_pos 99-169 (MHR / a-determinant).
- Loop region: s_protein_pos 41-70 ("a" determinant disordered loop).

Run from project root:
    python construct_v2/phase2/alphafold/analyze_alphafold_v2.py
"""

import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB import Superimposer

ROOT    = Path(__file__).parent.parent.parent.parent   # project root
AF_DIR  = ROOT / "construct_v2/alphafold_input"
OUT_DIR = ROOT / "construct_v2/phase2/alphafold/results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MUTANTS = ["WT", "T57A", "H60T", "I68S", "S55T", "T114S", "V177L"]

# Direct s_protein_pos (= AF residue number) for each mutant's changed site
MUTATION_SITES = {
    "WT":    [],
    "T57A":  [57],
    "H60T":  [60],
    "I68S":  [68],
    "S55T":  [55],
    "T114S": [114],
    "V177L": [177],
}

COLORS = {
    "WT":    "#444444",
    "T57A":  "#3a7ebf",
    "H60T":  "#e05555",
    "I68S":  "#4caf50",
    "S55T":  "#9b59b6",
    "T114S": "#f39c12",
    "V177L": "#1abc9c",
}

# Regions (1-indexed, inclusive)
EPITOPE_START, EPITOPE_END = 99, 169   # MHR / a-determinant
LOOP_START,    LOOP_END    = 41,  70   # disordered loop


# ── file finders ────────────────────────────────────────────────────────────

def find_file(mut, pattern):
    hits = sorted(glob.glob(str(AF_DIR / mut / pattern)))
    if not hits:
        raise FileNotFoundError(f"No file matching {pattern!r} in {AF_DIR / mut}")
    return hits[0]

def find_rank1_pdb(mut):
    return find_file(mut, "*unrelaxed_rank_001*.pdb")

def find_rank1_json(mut):
    return find_file(mut, "*scores_rank_001*.json")

def find_pae_json(mut):
    return find_file(mut, "*predicted_aligned_error_v1.json")


# ── loaders ─────────────────────────────────────────────────────────────────

def load_plddt(mut):
    with open(find_rank1_json(mut)) as f:
        d = json.load(f)
    return np.array(d["plddt"]), float(d["ptm"])

def load_pae(mut):
    with open(find_pae_json(mut)) as f:
        d = json.load(f)
    return np.array(d["predicted_aligned_error"])

def get_ca_atoms(pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    cas = []
    for model in struct:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    cas.append(res["CA"])
        break
    return cas

def compute_rmsd(fixed_cas, moving_cas):
    sup = Superimposer()
    n = min(len(fixed_cas), len(moving_cas))
    sup.set_atoms(fixed_cas[:n], moving_cas[:n])
    sup.apply(moving_cas[:n])
    return round(float(sup.rms), 3)

def region_rmsd(fixed_cas, moving_cas, start, end):
    """CA-RMSD over a residue window (1-indexed, inclusive)."""
    lo, hi = start - 1, end  # convert to 0-indexed slice
    f = fixed_cas[lo:hi]
    m = moving_cas[lo:hi]
    sup = Superimposer()
    n = min(len(f), len(m))
    if n < 3:
        return float("nan")
    sup.set_atoms(f[:n], m[:n])
    sup.apply(m[:n])
    return round(float(sup.rms), 3)

def local_plddt(plddt_arr, sites, half_win=5):
    vals = []
    for s in sites:
        idx = s - 1
        lo = max(0, idx - half_win)
        hi = min(len(plddt_arr), idx + half_win + 1)
        vals.extend(plddt_arr[lo:hi].tolist())
    return float(np.mean(vals)) if vals else float("nan")

def region_plddt(plddt_arr, start, end):
    return float(np.mean(plddt_arr[start - 1:end]))


# ── load all data ────────────────────────────────────────────────────────────

print("Loading AlphaFold data...")
plddt_data, ptm_data, pae_data, ca_data = {}, {}, {}, {}
for mut in MUTANTS:
    plddt_data[mut], ptm_data[mut] = load_plddt(mut)
    pae_data[mut]  = load_pae(mut)
    ca_data[mut]   = get_ca_atoms(find_rank1_pdb(mut))
    print(f"  {mut}: {len(plddt_data[mut])} residues, pTM={ptm_data[mut]:.3f}, "
          f"mean pLDDT={np.mean(plddt_data[mut]):.1f}")

n_res    = len(plddt_data["WT"])
residues = np.arange(1, n_res + 1)
wt_cas   = ca_data["WT"]

muts_only = [m for m in MUTANTS if m != "WT"]
rmsd_full    = {m: compute_rmsd(wt_cas, ca_data[m]) for m in muts_only}
rmsd_epitope = {m: region_rmsd(wt_cas, ca_data[m], EPITOPE_START, EPITOPE_END)
                for m in muts_only}
rmsd_loop    = {m: region_rmsd(wt_cas, ca_data[m], LOOP_START, LOOP_END)
                for m in muts_only}


# ── Figure 1: per-residue pLDDT ─────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(14, 5))
for mut in MUTANTS:
    lw = 2.4 if mut == "WT" else 1.5
    alpha = 1.0 if mut == "WT" else 0.85
    ax.plot(residues, plddt_data[mut], label=mut, color=COLORS[mut],
            linewidth=lw, alpha=alpha)

ax.axvspan(LOOP_START,    LOOP_END,    alpha=0.07, color="#f4a261", label="a-det loop (41–70)")
ax.axvspan(EPITOPE_START, EPITOPE_END, alpha=0.07, color="#a8dadc", label="MHR (99–169)")
ax.axhline(70, ls="--", color="gray", lw=1, alpha=0.6, label="pLDDT=70")
ax.axhline(50, ls=":",  color="gray", lw=1, alpha=0.4, label="pLDDT=50")

ax.set_xlabel("Residue (s_protein_pos)", fontsize=12)
ax.set_ylabel("pLDDT", fontsize=12)
ax.set_title("Per-residue pLDDT: WT vs shortlisted mutants (AlphaFold2)", fontsize=13)
ax.set_xlim(1, n_res)
ax.set_ylim(0, 100)
ax.legend(loc="upper right", fontsize=8, ncol=3)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
fig.savefig(OUT_DIR / "plddt_comparison.png", dpi=150)
plt.close(fig)
print("Saved: plddt_comparison.png")


# ── Figure 2: PAE heatmaps (3 x 3 grid, last cell blank) ────────────────────

ncols = 3
nrows = (len(MUTANTS) + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(13, nrows * 4.2))
axes_flat = axes.flat

for ax, mut in zip(axes_flat, MUTANTS):
    pae = pae_data[mut]
    im  = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=30, aspect="auto")
    ax.set_title(f"{mut}  (pTM={ptm_data[mut]:.3f})", fontsize=10)
    ax.set_xlabel("Scored residue", fontsize=8)
    ax.set_ylabel("Aligned residue", fontsize=8)
    for site in MUTATION_SITES.get(mut, []):
        ax.axvline(site - 1, color="red", lw=1.0, alpha=0.8)
        ax.axhline(site - 1, color="red", lw=1.0, alpha=0.8)
    plt.colorbar(im, ax=ax, fraction=0.046, label="PAE (Å)")

# hide leftover axes
for ax in list(axes_flat)[len(MUTANTS):]:
    ax.set_visible(False)

fig.suptitle("Predicted Aligned Error (AlphaFold2) — v2 construct", fontsize=13, y=1.01)
plt.tight_layout()
fig.savefig(OUT_DIR / "pae_heatmaps.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: pae_heatmaps.png")


# ── Figure 3: summary bar charts ─────────────────────────────────────────────

metrics = {
    "pTM":                      {m: ptm_data[m] for m in MUTANTS},
    "Mean pLDDT / 100":         {m: float(np.mean(plddt_data[m])) / 100 for m in MUTANTS},
    "Epitope pLDDT (99–169)/100": {m: region_plddt(plddt_data[m], EPITOPE_START, EPITOPE_END) / 100
                                   for m in MUTANTS},
    "Loop pLDDT (41–70)/100":   {m: region_plddt(plddt_data[m], LOOP_START, LOOP_END) / 100
                                  for m in MUTANTS},
}

fig, axes = plt.subplots(1, 4, figsize=(16, 5))
for ax, (title, vals) in zip(axes, metrics.items()):
    muts_list = list(vals.keys())
    ys = [vals[m] for m in muts_list]
    bars = ax.bar(muts_list, ys, color=[COLORS[m] for m in muts_list],
                  width=0.6, edgecolor="white")
    for bar, y in zip(bars, ys):
        ax.text(bar.get_x() + bar.get_width() / 2, y + 0.003,
                f"{y:.3f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_title(title, fontsize=9)
    ax.set_ylim(0, min(1.0, max(ys) * 1.22))
    ax.set_xticks(range(len(muts_list)))
    ax.set_xticklabels(muts_list, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Score", fontsize=9)
    ax.grid(axis="y", alpha=0.3)

fig.suptitle("AlphaFold2 confidence summary — v2 construct mutants", fontsize=12)
plt.tight_layout()
fig.savefig(OUT_DIR / "summary_metrics.png", dpi=150)
plt.close(fig)
print("Saved: summary_metrics.png")


# ── Figure 4: CA-RMSD vs WT ─────────────────────────────────────────────────

x = np.arange(len(muts_only))
width = 0.28
fig, ax = plt.subplots(figsize=(9, 5))
bars1 = ax.bar(x - width, [rmsd_full[m]    for m in muts_only], width,
               label="Full structure", color=[COLORS[m] for m in muts_only],
               alpha=0.9, edgecolor="white")
bars2 = ax.bar(x,          [rmsd_epitope[m] for m in muts_only], width,
               label="Epitope (99–169)", color=[COLORS[m] for m in muts_only],
               alpha=0.55, edgecolor="white", hatch="//")
bars3 = ax.bar(x + width, [rmsd_loop[m]    for m in muts_only], width,
               label="Loop (41–70)", color=[COLORS[m] for m in muts_only],
               alpha=0.35, edgecolor="white", hatch="xx")

for bars in [bars1, bars2, bars3]:
    for bar in bars:
        h = bar.get_height()
        if not np.isnan(h):
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.002,
                    f"{h:.2f}", ha="center", va="bottom", fontsize=7)

ax.set_xticks(x)
ax.set_xticklabels(muts_only, fontsize=10)
ax.set_ylabel("CA-RMSD vs WT (Å)", fontsize=11)
ax.set_title("Backbone deviation from WT — full / epitope / loop regions", fontsize=11)
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
fig.savefig(OUT_DIR / "rmsd_vs_wt.png", dpi=150)
plt.close(fig)
print("Saved: rmsd_vs_wt.png")


# ── CSV summary + final shortlist ────────────────────────────────────────────

rows = []
for mut in MUTANTS:
    rows.append({
        "mutant":               mut,
        "pTM":                  round(ptm_data[mut], 4),
        "mean_pLDDT":           round(float(np.mean(plddt_data[mut])), 2),
        "epitope_pLDDT":        round(region_plddt(plddt_data[mut], EPITOPE_START, EPITOPE_END), 2),
        "loop_pLDDT":           round(region_plddt(plddt_data[mut], LOOP_START, LOOP_END), 2),
        "site_pLDDT":           round(local_plddt(plddt_data[mut], MUTATION_SITES.get(mut, []), 5), 2),
        "CA_RMSD_full":         rmsd_full.get(mut, 0.0),
        "CA_RMSD_epitope":      rmsd_epitope.get(mut, 0.0),
        "CA_RMSD_loop":         rmsd_loop.get(mut, 0.0),
    })

af_df = pd.DataFrame(rows)
af_df.to_csv(OUT_DIR / "alphafold_summary.csv", index=False)
print("Saved: alphafold_summary.csv")

# --- merge with Phase 1 + MPNN scores ---
ph1 = pd.read_csv(ROOT / "construct_v2/results/final_ranked_simple.csv")
mpnn_scores = {
    "T57A":  -0.034376, "H60T":  -0.007140, "I68S":  -0.092927,
    "S55T":   0.015971, "T114S":  0.041130, "V177L":  0.024480,
}
ph1_sub = ph1[ph1["mutation"].isin(muts_only)].set_index("mutation")

shortlist_rows = []
for mut in muts_only:
    p = ph1_sub.loc[mut] if mut in ph1_sub.index else {}
    a = af_df[af_df["mutant"] == mut].iloc[0]
    wt_a = af_df[af_df["mutant"] == "WT"].iloc[0]

    # pLDDT delta vs WT at the epitope and loop regions
    epi_delta  = round(a["epitope_pLDDT"]  - wt_a["epitope_pLDDT"], 2)
    loop_delta = round(a["loop_pLDDT"]     - wt_a["loop_pLDDT"],    2)
    site_delta = round(a["site_pLDDT"]     - wt_a.get("mean_pLDDT", wt_a["mean_pLDDT"]), 2)

    shortlist_rows.append({
        "mutation":          mut,
        "region":            ("MHR" if 99 <= int(p.get("s_protein_pos", 0)) <= 169
                              else "a-det loop" if 41 <= int(p.get("s_protein_pos", 0)) <= 70
                              else "outside epitope"),
        "s_protein_pos":     p.get("s_protein_pos", ""),
        "esm_score":         round(float(p.get("esm_score", 0)), 3),
        "rosetta_ddg":       round(float(p.get("ddg", 0)), 3),
        "mpnn_delta":        mpnn_scores.get(mut, float("nan")),
        "ph1_final_score":   round(float(p.get("final_score", 0)), 4),
        "pTM":               a["pTM"],
        "mean_pLDDT":        a["mean_pLDDT"],
        "epitope_pLDDT":     a["epitope_pLDDT"],
        "loop_pLDDT":        a["loop_pLDDT"],
        "site_pLDDT":        a["site_pLDDT"],
        "epi_pLDDT_vs_WT":   epi_delta,
        "loop_pLDDT_vs_WT":  loop_delta,
        "CA_RMSD_full":      a["CA_RMSD_full"],
        "CA_RMSD_epitope":   a["CA_RMSD_epitope"],
        "CA_RMSD_loop":      a["CA_RMSD_loop"],
        "ph1_recommendation": p.get("recommendation", ""),
    })

sl_df = pd.DataFrame(shortlist_rows)

# Final combined score (higher = better):
#   40% Phase1 score | 20% pTM | 20% mean pLDDT (norm) | 20% MPNN fit
pTM_norm     = (sl_df["pTM"] - sl_df["pTM"].min()) / (sl_df["pTM"].max() - sl_df["pTM"].min() + 1e-9)
plddt_norm   = (sl_df["mean_pLDDT"] - sl_df["mean_pLDDT"].min()) / (sl_df["mean_pLDDT"].max() - sl_df["mean_pLDDT"].min() + 1e-9)
mpnn_arr     = sl_df["mpnn_delta"].values
mpnn_norm    = 1.0 - (mpnn_arr - mpnn_arr.min()) / (mpnn_arr.max() - mpnn_arr.min() + 1e-9)

sl_df["combined_score"] = (
    0.40 * sl_df["ph1_final_score"]
    + 0.20 * pTM_norm
    + 0.20 * plddt_norm
    + 0.20 * mpnn_norm
).round(4)

sl_df = sl_df.sort_values("combined_score", ascending=False).reset_index(drop=True)
sl_df.insert(0, "rank", range(1, len(sl_df) + 1))
sl_df.to_csv(OUT_DIR / "final_shortlist.csv", index=False)
print("Saved: final_shortlist.csv")

print("\n=== AlphaFold v2 Summary (vs WT) ===")
print(af_df.to_string(index=False))

print("\n=== FINAL SHORTLIST (combined ranking) ===")
display_cols = ["rank", "mutation", "region", "rosetta_ddg", "mpnn_delta",
                "ph1_final_score", "pTM", "mean_pLDDT", "epitope_pLDDT",
                "loop_pLDDT", "CA_RMSD_full", "CA_RMSD_loop", "combined_score"]
print(sl_df[display_cols].to_string(index=False))
