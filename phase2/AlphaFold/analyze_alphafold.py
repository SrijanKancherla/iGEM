"""
AlphaFold2 structural quality analysis for HBsAg WT and mutants.

Produces in phase2/AlphaFold/results/:
  plddt_comparison.png   -- per-residue pLDDT, all structures overlaid
  pae_heatmaps.png       -- 2x2 PAE matrices
  summary_metrics.png    -- pTM + mean pLDDT + local pLDDT bar chart
  rmsd_vs_wt.png         -- CA-RMSD of each mutant vs WT
  alphafold_summary.csv  -- all numeric metrics
"""

import glob
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB import Superimposer

AF_DIR = Path(__file__).parent
OUT_DIR = AF_DIR / "results"
OUT_DIR.mkdir(exist_ok=True)

MUTANTS = ["WT", "A76D", "A76K", "A90S"]

# AlphaFold renumbers residues 1-113 because the original PDB has a gap at 42-71.
# Mapping: original residue → AF residue
#   28-41  → 1-14  (no gap)
#   72-170 → 15-113 (original 72 = AF 15)
# So: original 76 = AF 19,  original 90 = AF 33
AF_RES = {"A76": 19, "A90": 33}

# 1-based AF residue numbers of mutation sites (for highlighting)
MUTATION_SITES = {
    "A76D": [AF_RES["A76"]], "A76K": [AF_RES["A76"]], "A90S": [AF_RES["A90"]],
    "WT": [],
}
COLORS = {"WT": "#444444", "A76D": "#e05555", "A76K": "#3a7ebf", "A90S": "#4caf50"}
HIGHLIGHT_COLORS = ["#f4a261", "#a8dadc"]


# ── helpers ─────────────────────────────────────────────────────────────────

def find_rank1_json(mut):
    pat = str(AF_DIR / mut / f"*scores_rank_001*.json")
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(f"No rank-001 scores JSON for {mut}")
    return hits[0]

def find_rank1_pdb(mut):
    pat = str(AF_DIR / mut / f"*relaxed_rank_001*.pdb")
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(f"No rank-001 relaxed PDB for {mut}")
    return hits[0]

def find_pae_json(mut):
    pat = str(AF_DIR / mut / f"*predicted_aligned_error_v1.json")
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(f"No PAE JSON for {mut}")
    return hits[0]

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


# ── load data ────────────────────────────────────────────────────────────────

plddt_data, ptm_data, pae_data, ca_data = {}, {}, {}, {}
for mut in MUTANTS:
    plddt_data[mut], ptm_data[mut] = load_plddt(mut)
    pae_data[mut] = load_pae(mut)
    ca_data[mut] = get_ca_atoms(find_rank1_pdb(mut))

n_res = len(plddt_data["WT"])
residues = np.arange(1, n_res + 1)
wt_cas = ca_data["WT"]

# RMSD vs WT for each mutant
rmsd_vs_wt = {mut: compute_rmsd(wt_cas, ca_data[mut]) for mut in MUTANTS if mut != "WT"}

# Local pLDDT around mutation sites (±5 residues, 0-indexed)
def local_plddt(mut, sites, half_win=5):
    arr = plddt_data[mut]
    vals = []
    for s in sites:
        idx = s - 1  # convert to 0-indexed
        lo, hi = max(0, idx - half_win), min(n_res, idx + half_win + 1)
        vals.extend(arr[lo:hi].tolist())
    return float(np.mean(vals)) if vals else float("nan")

local_plddt_data = {mut: local_plddt(mut, MUTATION_SITES.get(mut, [])) for mut in MUTANTS}
# For WT, average over both AF residue positions
local_plddt_data["WT"] = local_plddt("WT", [AF_RES["A76"], AF_RES["A90"]])


# ── Figure 1: per-residue pLDDT comparison ───────────────────────────────────

fig, ax = plt.subplots(figsize=(13, 5))
for mut in MUTANTS:
    lw = 2.2 if mut == "WT" else 1.6
    ax.plot(residues, plddt_data[mut], label=mut, color=COLORS[mut],
            linewidth=lw, alpha=0.9)

# Shade mutation site windows
for i, (site, hcol) in enumerate(zip([76, 90], HIGHLIGHT_COLORS)):
    ax.axvspan(site - 5, site + 5, alpha=0.15, color=hcol,
               label=f"Site ±5 (pos {site})")

ax.axhline(70, ls="--", color="gray", lw=1, alpha=0.7, label="pLDDT=70 threshold")
ax.axhline(50, ls=":", color="gray", lw=1, alpha=0.5, label="pLDDT=50 threshold")
ax.set_xlabel("Residue position", fontsize=12)
ax.set_ylabel("pLDDT", fontsize=12)
ax.set_title("Per-residue pLDDT: WT vs HBsAg mutants (AlphaFold2)", fontsize=13)
ax.set_xlim(1, n_res)
ax.set_ylim(0, 100)
ax.legend(loc="upper right", fontsize=9, ncol=2)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
fig.savefig(OUT_DIR / "plddt_comparison.png", dpi=150)
plt.close(fig)
print("Saved: plddt_comparison.png")


# ── Figure 2: PAE heatmaps (2x2) ─────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(11, 10))
for ax, mut in zip(axes.flat, MUTANTS):
    pae = pae_data[mut]
    im = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=30, aspect="auto")
    ax.set_title(f"{mut}  (pTM={ptm_data[mut]:.3f})", fontsize=11)
    ax.set_xlabel("Scored residue", fontsize=9)
    ax.set_ylabel("Aligned residue", fontsize=9)
    # Mark mutation sites
    for site in MUTATION_SITES.get(mut, []):
        ax.axvline(site - 1, color="red", lw=1.0, alpha=0.7)
        ax.axhline(site - 1, color="red", lw=1.0, alpha=0.7)
    plt.colorbar(im, ax=ax, fraction=0.046, label="PAE (Å)")

fig.suptitle("Predicted Aligned Error matrices (AlphaFold2)", fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(OUT_DIR / "pae_heatmaps.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: pae_heatmaps.png")


# ── Figure 3: summary bar chart ───────────────────────────────────────────────

metrics = {
    "pTM": {m: ptm_data[m] for m in MUTANTS},
    "Mean pLDDT / 100": {m: float(np.mean(plddt_data[m])) / 100 for m in MUTANTS},
    "Local pLDDT (site ±5) / 100": {m: local_plddt_data[m] / 100 for m in MUTANTS},
}

fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=False)
for ax, (title, vals) in zip(axes, metrics.items()):
    muts = list(vals.keys())
    ys = [vals[m] for m in muts]
    bars = ax.bar(muts, ys, color=[COLORS[m] for m in muts], width=0.55, edgecolor="white")
    for bar, y in zip(bars, ys):
        ax.text(bar.get_x() + bar.get_width() / 2, y + 0.005,
                f"{y:.3f}", ha="center", va="bottom", fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.set_ylim(0, max(ys) * 1.25)
    ax.set_ylabel("Score", fontsize=10)
    ax.grid(axis="y", alpha=0.3)

fig.suptitle("AlphaFold2 confidence summary: WT vs mutants", fontsize=13)
plt.tight_layout()
fig.savefig(OUT_DIR / "summary_metrics.png", dpi=150)
plt.close(fig)
print("Saved: summary_metrics.png")


# ── Figure 4: CA-RMSD vs WT ──────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 4))
muts_only = [m for m in MUTANTS if m != "WT"]
rmsds = [rmsd_vs_wt[m] for m in muts_only]
bars = ax.bar(muts_only, rmsds, color=[COLORS[m] for m in muts_only],
              width=0.5, edgecolor="white")
for bar, r in zip(bars, rmsds):
    ax.text(bar.get_x() + bar.get_width() / 2, r + 0.005,
            f"{r:.3f} Å", ha="center", va="bottom", fontsize=10)
ax.set_ylabel("CA-RMSD vs WT (Å)", fontsize=11)
ax.set_title("Backbone deviation from WT (AlphaFold2 structures)", fontsize=11)
ax.set_ylim(0, max(rmsds) * 1.35 if rmsds else 1)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
fig.savefig(OUT_DIR / "rmsd_vs_wt.png", dpi=150)
plt.close(fig)
print("Saved: rmsd_vs_wt.png")


# ── CSV summary ───────────────────────────────────────────────────────────────

rows = []
for mut in MUTANTS:
    rows.append({
        "mutant": mut,
        "pTM": ptm_data[mut],
        "mean_pLDDT": round(float(np.mean(plddt_data[mut])), 2),
        "local_pLDDT_site": round(local_plddt_data[mut], 2),
        "CA_RMSD_vs_WT": rmsd_vs_wt.get(mut, 0.0),
    })
df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "alphafold_summary.csv", index=False)
print("Saved: alphafold_summary.csv")

print("\n=== AlphaFold2 Summary ===")
print(df.to_string(index=False))
print(f"\nRMSD vs WT: {rmsd_vs_wt}")
