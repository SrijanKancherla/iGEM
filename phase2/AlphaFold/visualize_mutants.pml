# ============================================================
# HBsAg WT vs A76D / A76K / A90S  —  AlphaFold2 structural comparison
#
# Run from the repo root:
#   pymol phase2/AlphaFold/visualize_mutants.pml
# Or paste into an open PyMOL session with: @phase2/AlphaFold/visualize_mutants.pml
#
# Residue numbering note:
#   AlphaFold renumbers 1-113 (the original PDB has a gap at 42-71).
#   A76 (original) = residue 19 in AF structures
#   A90 (original) = residue 33 in AF structures
#   B-factor column = pLDDT confidence score
# ============================================================


# ── 1. Load structures ───────────────────────────────────────────────────────

load phase2/AlphaFold/WT/WT_9caa6_relaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb,   wt
load phase2/AlphaFold/A76D/A76D_ca4da_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb, a76d
load phase2/AlphaFold/A76K/A76K_f692f_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb, a76k
load phase2/AlphaFold/A90S/A90S_ed36d_relaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb, a90s


# ── 2. Superpose all mutants onto WT (CA atoms) ──────────────────────────────

super a76d, wt
super a76k, wt
super a90s, wt


# ── 3. Representation — cartoon backbone ─────────────────────────────────────

hide everything
show cartoon, all
set cartoon_transparency, 0.0

# Custom colours per structure
set_color col_wt,   [0.27, 0.27, 0.27]   # dark grey
set_color col_a76d, [0.88, 0.33, 0.33]   # red
set_color col_a76k, [0.23, 0.49, 0.75]   # blue
set_color col_a90s, [0.30, 0.69, 0.31]   # green

color col_wt,   wt
color col_a76d, a76d
color col_a76k, a76k
color col_a90s, a90s


# ── 4. Mutation-site selections (AF numbering) ────────────────────────────────

# Position A76 (AF residue 19) — mutated in A76D and A76K
select site76, (wt or a76d or a76k) and resi 19
# Position A90 (AF residue 33) — mutated in A90S
select site90, (wt or a90s) and resi 33
# Both sites together (for WT reference)
select both_sites, wt and (resi 19 or resi 33)

# Window of ±5 residues around each site
select loop76, all and resi 14-24
select loop90, all and resi 28-38


# ── 5. Scene 1: Full overlay coloured by structure ────────────────────────────

zoom all
orient all
set_view [\
     0.97,  0.14, -0.20,\
    -0.13,  0.99,  0.08,\
     0.21, -0.05,  0.98,\
     0.00,  0.00, -180.0,\
     0.00,  0.00,   0.00,\
   100.0, 300.0,   0.0]

scene S1_full_overlay, store, "Full backbone overlay: WT (grey) A76D (red) A76K (blue) A90S (green)"


# ── 6. Scene 2: pLDDT confidence colouring (B-factor) ────────────────────────

# Colour all structures by B-factor (= pLDDT); blue=low, white=70, red=100
spectrum b, blue_white_red, all, minimum=50, maximum=100
scene S2_pLDDT, store, "All structures coloured by pLDDT (blue<50 → white=70 → red=100)"

# Reset to per-structure colours for subsequent scenes
color col_wt,   wt
color col_a76d, a76d
color col_a76k, a76k
color col_a90s, a90s


# ── 7. Scene 3: Zoom on A76 site — sticks for all four structures ─────────────

# Show sticks only at the mutation site and ±2 neighbours
select zoom76, all and resi 17-21
show sticks, zoom76
set stick_transparency, 0.0
zoom zoom76, 8
orient zoom76

# Label the mutant residues with their identity
label (zoom76 and name CA and a76d), "A76D (res19)"
label (zoom76 and name CA and a76k), "A76K (res19)"
label (zoom76 and name CA and wt),   "WT-Ala (res19)"
set label_size, -0.4

scene S3_site76, store, "Close-up: A76 mutation site (AF res 19). WT=grey A76D=red A76K=blue"

# Clear labels before next scene
label (zoom76 and name CA), ""
hide sticks, zoom76


# ── 8. Scene 4: Zoom on A90 site ─────────────────────────────────────────────

select zoom90, all and resi 31-35
show sticks, zoom90
zoom zoom90, 8
orient zoom90

label (zoom90 and name CA and a90s), "A90S (res33)"
label (zoom90 and name CA and wt),   "WT-Ala (res33)"

scene S4_site90, store, "Close-up: A90 mutation site (AF res 33). WT=grey A90S=green"

label (zoom90 and name CA), ""
hide sticks, zoom90


# ── 9. Scene 5: WT alone coloured by pLDDT ───────────────────────────────────

hide everything, a76d
hide everything, a76k
hide everything, a90s
show cartoon, wt
spectrum b, blue_white_red, wt, minimum=50, maximum=100

# Highlight both mutation sites as sticks
show sticks, wt and (resi 19 or resi 33)
zoom wt
orient wt

scene S5_WT_pLDDT, store, "WT pLDDT map. Sticks at res 19 (A76) and res 33 (A90)"

# Restore all
show cartoon, all
color col_wt,   wt
color col_a76d, a76d
color col_a76k, a76k
color col_a90s, a90s
hide sticks, all


# ── 10. Final view: back to full overlay ─────────────────────────────────────

scene S1_full_overlay, recall
zoom all

# Print reminder
print ""
print "=== Scenes ==="
print "  S1_full_overlay  : backbone overlay, one colour per structure"
print "  S2_pLDDT         : all structures coloured by pLDDT"
print "  S3_site76        : zoom on A76 mutation site (AF res 19)"
print "  S4_site90        : zoom on A90 mutation site (AF res 33)"
print "  S5_WT_pLDDT      : WT pLDDT map with mutation sites as sticks"
print ""
print "Recall a scene with:  scene <name>, recall"
print "Cycle scenes with:    F1 / F2 / F3 ..."
