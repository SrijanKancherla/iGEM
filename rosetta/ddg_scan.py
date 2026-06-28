import os
import pandas as pd
import numpy as np

from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.toolbox.mutants import mutate_residue

from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.core.kinematics import MoveMap

os.makedirs("results", exist_ok=True)

init("-relax:default_repeats 1 -ex1 -ex2aro -mute all")

scorefxn = get_fa_scorefxn()

# ----------------------------
# LOAD + RELAX WT
# ----------------------------
wt_pose = pose_from_pdb("data/cleaned/hbsag.pdb")

print("Relaxing WT...")
relax = FastRelax()
relax.set_scorefxn(scorefxn)

wt_relaxed = wt_pose.clone()
relax.apply(wt_relaxed)
print(f"WT score after FastRelax: {scorefxn(wt_relaxed):.3f}")

# ----------------------------
# REUSABLE MOVERS
# ----------------------------
tf = TaskFactory()
tf.push_back(RestrictToRepacking())

mm = MoveMap()
mm.set_bb(False)
mm.set_chi(True)

pack_mover = PackRotamersMover(scorefxn)
pack_mover.task_factory(tf)

min_mover = MinMover()
min_mover.movemap(mm)
min_mover.score_function(scorefxn)
min_mover.min_type("lbfgs_armijo_nonmonotone")
min_mover.tolerance(0.1)

# Compute a fair WT reference by applying the same pack+min as mutants.
# Without this, WT (FastRelax) and mutants (pack+min) are on different energy
# scales and ddG is systematically inflated by ~15-30 REU.
wt_ref = wt_relaxed.clone()
pack_mover.apply(wt_ref)
min_mover.apply(wt_ref)
wt_ref_score = scorefxn(wt_ref)
print(f"WT pack+min reference score: {wt_ref_score:.3f}")

# ----------------------------
# LOAD MUTATIONS
# ----------------------------
mutation_csv = "data/mutations/esm2_mutations.csv"
if not os.path.exists(mutation_csv):
    mutation_csv = "data/mutations/esm_candidates.csv"

df = pd.read_csv(mutation_csv)
print(f"Loaded {len(df)} mutations from {mutation_csv}")

# ----------------------------
# DDG SCAN
# ----------------------------
print("Running mutation scan...")
results = []

for i, row in df.iterrows():
    mut_pose = wt_relaxed.clone()
    pos = int(row["pos"]) + 1  # Rosetta is 1-indexed
    mutate_residue(mut_pose, pos, row["mut"])
    pack_mover.apply(mut_pose)
    min_mover.apply(mut_pose)
    ddg = scorefxn(mut_pose) - wt_ref_score
    results.append([row["pos"], row["wt"], row["mut"], row.get("esm_score", np.nan), ddg])
    if i % 100 == 0:
        print(f"  {i}/{len(df)} mutations processed")

# ----------------------------
# SAVE
# ----------------------------
out = pd.DataFrame(results, columns=["pos", "wt", "mut", "esm_score", "ddg"])
out.to_csv("results/rosetta_ddg.csv", index=False)
print(f"DDG scan complete. Range: [{out['ddg'].min():.2f}, {out['ddg'].max():.2f}]")
