import os
import pandas as pd
import numpy as np

from pyrosetta import *
from pyrosetta.toolbox.mutants import mutate_residue

from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.minimization_packing import (
    PackRotamersMover,
    MinMover,
)

from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.kinematics import MoveMap

os.makedirs("results", exist_ok=True)

# ----------------------------------------
# INIT
# ----------------------------------------

init("-relax:default_repeats 1 -ex1 -ex2aro -mute all")

scorefxn = get_fa_scorefxn()

# ----------------------------------------
# LOAD WT
# ----------------------------------------

wt_pose = pose_from_pdb("data/cleaned/hbsag.pdb")

# ----------------------------------------
# WT ENSEMBLE
# ----------------------------------------

print("Generating WT ensemble...")

relax = FastRelax()
relax.set_scorefxn(scorefxn)

wt_scores = []
wt_relaxed_models = []

N_ENSEMBLE = 20

for i in range(N_ENSEMBLE):

    p = wt_pose.clone()

    relax.apply(p)

    wt_relaxed_models.append(p)
    wt_scores.append(scorefxn(p))

wt_score_mean = np.mean(wt_scores)
wt_score_std = np.std(wt_scores)

print(
    f"WT ensemble: mean={wt_score_mean:.2f} "
    f"std={wt_score_std:.2f}"
)

# ----------------------------------------
# LOCAL MINIMIZATION SETUP
# ----------------------------------------

mm = MoveMap()
mm.set_bb(False)
mm.set_chi(True)

min_mover = MinMover()
min_mover.movemap(mm)
min_mover.score_function(scorefxn)
min_mover.min_type("lbfgs_armijo_nonmonotone")

# ----------------------------------------
# LOAD MUTATIONS
# ----------------------------------------

df = pd.read_csv("data/mutations/esm_candidates.csv")

results = []

# ----------------------------------------
# MUTATION SCAN
# ----------------------------------------

print("Scanning mutations...")

for i, row in df.iterrows():

    # start from one relaxed WT model
    mut_pose = wt_relaxed_models[i % N_ENSEMBLE].clone()

    pos = int(row["pos"]) + 1
    mut = row["mut"]

    mutate_residue(mut_pose, pos, mut)

    tf = TaskFactory()
    pack_task = tf.create_packer_task(mut_pose)

    pack = PackRotamersMover(scorefxn, pack_task)
    pack.apply(mut_pose)

    min_mover.apply(mut_pose)

    mut_score = scorefxn(mut_pose)

    ddg = mut_score - wt_score_mean

    results.append(
        [
            row["pos"],
            row["wt"],
            row["mut"],
            ddg,
        ]
    )

    if i % 100 == 0:
        print(f"{i}/{len(df)} complete")

# ----------------------------------------
# SAVE
# ----------------------------------------

out = pd.DataFrame(
    results,
    columns=["pos", "wt", "mut", "ddg"]
)

out.to_csv("results/rosetta_ddg.csv", index=False)

print("Rosetta scan complete.")