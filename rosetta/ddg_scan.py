from pyrosetta import *
from pyrosetta.toolbox.mutants import mutate_residue
import pandas as pd
import os

os.makedirs("results", exist_ok=True)

init()

pose = pose_from_pdb("data/cleaned/hbsag.pdb")

scorefxn = get_fa_scorefxn()

wt_pose = pose.clone()
scorefxn(wt_pose)
wt_energy = wt_pose.energies().total_energy()

df = pd.read_csv("data/mutations/esm_candidates.csv")

results = []

for i, row in df.iterrows():

    mut_pose = pose.clone()

    # Rosetta is 1-indexed
    mutate_residue(mut_pose, int(row["pos"]) + 1, row["mut"])

    scorefxn(mut_pose)
    mut_energy = mut_pose.energies().total_energy()

    ddg = mut_energy - wt_energy

    results.append([row["pos"], row["wt"], row["mut"], ddg])

out = pd.DataFrame(results, columns=["pos", "wt", "mut", "ddg"])
out.to_csv("results/rosetta_ddg.csv", index=False)

print("Rosetta ddG scan complete.")