"""
Rosetta ddG scanning for the v2 construct-based pipeline.

Requires: construct_v2/data/cleaned/alphafold_construct.pdb
  - Submit construct_v2/alphafold_input/hbsag_for_alphafold.fasta to AlphaFold
    (AlphaFold Server: https://alphafoldserver.com, or ColabFold, or local AF2).
  - Save the best-ranked relaxed model as:
      construct_v2/data/cleaned/alphafold_construct.pdb
  - Chain A must contain residues numbered 1-226 matching the 226 aa HBsAg sequence.

Numbering carried through from esm2_mutations_v2:
  pos (0-indexed in 226 aa) → Rosetta residue = pos + 1 (Rosetta is 1-indexed, AF2 PDB starts at 1)
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

PDB_V2       = Path("construct_v2/data/cleaned/alphafold_construct.pdb")
MUTATION_CSV = Path("construct_v2/mutations/esm2_mutations.csv")
OUTPUT_CSV   = Path("construct_v2/results/rosetta_ddg.csv")


def check_prerequisites():
    missing = []
    if not PDB_V2.exists():
        missing.append(
            f"\n  ✗  {PDB_V2}\n"
            "     → Run AlphaFold on construct_v2/alphafold_input/hbsag_for_alphafold.fasta\n"
            "       and save the relaxed rank-1 model to the path above."
        )
    if not MUTATION_CSV.exists():
        missing.append(
            f"\n  ✗  {MUTATION_CSV}\n"
            "     → Run:  python esm/esm2_mutations_v2.py"
        )
    if missing:
        print("ddg_scan_v2.py: missing prerequisites:" + "".join(missing))
        sys.exit(1)


def main():
    check_prerequisites()

    from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
    from pyrosetta.toolbox.mutants import mutate_residue
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
    from pyrosetta.rosetta.core.kinematics import MoveMap

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    init("-relax:default_repeats 1 -ex1 -ex2aro -mute all")
    scorefxn = get_fa_scorefxn()

    print(f"Loading WT structure: {PDB_V2}")
    wt_pose = pose_from_pdb(str(PDB_V2))

    print("Relaxing WT...")
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    wt_relaxed = wt_pose.clone()
    relax.apply(wt_relaxed)

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

    wt_ref = wt_relaxed.clone()
    pack_mover.apply(wt_ref)
    min_mover.apply(wt_ref)
    wt_ref_score = scorefxn(wt_ref)
    print(f"WT pack+min reference score: {wt_ref_score:.3f}")

    df = pd.read_csv(MUTATION_CSV)
    print(f"Loaded {len(df)} mutations from {MUTATION_CSV}")

    results = []
    for i, row in df.iterrows():
        mut_pose = wt_relaxed.clone()
        rosetta_pos = int(row["pos"]) + 1  # Rosetta is 1-indexed; AF2 PDB residues start at 1
        mutate_residue(mut_pose, rosetta_pos, row["mut"])
        pack_mover.apply(mut_pose)
        min_mover.apply(mut_pose)
        ddg = scorefxn(mut_pose) - wt_ref_score
        results.append({
            "pos":           int(row["pos"]),
            "s_protein_pos": int(row["s_protein_pos"]),
            "construct_pos": int(row["construct_pos"]),
            "wt":            row["wt"],
            "mut":           row["mut"],
            "esm_score":     row.get("esm_score", np.nan),
            "ddg":           ddg,
        })
        if i % 100 == 0:
            print(f"  {i}/{len(df)} mutations processed")

    out = pd.DataFrame(results)
    out.to_csv(OUTPUT_CSV, index=False)
    print(f"DDG scan complete → {OUTPUT_CSV}")
    print(f"DDG range: [{out['ddg'].min():.2f}, {out['ddg'].max():.2f}]")


if __name__ == "__main__":
    main()
