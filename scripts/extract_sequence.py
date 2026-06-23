from Bio.PDB import PDBParser, PPBuilder
import os

INPUT = "data/cleaned/hbsag.pdb"

if not os.path.exists(INPUT):
    raise FileNotFoundError("Cleaned PDB not found. Run clean_structure.py first.")

parser = PDBParser(QUIET=True)
structure = parser.get_structure("hb", INPUT)

ppb = PPBuilder()

sequence = ""

for pp in ppb.build_peptides(structure):
    sequence += str(pp.get_sequence())

if len(sequence) == 0:
    raise ValueError("No sequence extracted — check chain selection.")

os.makedirs("data/cleaned", exist_ok=True)

with open("data/cleaned/hbsag.fasta", "w") as f:
    f.write(">HBsAg\n")
    f.write(sequence)

print("Sequence length:", len(sequence))
print("Saved FASTA → data/cleaned/hbsag.fasta")