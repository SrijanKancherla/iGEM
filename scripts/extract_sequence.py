from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser(QUIET=True)
structure = parser.get_structure("hb", "data/cleaned/hbsag.pdb")

ppb = PPBuilder()

sequence = ""
for pp in ppb.build_peptides(structure):
    sequence += str(pp.get_sequence())

with open("data/cleaned/hbsag.fasta", "w") as f:
    f.write(">HBsAg\n")
    f.write(sequence)

print(sequence)
print("Sequence saved.")