from Bio.PDB import PDBParser, PDBIO
import os

INPUT = "data/raw/9UBQ.pdb"
OUTPUT = "data/cleaned/hbsag.pdb"

os.makedirs("data/cleaned", exist_ok=True)

KEEP_CHAINS = ["A"]


class SelectChains:
    def __init__(self, keep):
        self.keep = keep

    # IMPORTANT: must include ALL levels
    def accept_model(self, model):
        return True

    def accept_chain(self, chain):
        return chain.id in self.keep

    def accept_residue(self, residue):
        return True

    def accept_atom(self, atom):
        return True


parser = PDBParser(QUIET=True)
structure = parser.get_structure("hb", INPUT)

io = PDBIO()
io.set_structure(structure)
io.save(OUTPUT, SelectChains(KEEP_CHAINS))

print("Cleaned structure saved →", OUTPUT)