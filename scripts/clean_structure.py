from Bio.PDB import PDBParser, PDBIO
import os

INPUT = "data/raw/9UBQ.pdb"
OUTPUT = "data/cleaned/hbsag.pdb"

os.makedirs("data/cleaned", exist_ok=True)

KEEP_CHAINS = ["A"]  # you may adjust after inspection


class ChainSelect:
    def __init__(self, keep):
        self.keep = keep

    def accept_chain(self, chain):
        return chain.id in self.keep


parser = PDBParser(QUIET=True)
structure = parser.get_structure("hb", INPUT)

io = PDBIO()
io.set_structure(structure)
io.save(OUTPUT, ChainSelect(KEEP_CHAINS))

print("Cleaned structure saved.")