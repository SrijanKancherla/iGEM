import esm
import torch
import pandas as pd
import os

os.makedirs("data/mutations", exist_ok=True)

# load model
model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
model.eval()

batch_converter = alphabet.get_batch_converter()

# load sequence
seq = open("data/cleaned/hbsag.fasta").read().split("\n")[1]

data = [("hb", seq)]
labels, strs, tokens = batch_converter(data)

with torch.no_grad():
    out = model(tokens, repr_layers=[12])

logits = out["logits"][0]

aas = list("ACDEFGHIKLMNPQRSTVWY")

rows = []

for i in range(len(seq)):
    wt = seq[i]

    for aa in aas:
        if aa == wt:
            continue

        score = logits[i, alphabet.get_idx(aa)].item()
        rows.append([i, wt, aa, score])

df = pd.DataFrame(rows, columns=["pos", "wt", "mut", "esm_score"])
df.to_csv("data/mutations/esm_candidates.csv", index=False)

print("ESM mutation list saved.")