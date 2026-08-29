"""
Build MSA for the v2 construct-based pipeline.

Query: construct_v2/data/cleaned/hbsag_construct.fasta (226 aa HBsAg)
Output: construct_v2/data/msa/hbsag_construct_msa.fasta

Fetches diverse HBsAg sequences from NCBI and projects each into the 226 aa
reference frame so every row is exactly 226 columns (no gaps in reference row).
This format is required by analysis/conservation_msa_score_v2.py.

Run once before the conservation scoring step:
    python scripts/build_msa_v2.py
"""

import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

Entrez.email = "srijankk@uio.no"

WT_FASTA   = Path("construct_v2/data/cleaned/hbsag_construct.fasta")
OUTPUT_MSA = Path("construct_v2/data/msa/hbsag_construct_msa.fasta")

MIN_IDENTITY = 0.25   # lower than v1 because we include more of the protein
MIN_SEQS     = 10


def read_wt():
    for rec in SeqIO.parse(WT_FASTA, "fasta"):
        return str(rec.seq).upper().replace("-", "")
    raise ValueError(f"No sequence in {WT_FASTA}")


def fetch_ids(max_ids=600):
    print("Searching NCBI protein for HBsAg sequences...")
    handle = Entrez.esearch(
        db="protein",
        term=(
            'hepatitis B virus[Organism] AND surface[All Fields]'
            ' AND 150:400[SLEN]'   # wider range: 226 aa + room for preS2 variants
        ),
        retmax=max_ids,
    )
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    print(f"  Found {len(ids)} candidate accessions")
    return ids


def fetch_sequences(ids, batch_size=50):
    sequences = []
    total = len(ids)
    for start in range(0, total, batch_size):
        batch = ids[start:start + batch_size]
        handle = Entrez.efetch(db="protein", id=",".join(batch),
                               rettype="fasta", retmode="text")
        for rec in SeqIO.parse(handle, "fasta"):
            seq = (str(rec.seq).upper()
                   .replace("-", "").replace("X", "").replace("*", ""))
            if 150 <= len(seq) <= 400:
                sequences.append((rec.id, seq))
        handle.close()
        time.sleep(0.34)
        fetched = min(start + batch_size, total)
        print(f"  Fetched {fetched}/{total} ({len(sequences)} valid)")
    return sequences


def align_to_ref_frame(ref, query, aligner):
    try:
        aln = next(aligner.align(ref, query))
        result = "".join(q for r, q in zip(str(aln[0]), str(aln[1])) if r != "-")
        return result if len(result) == len(ref) else None
    except Exception:
        return None


def build_msa(wt_seq, fetched):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    msa = [("WT_HBsAg_construct_v2", wt_seq)]
    seen = {wt_seq}
    skipped_identity = skipped_align = 0

    for seq_id, seq in fetched:
        projected = align_to_ref_frame(wt_seq, seq, aligner)
        if projected is None:
            skipped_align += 1
            continue
        identity = sum(a == b for a, b in zip(wt_seq, projected)) / len(wt_seq)
        if identity < MIN_IDENTITY:
            skipped_identity += 1
            continue
        if projected in seen:
            continue
        seen.add(projected)
        msa.append((seq_id, projected))

    print(f"  Skipped {skipped_align} (alignment failed), "
          f"{skipped_identity} (identity < {MIN_IDENTITY:.0%})")
    return msa


def main():
    OUTPUT_MSA.parent.mkdir(parents=True, exist_ok=True)
    wt_seq = read_wt()
    print(f"WT sequence: {len(wt_seq)} aa")

    ids = fetch_ids(max_ids=600)
    if not ids:
        print("No accessions — check internet / NCBI status.")
        sys.exit(1)

    fetched = fetch_sequences(ids)
    print(f"\nAligning {len(fetched)} sequences to WT reference frame...")
    msa = build_msa(wt_seq, fetched)
    depth = len(msa)
    print(f"MSA depth: {depth} sequences  (alignment length: {len(wt_seq)} columns)")

    if depth < MIN_SEQS:
        print(f"WARNING: only {depth} sequences — conservation scores will be unreliable.")

    with open(OUTPUT_MSA, "w") as fh:
        for seq_id, seq in msa:
            fh.write(f">{seq_id}\n{seq}\n")
    print(f"Saved: {OUTPUT_MSA}")


if __name__ == "__main__":
    main()
