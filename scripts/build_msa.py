#!/usr/bin/env python3
"""
Fetch diverse HBsAg sequences from NCBI and build an MSA in WT reference frame.

Strategy: pairwise-align each fetched sequence to the WT reference, then
project it into the reference coordinate frame (no gaps in WT row). The
result is a FASTA where every sequence has exactly len(WT) columns —
the format expected by analysis/conservation_msa_score.py.
"""

import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

Entrez.email = "srijankk@uio.no"

WT_FASTA = Path("data/cleaned/hbsag.fasta")
OUTPUT_MSA = Path("data/msa/hbsag_msa.fasta")

MIN_IDENTITY = 0.30   # drop sequences <30% identical to WT
MIN_SEQS = 10         # warn if fewer than this end up in MSA


def read_wt():
    for rec in SeqIO.parse(WT_FASTA, "fasta"):
        return str(rec.seq).upper().replace("-", "")
    raise ValueError(f"No sequence found in {WT_FASTA}")


def fetch_ids(max_ids=500):
    """Search NCBI protein for HBsAg sequences in a suitable length range."""
    print("Searching NCBI protein database for HBsAg sequences...")
    handle = Entrez.esearch(
        db="protein",
        term=(
            'hepatitis B virus[Organism] AND surface[All Fields]'
            ' AND 100:300[SLEN]'
        ),
        retmax=max_ids,
    )
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    print(f"  Found {len(ids)} candidate accessions")
    return ids


def fetch_sequences(ids, batch_size=50):
    """Fetch protein sequences from NCBI in batches, return (id, seq) pairs."""
    sequences = []
    total = len(ids)
    for start in range(0, total, batch_size):
        batch = ids[start : start + batch_size]
        handle = Entrez.efetch(
            db="protein",
            id=",".join(batch),
            rettype="fasta",
            retmode="text",
        )
        for rec in SeqIO.parse(handle, "fasta"):
            seq = (
                str(rec.seq)
                .upper()
                .replace("-", "")
                .replace("X", "")
                .replace("*", "")
            )
            if 80 <= len(seq) <= 300:
                sequences.append((rec.id, seq))
        handle.close()
        time.sleep(0.34)  # NCBI rate limit: max 3 requests/sec
        fetched = min(start + batch_size, total)
        print(f"  Fetched {fetched}/{total} accessions  ({len(sequences)} valid sequences)")
    return sequences


def align_to_ref_frame(ref, query, aligner):
    """
    Align query to ref globally, then return query projected into the reference
    coordinate frame — one character per non-gap position in the reference.
    Result is always exactly len(ref) characters.
    """
    try:
        aln = next(aligner.align(ref, query))
        aligned_ref = str(aln[0])
        aligned_query = str(aln[1])
        result = "".join(q for r, q in zip(aligned_ref, aligned_query) if r != "-")
        if len(result) != len(ref):
            return None
        return result
    except Exception:
        return None


def build_msa(wt_seq, fetched):
    """Align all fetched sequences to WT and return the MSA list."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    # WT is always the first (reference) row
    msa = [("WT_HBsAg_reference", wt_seq)]
    seen = {wt_seq}
    skipped_identity = 0
    skipped_align = 0

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

    print(f"  Skipped {skipped_align} sequences (alignment failed)")
    print(f"  Skipped {skipped_identity} sequences (identity < {MIN_IDENTITY:.0%})")
    return msa


def main():
    OUTPUT_MSA.parent.mkdir(parents=True, exist_ok=True)

    wt_seq = read_wt()
    print(f"WT sequence: {len(wt_seq)} aa")

    ids = fetch_ids(max_ids=500)
    if not ids:
        print("No accessions returned — check internet connection or NCBI status.")
        sys.exit(1)

    fetched = fetch_sequences(ids)
    print(f"\nAligning {len(fetched)} sequences to WT reference frame...")
    msa = build_msa(wt_seq, fetched)

    depth = len(msa)
    print(f"\nMSA depth: {depth} sequences  (alignment length: {len(wt_seq)} columns)")

    if depth < MIN_SEQS:
        print(f"WARNING: only {depth} sequences — conservation scores will be unreliable.")

    with open(OUTPUT_MSA, "w") as fh:
        for seq_id, seq in msa:
            fh.write(f">{seq_id}\n{seq}\n")

    print(f"Saved: {OUTPUT_MSA}")


if __name__ == "__main__":
    main()
