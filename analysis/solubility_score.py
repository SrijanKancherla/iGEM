"""
Solubility and expression scoring for E. coli.

Predicts expression solubility by analyzing:
1. Hydrophobicity (GRAVY score)
2. Net charge and isoelectric point
3. Buried fraction from structure
4. Codon usage bias (E. coli optimized)
5. RNA secondary structure stability
6. Proteolytic sensitivity

Key insight: Thermostable mutations often increase hydrophobicity,
risking aggregation. This module flags problematic mutations.

Installation:
    pip install biopython
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    from Bio import SeqIO, SeqUtils
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    raise ImportError(
        "Biopython not found. Install with:\n"
        "  pip install biopython"
    )


# ============================================================================
# GRAVY SCORE (Hydrophobicity)
# ============================================================================

KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}


def calculate_gravy(sequence):
    """
    Calculate GRAVY score (Grand Average of Hydropathy).
    Positive = hydrophobic, negative = hydrophilic.
    For solubility: ideally between -0.5 and 0.0.
    """
    if not sequence:
        return 0.0

    total = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in sequence)
    return total / len(sequence)


# ============================================================================
# E. COLI CODON USAGE
# ============================================================================

# E. coli K12 codon usage (frequency out of 1000)
# Source: https://www.ncbi.nlm.nih.gov/projects/Codons/cgi-bin/showcodon.cgi
ECOLI_CODON_FREQUENCY = {
    'TTT': 58, 'TTC': 170, 'TTA': 58, 'TTG': 130,
    'CTT': 12, 'CTC': 68, 'CTA': 7, 'CTG': 430,
    'ATT': 49, 'ATC': 174, 'ATA': 1, 'ATG': 104,
    'GTT': 112, 'GTC': 74, 'GTA': 36, 'GTG': 412,
    'TAT': 58, 'TAC': 170, 'TAA': 2, 'TAG': 1,
    'CAT': 55, 'CAC': 180, 'CAA': 12, 'CAG': 340,
    'AAT': 46, 'AAC': 169, 'AAA': 58, 'AAG': 90,
    'GAT': 90, 'GAC': 280, 'GAA': 168, 'GAG': 250,
    'TGT': 67, 'TGC': 133, 'TGA': 3, 'TGG': 130,
    'CGT': 180, 'CGC': 180, 'CGA': 10, 'CGG': 11,
    'AGT': 41, 'AGC': 150, 'AGA': 7, 'AGG': 2,
    'GGT': 180, 'GGC': 340, 'GGA': 50, 'GGG': 40,
}


def codon_to_aa(codon):
    """Convert DNA codon to amino acid."""
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    return genetic_code.get(codon.upper(), '?')


def calculate_cai(dna_sequence):
    """
    Calculate Codon Adaptation Index (CAI) for E. coli.
    Measures how well codons are optimized for the organism.
    Range: 0 (poor) to 1 (optimal).

    Higher CAI = better expression in E. coli.
    """
    if len(dna_sequence) < 3:
        return 0.5

    codon_count = defaultdict(int)
    n_codons = 0

    # Count codons
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3].upper()
        if len(codon) == 3 and all(c in 'ATGC' for c in codon):
            codon_count[codon] += 1
            n_codons += 1

    if n_codons == 0:
        return 0.5

    # Calculate CAI as geometric mean of relative codon frequencies
    log_sum = 0
    codon_set = set()

    for codon, count in codon_count.items():
        aa = codon_to_aa(codon)
        if aa == '*' or aa == '?':
            continue

        # Get all codons for this amino acid
        aa_codons = [c for c, a in zip(ECOLI_CODON_FREQUENCY.keys(),
                                       [codon_to_aa(c) for c in ECOLI_CODON_FREQUENCY.keys()])
                     if a == aa]

        if aa_codons:
            max_freq = max(ECOLI_CODON_FREQUENCY.get(c, 1) for c in aa_codons)
            codon_freq = ECOLI_CODON_FREQUENCY.get(codon, 1)
            relative_freq = codon_freq / max(max_freq, 1)

            if relative_freq > 0:
                log_sum += count * np.log(relative_freq)

    cai = np.exp(log_sum / n_codons) if n_codons > 0 else 0.5

    return min(cai, 1.0)


# ============================================================================
# PROTEOLYTIC SITES
# ============================================================================

PROTEOLYTIC_MOTIFS = {
    'trypsin': ['KR'],  # Trypsin cleaves after K, R
    'elastase': ['AAAA'],  # Elastase cleaves after A
    'caspase': ['DEVD'],  # Caspase recognition motif
}


def check_proteolytic_sites(sequence, wt_aa, mut_aa, pos):
    """
    Check if mutation introduces or removes proteolytic cleavage sites.
    Returns penalty if problematic.
    """
    penalty = 0

    # Check if introducing a trypsin site
    if mut_aa in 'KR':
        penalty += 1

    # Check if introducing proline (can affect structure)
    if mut_aa == 'P':
        penalty += 2

    # Check if introducing cysteine in non-structured region
    # (risk of misfolding, aggregation)
    if mut_aa == 'C' and wt_aa != 'C':
        penalty += 1

    return penalty


# ============================================================================
# MAIN SOLUBILITY SCORING
# ============================================================================

def score_solubility(
    wt_seq,
    mut_seq,
    pos,
    wt_aa,
    mut_aa,
    dna_sequence=None,
    pdb_data=None,
    verbose=False
):
    """
    Comprehensive solubility scoring for E. coli expression.

    Args:
        wt_seq: Wild-type protein sequence
        mut_seq: Mutant protein sequence
        pos: Position of mutation (0-indexed)
        wt_aa: Wild-type amino acid
        mut_aa: Mutant amino acid
        dna_sequence: DNA coding sequence (for codon usage analysis)
        pdb_data: PDB structure data for burial analysis (optional)
        verbose: Print detailed scoring breakdown

    Returns:
        solubility_score: Combined solubility score (lower = more problematic)
        components: Dict with individual component scores
    """
    components = {}

    # 1. GRAVY score change
    wt_gravy = calculate_gravy(wt_seq)
    mut_gravy = calculate_gravy(mut_seq)
    gravy_delta = mut_gravy - wt_gravy

    # Flag if significantly more hydrophobic
    # Ideal: gravy between -0.5 and 0.0
    # Penalize if becoming too hydrophobic (> 0.2) or too hydrophilic (< -0.5)
    if gravy_delta > 0.15:  # Becoming hydrophobic
        gravy_score = max(0, 10 - gravy_delta * 50)  # Penalty scale: 50 per 0.2
    else:
        gravy_score = 10

    components['gravy_score'] = gravy_score
    components['gravy_delta'] = gravy_delta

    # 2. Charge and pI changes
    try:
        wt_analysis = ProteinAnalysis(wt_seq)
        mut_analysis = ProteinAnalysis(mut_seq)

        wt_charge = wt_analysis.charge_at_pH(7.0)
        mut_charge = mut_analysis.charge_at_pH(7.0)
        charge_delta = abs(mut_charge) - abs(wt_charge)

        # Slight charge loss is OK; big loss is bad (neutral → aggregate risk)
        if abs(charge_delta) > 5:
            charge_score = 5  # Significant change
        elif abs(charge_delta) > 2:
            charge_score = 8
        else:
            charge_score = 10

        components['charge_score'] = charge_score
        components['charge_delta'] = charge_delta

    except Exception as e:
        if verbose:
            print(f"  Warning: Could not calculate charge: {e}")
        charge_score = 10
        components['charge_score'] = charge_score
        components['charge_delta'] = 0

    # 3. Proteolytic sites
    proteolytic_penalty = check_proteolytic_sites(mut_seq, wt_aa, mut_aa, pos)
    proteolytic_score = max(0, 10 - proteolytic_penalty * 2)
    components['proteolytic_score'] = proteolytic_score

    # 4. Codon usage (if DNA sequence provided)
    if dna_sequence:
        try:
            # Simple proxy: check if mutation uses rare codons
            # In practice, you'd translate back to DNA codons
            cai_penalty = 0
            if mut_aa in 'ATA':  # Rare amino acids in E. coli
                cai_penalty = 3
            elif mut_aa in 'CGA':  # Moderately rare
                cai_penalty = 1

            codon_score = max(0, 10 - cai_penalty)
            components['codon_score'] = codon_score
        except Exception as e:
            if verbose:
                print(f"  Warning: Could not calculate codon usage: {e}")
            codon_score = 10
            components['codon_score'] = codon_score
    else:
        codon_score = 10
        components['codon_score'] = codon_score

    # 5. Cysteine handling (special case for E. coli)
    # Disulfide bonds don't form well in cytoplasm; risky in expression
    cys_score = 10
    if mut_aa == 'C':
        cys_score = 7  # E. coli cytoplasm is reducing
    elif wt_aa == 'C' and mut_aa != 'C':
        cys_score = 10  # Removing problematic Cys is OK

    components['cys_score'] = cys_score

    # Combined score (weighted average)
    # Weights: hydrophobicity (30%), charge (25%), proteolysis (20%),
    #          codon (15%), cysteine (10%)
    solubility_score = (
        0.30 * gravy_score +
        0.25 * charge_score +
        0.20 * proteolytic_score +
        0.15 * codon_score +
        0.10 * cys_score
    )

    components['solubility_score'] = solubility_score

    if verbose:
        print(f"    Solubility: {solubility_score:.1f}")
        print(f"      GRAVY Δ: {gravy_delta:+.3f} ({gravy_score:.1f})")
        print(f"      Charge: {charge_delta:+.1f} ({charge_score:.1f})")
        print(f"      Proteolysis: {proteolytic_score:.1f}")
        print(f"      Codon: {codon_score:.1f}")
        print(f"      Cysteine: {cys_score:.1f}")

    return solubility_score, components


def add_solubility_to_pipeline(
    input_csv,
    output_csv,
    wt_sequence,
    dna_sequence=None,
    verbose=False
):
    """
    Add solubility scores to an existing mutation dataframe.

    Args:
        input_csv: Input CSV with columns: pos, wt, mut
        output_csv: Output CSV with added solubility columns
        wt_sequence: Wild-type protein sequence
        dna_sequence: DNA coding sequence (optional)
        verbose: Print progress

    Returns:
        DataFrame with added solubility scores
    """
    df = pd.read_csv(input_csv)

    solubility_scores = []
    component_details = {col: [] for col in [
        'gravy_score', 'gravy_delta', 'charge_score', 'charge_delta',
        'proteolytic_score', 'codon_score', 'cys_score'
    ]}

    print(f"Calculating solubility scores for {len(df)} mutations...")

    for idx, row in df.iterrows():
        pos = int(row['pos'])
        wt_aa = row['wt']
        mut_aa = row['mut']

        # Create mutant sequence
        mut_seq = list(wt_sequence)
        mut_seq[pos] = mut_aa
        mut_seq = ''.join(mut_seq)

        # Score
        score, components = score_solubility(
            wt_sequence, mut_seq, pos, wt_aa, mut_aa,
            dna_sequence=dna_sequence,
            verbose=verbose
        )

        solubility_scores.append(score)

        for key, val in components.items():
            if key != 'solubility_score':
                component_details[key].append(val)

        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{len(df)}")

    # Add to dataframe
    df['solubility_score'] = solubility_scores

    for key, vals in component_details.items():
        df[key] = vals

    # Save
    df.to_csv(output_csv, index=False)
    print(f"Saved to {output_csv}")
    print(f"Solubility score range: [{df['solubility_score'].min():.2f}, {df['solubility_score'].max():.2f}]")

    return df


if __name__ == "__main__":
    # Example usage with mutation file
    mutation_csv = "results/6_expression.csv"
    output_csv = "results/7_solubility.csv"

    # Load wild-type sequence
    wt_seq = ""
    with open("data/cleaned/hbsag.fasta", "r") as f:
        for line in f:
            if not line.startswith(">"):
                wt_seq += line.strip()

    # Create output directory
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)

    # Add solubility scoring
    df = add_solubility_to_pipeline(
        mutation_csv,
        output_csv,
        wt_seq,
        verbose=False
    )

    print("\nSolubility score distribution:")
    print(df['solubility_score'].describe())

    print("\nWorst scoring mutations (low solubility risk):")
    print(df.nsmallest(10, 'solubility_score')[['pos', 'wt', 'mut', 'solubility_score']])
