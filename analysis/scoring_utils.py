# Known HBsAg antigenic region

MHR_START = 99
MHR_END = 169

MAJOR_EPITOPE_START = 124
MAJOR_EPITOPE_END = 147


def tcell_penalty(pos):

    if MAJOR_EPITOPE_START <= pos <= MAJOR_EPITOPE_END:
        return 25

    return 0


def bcell_penalty(pos):

    if MHR_START <= pos <= MHR_END:
        return 10

    return 0


def expression_penalty(wt, mut):

    penalty = 0

    if wt == "C" or mut == "C":
        penalty += 10

    if mut == "P":
        penalty += 5

    return penalty


def conservation_penalty(pos):

    # placeholder until real MSA

    highly_conserved = {
        48,49,50,
        121,122,123,
        124,125,126,
        145
    }

    return 20 if pos in highly_conserved else 0