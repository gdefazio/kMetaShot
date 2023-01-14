# generate the complement or reverse_complement of DNA string


def complement(dna):
    """Return complement's dna string."""
    dna = dna.upper()
    mode = type(dna)
    if mode is str:
        compl = {"A": "T",
                 "T": "A",
                 "C": "G",
                 "G": "C",
                 "N": "N",
                 "R": "Y",
                 "Y": "R",
                 "S": "S",
                 "W": "W",
                 "K": "M",
                 "M": "K",
                 "B": "V",
                 "D": "H",
                 "H": "D",
                 "V": "B"}
    elif mode == bytes:
        compl = {b"A": b"T",
                 65: b"T",
                 b"T": b"A",
                 84: b"A",
                 b"C": b"G",
                 67: b"G",
                 b"G": b"C",
                 71: b"C",
                 b"N": b"N",
                 78: b"N",
                 b"R": b"Y",
                 b"Y": b"R",
                 b"S": b"S",
                 b"W": b"W",
                 b"K": b"M",
                 b"M": b"K",
                 b"B": b"V",
                 b"D": b"H",
                 b"H": b"D",
                 b"V": b"B"}
    else:
        raise ValueError('Mode must be "s" for string or "b" for byte string')
    for nucleotide in dna:
        if nucleotide not in compl.keys():
            raise NameError("not nucleotide in dna string: %s" % nucleotide)
    comp = []
    for nt in dna:
        comp.append(compl[nt])
    return comp


def complrev(dna):
    "Return complement reverse of dna string."
    dna = complement(dna)
    dna.reverse()
    try:
        return "".join(dna)
    except TypeError:
        return b''.join(dna)
