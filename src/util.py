def nucleotide_freqs(sequences: list[str]):
    """
    Calculates nucleotide frequencies from DNA sequences.

    >>> nucleotide_freqs(['ACGT', 'TGCA'])
    >>> {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    """
    freqs = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0
    for seq in sequences:
        for base in seq:
            base = base.upper()
            if base in freqs:
                freqs[base] += 1
                total += 1
    return {base: freq / total for base, freq in freqs.items()}


def reverse_compliment(sequence: str):
    """
    Returns the reverse compliment of the sequence.

    >>> reverse_compliment('ATGGGTAGCG')
    >>> 'CGCTACCCAT'
    """
    COMPLIMENT = {'A': 'T', 'T': 'A',
                  'G': 'C', 'C': 'G',
                  'N': 'N'}
    compliment_strand = ''
    for base in sequence:
        compliment_strand += COMPLIMENT[base]

    return compliment_strand[::-1]
