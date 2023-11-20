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
