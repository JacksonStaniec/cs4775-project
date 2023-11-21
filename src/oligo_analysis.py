"""
Oligo-analysis motif finder. Based on RSAT oligo-analysis script: https://github.com/rsa-tools/rsat-code/blob/master/perl-scripts/oligo-analysis

Usage:
>>> python oligo_analysis.py -f motif.fasta -k 10
"""
import logging
import argparse
import numpy as np
from tabulate import tabulate
import util
import csv


# fasta reader from class
def read_fasta(filename):
    ''' Reads in the sequences from the motif files.

    Arguments:
        filename: which filename to read sequences from
    Returns:
        output: list of sequences in motif file
    '''
    with open(filename, "r") as f:
        output = []
        s = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                # skip the line that begins with ">"
                if s == "":
                    continue
                output.append(s)
                s = ""
            # keep appending to the current sequence when the line doesn't begin
            # with ">"
            else:
                s += l.strip()
        output.append(s)
        return output


class Pattern():
    """Maintains counts of oligomer patterns and which sequences they occurred in."""

    def __init__(self):
        self.occ = 0  # total number of occurrences
        self.mseq = set()  # sequence ids that contained a pattern
        self.obs_freq = 0  # observed relative frequency
        self.exp_freq = 0  # expected relative frequency
        self.ms_freq = 0  # observed frequency in sequences
        self.occ_P = 0  # p-vale of occurences

        # indexes where overlap is forbidden
        self.forbidden: set[int] = set()

    def updated_occurred(self, index, k, no_overlaps):
        """Increments the number of occurences. Accounts for overlapping matches."""
        if no_overlaps:
            if index not in self.forbidden:
                self.occ += 1

                # update forbidden region
                for i in range(0, k):
                    self.forbidden.add(index - i)
        else:
            self.occ += 1

    def reset_forbidden(self):
        # indexes where overlap is forbidden
        self.forbidden: set[int] = set()

    def update_contained_with(self, sequence):
        """Updates the sequences that contained this pattern."""
        self.mseq.add(sequence)

    def set_obs_freq(self, obs_freq):
        """Sets the observed relative frequency of this pattern."""
        self.obs_freq = obs_freq

    def set_ms_freq(self, ms_freq):
        """Sets the observed sequence frequency of this pattern."""
        self.ms_freq = ms_freq

    def set_exp_freq(self, exp_freq):
        """Sets the expected relative frequency of this pattern."""
        self.exp_freq = exp_freq


class OligoPatterns():
    """Maintains oligo patterns through analysis."""

    def __init__(self, num_sequences: int):
        self.patterns: dict[str, Pattern] = dict()
        self.__num_sequences = num_sequences
        self.__nucleotide_freq = None

        self.HEADERS = ['pattern', 'mseq', 'occ', 'exp', 'p-value']

    def get_pattern(self, pattern: str) -> Pattern:
        """Gets the corresponding Pattern object if it exists, 
        or creates a new one if not encountered before."""
        pattern = pattern.upper()
        if pattern not in self.patterns:
            self.patterns[pattern] = Pattern()

        return self.patterns[pattern]

    def get_results(self, top=20):
        """Gets `top` counts ranked by signifigance."""
        logging.info(f'Getting results')

        # keep track of pattern and p-value for comparison
        top_patterns = [(pattern, self.patterns[pattern].occ_P)
                        for pattern in self.patterns.keys()]
        top_patterns.sort(key=lambda tup: tup[-1])
        top_patterns = top_patterns[:top]

        results = [
            [pattern,
             len(self.patterns[pattern].mseq),
             self.patterns[pattern].occ,
             self.patterns[pattern].exp_freq * self.total_occurences,
             self.patterns[pattern].occ_P] for pattern, _ in top_patterns
        ]

        return results

    def print_results(self, top=20):
        """Prints `top` counts ranked by signifigance."""
        logging.info(f'Printing results')

        print(tabulate(self.get_results(top=top), headers=self.HEADERS))

    def write_results(self, outfile, top=20):
        """Writes `top` results to `outfile`.csv"""
        logging.info(f'Writing results')
        with open(outfile, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(self.HEADERS)

            # write data
            csvwriter.writerows(self.get_results(top=top))

    def calibrate_nucleotide_freq(self, nucleotide_freq: dict[str, int]):
        """Sets the nucleotide frequency to use in significance calculations."""
        self.__nucleotide_freq = nucleotide_freq

        # ambiguous nucleotide code
        self.__nucleotide_freq['N'] = 1
        logging.info(f'Nucleotide frequencies used: {self.__nucleotide_freq}')

    def remove_repeat_patterns(self):
        """Removes all patterns from consideration that are just repeat nucleotides. (Ex. TTTTTT)"""
        # repeat if first character is all characters
        all_patterns = list(self.patterns.keys())
        for pattern in all_patterns:
            if (pattern.count(pattern[0]) == len(pattern)):
                self.patterns.pop(pattern)

    def reset_overlap_regions(self):
        """Resets overlapping regions in preparation for new sequence."""
        for pattern in self.patterns.values():
            pattern.reset_forbidden()

    @property
    def total_occurences(self):
        """The total number of occurrences across all patterns."""
        sum = 0
        for _, pattern in self.patterns.items():
            sum += pattern.occ
        return sum

    @property
    def num_sequences(self):
        """The number of sequences used in this analysis"""
        return self.__num_sequences

    @property
    def nucleotide_freq(self):
        """The nucleotide frequency used in this analysis: Of the form `{A: 0.89, G: 0.11, ...}`. """
        return self.__nucleotide_freq


def count_oligos(sequences: list[str], oligo_length: int, no_overlaps: bool):
    """Discover and count occurences for all k-mers."""
    logging.info('Counting oligo frequencies')

    oligo = OligoPatterns(len(sequences))
    for seq_id, sequence in enumerate(sequences):
        last_pos = len(sequence) - oligo_length

        # count oligomers
        current_pos = last_pos

        while (current_pos >= 0):
            # occurences
            pattern_seq = sequence[current_pos: current_pos + oligo_length]
            curr_pattern = oligo.get_pattern(pattern_seq)

            curr_pattern.updated_occurred(
                current_pos, oligo_length, no_overlaps)
            curr_pattern.update_contained_with(seq_id)

            current_pos -= 1

        oligo.reset_overlap_regions()

    logging.info(f'{len(oligo.patterns)} oligomers to consider')

    return oligo


def calc_frequencies(oligo: OligoPatterns):
    """Calculates relative frequencies of each oligomer."""
    logging.info("Calculating relative frequencies")

    # relative frequencies of each oligo pattern
    total_occ = oligo.total_occurences
    for pattern_seq in oligo.patterns.keys():
        pattern = oligo.get_pattern(pattern_seq)
        pattern.set_obs_freq(pattern.occ / total_occ)

    # frequency of matching sequences
    total_seq = oligo.num_sequences
    for pattern_seq in oligo.patterns.keys():
        pattern = oligo.get_pattern(pattern_seq)
        pattern.set_ms_freq(len(pattern.mseq) / total_seq)


def calc_expected(oligo: OligoPatterns):
    """Calculates expected frequencies of each oligomer."""
    logging.info("Calculating expected frequencies")

    for pattern_seq in oligo.patterns.keys():
        pattern = oligo.get_pattern(pattern_seq)
        exp_freq = 1
        for base in pattern_seq:
            exp_freq *= oligo.nucleotide_freq[base]

        pattern.set_exp_freq(exp_freq)


def sum_of_binomials(proba, trials, a, b):
    # copied from RSAT tools. https://github.com/rsa-tools/rsat-code/blob/d1aa8ad2964d1b3afd72e85aa609be217c204356/python-scripts/lib/dist.py#L114
    logproba = np.log(proba)
    q = max(0, 1 - proba)
    logq = np.log(q)
    sum_of_bin = 0
    logbin = 0
    prev_value = -1
    expected = trials * proba
    logbin = trials * np.log(q)
    if a == 0:
        sum_of_bin = np.exp(logbin)
        loop_start = 1
    else:
        loop_start = a
        for x in range(1, a):
            logbin += np.log(trials - x + 1) - np.log(x) + logproba - logq

    for x in range(loop_start, b+1):
        logbin += np.log(trials - x + 1) - np.log(x) + logproba - logq
        sum_of_bin += np.exp(logbin)

        if (x > expected) and (sum_of_bin <= prev_value):
            break
        prev_value = sum_of_bin

    return max(min(sum_of_bin, 1.0), 0.0)


def calc_prob(oligo: OligoPatterns):
    """Calculates right-tailed p-values."""
    logging.info("Calculating p-values")
    total_occ = oligo.total_occurences

    # report progress every 5%
    quantile = round(len(oligo.patterns.keys()) / 20)
    thresholds = [i * quantile for i in range(1, 20)]

    # right tailed (overepresentation)
    i = 0
    for pattern_seq in oligo.patterns.keys():
        if i in thresholds:
            logging.info(
                f'Calculating p-values - {(thresholds.index(i) + 1) * 5}% complete')

        pattern = oligo.get_pattern(pattern_seq)
        pattern.occ_P = sum_of_binomials(
            pattern.exp_freq, total_occ, pattern.occ, total_occ)
        i += 1

    logging.info(f'Calculating p-values - 100% complete')


def main():
    parser = argparse.ArgumentParser(
        description='Oligo-analysis')
    parser.add_argument('-f', action="store", dest="f",
                        type=str, default='error',
                        help="The fasta file with the sequences")
    parser.add_argument('-k', action="store", dest="k",
                        type=int, default=10, help="Find k-mer motifs")
    LOG_LEVELS = {'info': logging.INFO,
                  'debug': logging.DEBUG, 'critical': logging.CRITICAL}
    parser.add_argument('-log', action="store", dest="loglevel", choices=list(LOG_LEVELS.keys()),
                        type=str, default="info", help="Display log messages during execution")
    parser.add_argument('-top', action="store", dest="top",
                        type=int, default=25, help="Number of results to display")
    parser.add_argument('-noov', action="store_true", dest="no_overlap",
                        default=False, help="Don't count overlapping matches")
    parser.add_argument('-print', action="store_true", dest="print_results",
                        default=False, help="Whether to print the results.")
    parser.add_argument('-out', action="store", dest="outfile",
                        default=None, help="File to write csv of results.")

    args = parser.parse_args()
    logging.basicConfig(encoding='utf-8',
                        level=LOG_LEVELS[args.loglevel],
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    sequences = read_fasta(args.f)
    k = args.k

    patterns = count_oligos(sequences, k, args.no_overlap)
    patterns.calibrate_nucleotide_freq(util.nucleotide_freqs(sequences))
    patterns.remove_repeat_patterns()
    calc_frequencies(patterns)
    calc_expected(patterns)
    calc_prob(patterns)

    if args.print_results:
        patterns.print_results(top=args.top)

    if args.outfile:
        patterns.write_results(args.outfile, top=args.top)


if __name__ == '__main__':
    main()
