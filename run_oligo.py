"""
Runs oligo-analysis on all pro-Seq data.
"""
from multiprocessing import Pool
import subprocess


def run_analysis(args):
    k, data = args
    outfile = f'{data}_k{k}'
    result = subprocess.call(["python", "./src/oligo_analysis.py", "-f",
                             f"./data/{data}_shortened.fasta", '-k', str(
                                 k), '--noov', '--norc', '-top', '10000', '-log', 'critical',
                              '-out', f'./results/oligo/{outfile}.csv'],  shell=False)


def main():
    k = [5, 6, 7]
    data = ['human', 'chimp', 'monkey', 'mouse', 'rat']

    args = []
    for oligo_length in k:
        for datafile in data:
            args.append([oligo_length, datafile])

    print("Running oligo-analysis on data...")
    pool = Pool(None)
    r = pool.map_async(run_analysis, args)
    r.wait()
    print("Done.")


if __name__ == '__main__':
    main()
