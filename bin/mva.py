"""
    ReGAIN Copyright 2024 University of Utah

Authors:
    Elijah R. Bring Horvath, Mathew G. Stein, Matthew A. Mulvey, E. Javier Hernandez, Jaclyn M. Winter

    Conducts Multidimensional Variance Analysis (MVA) on genomic data to explore gene distribution and relationships.
    
    This optional module supports various measures of distance to analyze the genomic data, providing insights
    into gene clustering and interaction. The analysis outputs visual representations in PNG and PDF formats,
    aiding in the interpretation of complex genomic data structures.
"""

import argparse
import subprocess
import os

def run_mva():
    parser = argparse.ArgumentParser(description="Perform Multivariate Analysis")
    parser.add_argument('-i', '--input', required=True, help='Input data file in CSV format')
    parser.add_argument('-m', '--method', default='euclidean', help='manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, mahalanobis, chisq, chord, hellinger, aitchison, or robust.aitchison')
    parser.add_argument('-c', '--num-centers', type=int, default=1, help='How many clusters to group into (max 10)')
    parser.add_argument('-C', '--confidence', type=float, default=0.95, help='Enter confidence value (suggested 0.95)')

    args = parser.parse_args()

    dir_of_executable = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executable, 'mva.sh')

    command = command = command = [path_to_shell_script, args.input, args.method, str(args.num-centers), str(args.confidence)]
    subprocess.run(command)

if __name__ == "__main__":
    run_mva()