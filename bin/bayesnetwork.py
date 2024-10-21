"""
Docstrings
"""

import argparse
import subprocess
import os

def run_bayesnetwork():
    parser = argparse.ArgumentParser(description='Visualize Bayesian Network')
    parser.add_argument('-i', '--input', required=True, help='Input Network .RDS File')
    parser.add_argument('-d', '--data', required=True, help='Input data matrix file')
    parser.add_argument('-M', '--metadata', required=True, help='Input metadata file')
    parser.add_arguemtn('-s', '--statistics_results', required=True, help='Input results file from bnL/bnS analysis')

    args = parser.parse_args()

    dir_of_executables = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executables, 'network.sh')

    command = [path_to_shell_script, args.input, args.data, args.metadata, args.statistics_results]
    subprocess.run(command)

if __name__ == '__main__':
    run_bayesnetwork