"""
    ReGAIN Copyright 2024 University of Utah
    
    Performs Bayesian network structure learning on a given dataset to identify relationships between genes.
    
    This function applies either a small-scale (bnS) or large-scale (bnL) approach based on the dataset size,
    utilizing bootstrapping and resampling techniques to validate the network's structure. It's designed
    to handle genomic data, specifically presence/absence matrices, facilitating the understanding of complex
    gene interactions and resistance mechanisms.
    
    Parameters:
        input_file (str): CSV file containing the presence/absence matrix.
        metadata_file (str): File containing gene names and descriptions.
        output_bootstrap (str): Output file for bootstrap results.
        threads (int): Number of CPU cores to use.
        number_of_bootstraps (int): How many bootstrap iterations to perform.
        number_of_resamples (int): Number of resamples for bootstrapping.
    
    Returns:
        None
"""

import argparse
import subprocess
import os

def run_bnL():
    parser = argparse.ArgumentParser(description='Run Bayesian Network analysis')
    parser.add_argument('-i', '--input', required=True, help='Input file in CSV format')
    parser.add_argument('-M', '--metadata', required=True, help='Input metadata file with genes to query')
    parser.add_argument('-o', '--output_boot', required=True, help='Output file name for Network')
    parser.add_argument('-T', '--threads', type=int, help='Dedicated threads')
    parser.add_argument('-n', '--number_of_bootstraps', type=int, help='Bootstrap number (ideally 300-500)')
    parser.add_argument('-r', '--number_of_resamples', type=int, required=True, help='Input number of data resamples')

    args = parser.parse_args()

    # Get the directory of the current script
    dir_of_executable = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executable, 'bnCPQuery.sh')

    # Construct the command to call the shell script
    command = [path_to_shell_script, args.input, args.output_boot, str(args.threads), str(args.number_of_bootstraps), str(args.number_of_resamples), args.metadata]
    subprocess.run(command)

if __name__ == "__main__":
    run_bnL()