#!/usr/bin/env python3

import argparse, logging, sys, subprocess, os
from datetime import datetime
from . import __version__

def run_bnL(args):
    dir_of_executable = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executable, 'bnCPQuery.sh')
    command = [
        path_to_shell_script,
        "--input", args.input,
        "--metadata", args.metadata,
        "--output_boot", args.output_boot,
        "--threads", str(args.threads) if args.threads is not None else "2",
        "--number_of_bootstraps", str(args.number_of_bootstraps),
        "--resamples", str(args.number_of_resamples),
        "--iss", str(args.iss if args.iss is not None else 10),
        "--cp-samples", str(args.cp_samples if args.cp_samples is not None else 10000),

    ]
    if getattr(args, "blacklist", None):
        command += ["--blacklist", args.blacklist]

    if getattr(args, "no_viz", False):
        command += ["--no-viz"]

    subprocess.run(command, check=True)

def run_bnS(args):
    dir_of_executable = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executable, 'bnQueryGrain.sh')
    command = [
        path_to_shell_script,
        "--input", args.input,
        "--metadata", args.metadata,
        "--output_boot", args.output_boot,
        "--threads", str(args.threads) if args.threads is not None else "2",
        "--number_of_bootstraps", str(args.number_of_bootstraps),
        "--resamples", str(args.number_of_resamples),
        "--iss", str(args.iss if args.iss is not None else 10),
    ]
    if getattr(args, "blacklist", None):
        command += ["--blacklist", args.blacklist]
    
    if getattr(args, "no_viz", False):
        command += ["--no-viz"]

    subprocess.run(command, check=True)

def run_mva(args):
    dir_of_executable = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executable, 'mva.sh')
    command = [path_to_shell_script, args.input, args.method, str(args.num_centers), str(args.confidence)]
    subprocess.run(command)

def run_network(args):
    dir_of_executables = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executables, 'bayesnetwork.sh')
    command = [path_to_shell_script, args.input, args.data, args.metadata, args.statistics_results]
    subprocess.run(command)

def check_curate_range(value):
    ivalue = int(value)
    if ivalue < 0 or ivalue > 100:
        raise argparse.ArgumentTypeError(f" \033[91m{value} is out of allowed range (0-100)\033[0m \n")
    return ivalue

def check_extract_range(value):
    value = float(value)
    if value < 0 or value > 100:
        raise argparse.ArgumentTypeError(f" \033[91mValue is out of allowed range (0-100)\033[0m \n")
    return value

def main():
    parser = argparse.ArgumentParser(prog='regain',
                                     description='ReGAIN: Resistance Gene Association and Inference Network',
                                     epilog=" \033[92mFor detailed instructions, visit: https://github.com/ERBringHorvath/regain_cl\033[0m")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {__version__}')

    subparsers = parser.add_subparsers(dest='command')
    
    amr_parser = subparsers.add_parser('AMR', help='Module 1: Resistance/Virulence Gene Identification')
    amr_parser.add_argument('-d', '--directory', required=True, help='The directory containing your FASTA files')
    amr_parser.add_argument('-O', '--organism', help='Specify an organism for the AMRfinder search')
    amr_parser.add_argument('-T', '--threads', type=int, help='Dedicated threads')
    amr_parser.add_argument('-o', '--output-dir', help='Output directory path')

    matrix_parser = subparsers.add_parser('matrix', help='Module 1.1: Dataset Creation. Use AMR Module 1 results to prepare dataset for Bayesian network analysis Module 2')
    matrix_parser.add_argument('-d', '--directory', required=True, help='Input directory to CSV files to search')
    matrix_parser.add_argument('--keep-gene-names', action='store_true', help='Replace special characters in gene names when invoked')
    matrix_parser.add_argument('--min', type=int, default=5, required=True, help='Minimum required gene occurrences')
    matrix_parser.add_argument('--max', type=int, help='Maximum allowed gene occurrences')
    matrix_parser.add_argument('--gene-type', required=True, choices=['resistance', 'virulence', 'all'], help='specify gene type: resistance, virulence, all')
    matrix_parser.add_argument('--report-all', action='store_true', help='Return actual gene count. Overrides binary output.')
    
    bnS_parser = subparsers.add_parser('bnS', help='Module 2: Bayesian Network Structure Learning. For use with datasets containing < 100 genes')
    bnS_parser.add_argument('-i', '--input', required=True, help='Input file in CSV format')
    bnS_parser.add_argument('-M', '--metadata', required=True, help='Input metadata file with genes to query')
    bnS_parser.add_argument('-o', '--output_boot', required=True, help='Output file name for Network')
    bnS_parser.add_argument('-T', '--threads', type=int, help='Dedicated threads')
    bnS_parser.add_argument('-n', '--number-of-bootstraps', required=True, type=int, help='Bootstrap number (ideally 300-500)')
    bnS_parser.add_argument('-r', '--number-of-resamples', type=int, required=True, help='Input number of data resamples')
    bnS_parser.add_argument('-b', '--blacklist', help='Optional blacklist CSV (no header): from,to')
    bnS_parser.add_argument('--iss', type=int, default=10, help='Imaginary sample size for BDe score (default: 10)')
    bnS_parser.add_argument('--no-viz', dest='no_viz', action='store_true',
                        help='Skip HTML/PDF visualization')

    bnL_parser = subparsers.add_parser('bnL', help='Module 2: Bayesian Network Structure Learning. For use with datasets containing ≥ 100 genes')
    bnL_parser.add_argument('-i', '--input', required=True, help='Input file in CSV format')
    bnL_parser.add_argument('-M', '--metadata', required=True, help='Input metadata file with genes to query')
    bnL_parser.add_argument('-o', '--output_boot', required=True, help='Output file name for Network')
    bnL_parser.add_argument('-T', '--threads', type=int, help='Dedicated threads')
    bnL_parser.add_argument('-n', '--number-of-bootstraps', required=True, type=int, help='Bootstrap number (ideally 300-500)')
    bnL_parser.add_argument('-r', '--number-of-resamples', required=True, type=int, help='Input number of data resamples')
    bnL_parser.add_argument('-b', '--blacklist', help='Optional blacklist CSV (no header): from,to')
    bnL_parser.add_argument('--iss', type=int, default=10, help='Imaginary sample size for BDe score (default: 10)')
    bnL_parser.add_argument('--no-viz', dest='no_viz', action='store_true',
                        help='Skip HTML/PDF visualization')
    bnL_parser.add_argument('--cp-samples', dest='cp_samples', type=int, default=10000,
                        help='Monte Carlo samples for cpquery (default: 10000)')

    network_parser = subparsers.add_parser('network', help='Visualize Bayesian network. For use if visualization step in bnS/bnL fails')
    network_parser.add_argument('-i', '--input', required=True, help='Input Network .RDS File')
    network_parser.add_argument('-d', '--data', required=True, help='Input data matrix file')
    network_parser.add_argument('-M', '--metadata', required=True, help='Input metadata file')
    network_parser.add_argument('-s', '--statistics_results', required=True, help='Input results file from bnL/bnS analysis')

    curate_parser = subparsers.add_parser('curate', help='ReGAIN Curate: Create a curated dataset for Bayesian network analysis Module using your own genes of interest')
    curate_parser.add_argument('-d', '--directory', required=True, help='Input directory to genomes in FASTA format')
    curate_parser.add_argument('-q', '--query', required=True, help='Input directory to gene query files in FASTA format')
    curate_parser.add_argument('-T', '--threads', type=int, help='How many cores to dedicate')
    curate_parser.add_argument('--nucleotide-query', action='store_true', help='Use blastn for nucleotide queries')
    curate_parser.add_argument('--report-all', action='store_true', help='Report all BLAST hits instead of only the strongest hit')
    curate_parser.add_argument('--perc', type=check_curate_range, help='Override internal percent identity threshold (default = 75)')
    curate_parser.add_argument('--cov', type=check_curate_range, help='Override internal query coverage threshold (default = 80)')
    curate_parser.add_argument('--min-seq-len', type=int, help='Minimum sequence length for database searches')
    curate_parser.add_argument('--keep-gene-names', action='store_true', help='Do not remove special characters from gene names')
    curate_parser.add_argument('--min', type=int, default=5, required=True, help='Minimum required gene occurrence in filtered data matrix')
    curate_parser.add_argument('--max', type=int, required=True, help='Maximum allowed gene occurrence in filtered data matrix')

    parser_extract = subparsers.add_parser('extract', help="ReGAIN Curate: extract sequences based ReGAIN Curate results")
    parser_extract.add_argument('-c', '--csv-path', required=True, help="Path to BLAST results files.")
    parser_extract.add_argument('-f', '--fasta-directory', required=True, help="Path to reference FASTA assemblies.")
    parser_extract.add_argument('-o', '--output-fasta', required=True, help="Output FASTA file name.")
    parser_extract.add_argument('-T', '--threads', help="How many CPUs to dedicate.")
    parser_extract.add_argument('--min-evalue', type=float, default=1e-5, help='Minimum e-value threshold for sequence extraction.')
    parser_extract.add_argument('--min-perc', type=check_extract_range, default=90.0, help='Minimum percent identity threshold for sequence extraction.')
    parser_extract.add_argument('--min-cov', type=check_extract_range, default=75.0, help='Minimum query coverage threshold for sequence extraction.')
    parser_extract.add_argument('--translate', action='store_true', help='Translate extracted nucleotide sequence using the standard genetic code.')

    parser_combine = subparsers.add_parser('combine', help='Combine datasets from ReGAIN AMR and ReGAIN Curate')
    parser_combine.add_argument('--matrix1', required=True, help='Path to ReGAIN AMR data matrix file')
    parser_combine.add_argument('--matrix2', required=True, help='Path to ReGAIN Curate data matrix file')
    parser_combine.add_argument('--metadata1', required=True, help='Path to ReGAIN AMR metadata file')
    parser_combine.add_argument('--metadata2', required=True, help='Path to ReGAIN Curate metadata file')
    parser_combine.add_argument('--delete-duplicates', action='store_true', help='Delete any duplicate values from combined dataset')

    mva_parser = subparsers.add_parser('MVA', help='Multivariate analysis. For use with datasets from Module 1.1 or ReGAIN Curate')
    mva_parser.add_argument('-i', '--input', required=True, help='Input data file in CSV format')
    mva_parser.add_argument('-m', '--method', default='euclidean', help='manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, mahalanobis, chisq, chord, hellinger, aitchison, or robust.aitchison')
    mva_parser.add_argument('-c', '--num_centers', type=int, default=1, help='How many clusters to group into (max 10)')
    mva_parser.add_argument('-C', '--confidence', type=float, default=0.95, help='Enter confidence value')

    args = parser.parse_args()

    # Set up logging for our program
    logging.basicConfig(filename=f'{args.command}_log.txt', level=logging.INFO)
    logging.info(f'Program started at {datetime.now()}')
    logging.info(f'User command: {" ".join(sys.argv)}')

    if args.command == 'AMR':
        import AMR
        AMR.run(args)
    elif args.command == 'bnL':
        run_bnL(args)
    elif args.command == 'bnS':
        run_bnS(args)
    elif args.command == 'MVA':
        run_mva(args)
    elif args.command == 'matrix':
        import matrix
        matrix.run(args)
    elif args.command == 'curate':
        import curate
        curate.run(args)
    elif args.command == 'network':
        import bayesnetwork
        run_network(args)
    elif args.command == 'extract':
        import extract_seqs
        extract_seqs.run(args)
    elif args.command == 'combine':
        import combine
        combine.run(args)
    else:
        parser.print_help()
    
    logging.info(f'Program ended at {datetime.now()}')

if __name__ == "__main__":
    main()
