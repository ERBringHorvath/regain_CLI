"""
    ReGAIN Copyright 2024 University of Utah
    
    Identifies resistance and virulence genes from FASTA files using the AMRfinderPlus tool.
    
    This function is part of the ReGAIN suite, designed to analyze genomic sequences to detect
    antimicrobial resistance (AMR) and virulence factors. Users can specify an organism to target
    the analysis or conduct a non-specific search across multiple organisms. The function supports
    multi-threading to expedite analysis.

    This module utilizes the NCBI AMRfinderPlus Software:
        Feldgarden M, Brover V, Gonzalez-Escalona N, et al. 
        AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. 
        Sci Rep. 2021;11(1):12728. Published 2021 Jun 16. doi:10.1038/s41598-021-91456-0
    
    Parameters:
        directory (str): Path to the directory containing FASTA files for analysis.
        organism (str, optional): Specific organism to analyze. If None, analysis is non-specific.
        threads (int, optional): Number of CPU cores to dedicate to the analysis. Default is 1.
        output_dir (str): Directory to store the AMRfinder results.
    
    Returns:
        None
"""

import os
import argparse
from tqdm import tqdm

###Expected time
print("\033[35m" + "\nIf working with a large number of FASTA files, this analysis may take hours to several days to complete, depending on your system" + "\033[0m")

###Include ARMfinder citation, green ANSI escape code
print("\033[92m" + "\nReGAIN utilizes AMRfinder Plus, so be sure to cite them! \n \nFeldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728. Published 2021 Jun 16. doi:10.1038/s41598-021-91456-0 \n" + "\033[0m")

def run(args):
    directory = args.directory
    fasta_extensions = [".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn"]
    organisms = {
        "Acinetobacter_baumannii": "Acinetobacter_baumannii",
        "Burkholderia_cepacia": "Burkholderia_cepacia",
        "Burkholderia_pseudomallei": "Burkholderia_pseudomallei",
        "Campylobacter": "Campylobacter",
        "Clostridioides_difficile": "Clostridioides_difficile",
        "Enterococcus_faecalis": "Enterococcus_faecalis",
        "Enterococcus_faecium": "Enterococcus_faecium",
        "Escherichia": "Escherichia",
        "Klebsiella": "Klebsiella",
        "Neisseria": "Neisseria",
        "Pseudomonas_aeruginosa": "Pseudomonas_aeruginosa",
        "Salmonella": "Salmonella",
        "Staphylococcus_aureus": "Staphylococcus_aureus",
        "Staphylococcus_pseudintermedius": "Staphylococcus_pseudintermedius",
        "Streptococcus_agalactiae": "Streptococcus_agalactiae",
        "Streptococcus_pneumoniae": "Streptococcus_pneumoniae",
        "Streptococcus_pyogenes": "Streptococcus_pyogenes",
        "Vibrio_cholerae": "Vibrio_cholerae"
    }

    ##Specific or general
    organism_flag = ""
    if args.organism is not None:  
        if args.organism in organisms:  
            organism_flag = f"-O {organisms[args.organism]}"
            print(f"\033[36m Starting {args.organism} analysis...\n" + "\033[0m")
        else:  # If provided but not valid
            print("\033[36m Sorry, not a valid organism. Hit CTRL + C to restart or continue without organism-specific analysis.\n" + "\033[0m")
    else:  # If the organism argument is not provided
        print("\033[36m Starting general analysis...\n" + "\033[0m")

    threads_flag = ""
    if args.threads:
        threads_flag = f"--threads {args.threads}"
        print(f"\033[36m Using {args.threads} cores.\n" + "\033[0m")
    else:
        print("\033[36m Using 1 core...\n" + "\033[0m")

    output_directory = os.path.join(directory, "AMRfinder_Results")
    if args.output_dir:
        output_directory = args.output_dir
        os.makedirs(output_directory, exist_ok=True)

    files = [filename for filename in os.listdir(directory) if any(filename.endswith(ext) for ext in fasta_extensions)]
    for filename in tqdm(files, desc="Processing files", unit="file"):
        output_file = os.path.join(output_directory, f"{os.path.splitext(filename)[0]}.amrfinder.csv")
        command = f"amrfinder --plus -n {os.path.join(directory, filename)} {organism_flag} {threads_flag} -o {output_file}"
        os.system(command)