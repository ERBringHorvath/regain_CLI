#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 University of Utah

import os
import subprocess
import csv
import tempfile
import time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def make_blast_db(args, temp_dir):
    """
    Create required BLAST databases from input genome FASTA files
    """
    # Use the temp_dir for storing the output databases
    output_dir = temp_dir

    input_dir = args.fasta_directory
    extensions = [".fasta", ".fna", ".fa", ".fas", ".faa", ".ffn", ".frn"]

    print("\n \033[93mBuilding databases....\033[0m")

    for filename in os.listdir(input_dir):
        if any(filename.endswith(ext) for ext in extensions):
            input_file = os.path.join(input_dir, filename)
            base_name = os.path.splitext(filename)[0]
            output_file = os.path.join(output_dir, base_name)

            dbtype = 'nucl' if filename.endswith(('.fasta', '.fna', '.fa', '.fas', '.ffn', '.frn')) else 'prot'
            cmd = f"makeblastdb -in {input_file} -out {output_file} -dbtype {dbtype}"
            subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print("\n \033[92mDatabases built\033[0m \n")

def execute_blast_query(data):
    blast, db_path, query_file_path, output_file, evalue_threshold, db_name, query_file_base_name, min_seq_len = data
    cmd = f"{blast} -query {query_file_path} -db {db_path} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out {output_file} -evalue {evalue_threshold}"
    
    if min_seq_len:
        cmd += f" -task {blast}-short -dust no"
    subprocess.run(cmd, shell=True)

    return output_file, db_name, query_file_base_name 

def run_multiblastp(args, temp_dir, report_only_lowest_evalue):
    """
    Runs BLAST queries in parallel, combines results, then filters output based on internal thresholds
    or user-defined input.
    """

    blast = 'blastn' if args.nucleotide_query else 'tblastn'
    db_dir = temp_dir  # temp_dir is a string representing the directory path
    query_path = args.query
    threads = args.threads if hasattr(args, 'threads') else 4
    evalue_threshold = args.evalue if hasattr(args, 'evalue') and args.evalue is not None else 0.00001
    perc_identity_threshold = args.perc if hasattr(args, 'perc') and args.perc is not None else 90
    query_coverage_threshold = args.cov if hasattr(args, 'cov') and args.cov is not None else 75

    results_output_dir = os.path.join(os.getcwd(), 'ReGAIN_Curate')
    if not os.path.exists(results_output_dir):
        os.makedirs(results_output_dir)

    if not os.path.exists(query_path):
        print(f" \033[91m{query_path} does not appear to be a valid directory\033[0m \n")
        return

    extensions = ['.fasta', '.fna', '.fa', '.fas', ".faa"]
    fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen']

    print(" \033[093mSearching databases....\033[0m\n")

    # Define tasks for BLAST execution
    tasks = []
    for query_file in os.listdir(query_path):
        if any(query_file.endswith(ext) for ext in extensions):
            query_file_path = os.path.join(query_path, query_file)
            query_file_base_name = os.path.splitext(query_file)[0]

            for file_name in os.listdir(db_dir):
                if file_name.endswith(".nhr"):
                    db_name = os.path.splitext(file_name)[0]
                    db_path = os.path.join(db_dir, db_name)
                    output_file_path = os.path.join(results_output_dir, f"{db_name}_{query_file_base_name}_results.txt")
                    # Append db_name and query_file_base_name to the task for direct access later
                    tasks.append((blast, db_path, query_file_path, output_file_path, evalue_threshold, db_name, query_file_base_name, args.min_seq_len))

    # Execute BLAST queries in parallel
    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(execute_blast_query, tasks)

    # Combine all results into a single DataFrame
    combined_results = []
    for output_file in os.listdir(results_output_dir):
        if output_file.endswith("_results.txt"):
            output_file_path = os.path.join(results_output_dir, output_file)
            with open(output_file_path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) == len(fieldnames): 
                        task_info = next((task[5:7] for task in tasks if task[3] == output_file_path), (None, None))
                        combined_results.append(row + list(task_info))  #Append db_name and query_file_base_name

    ##Check if there are any results to process
    if not combined_results:
        print("\n \033[91mNo BLAST results found. Please check the input files and parameters\033[0m \n")
        return
    
    df = pd.DataFrame(combined_results, columns=fieldnames + ['database', 'query_file_name'])

    ##Convert necessary columns to numeric for calculations
    df[['qstart', 'qend', 'qlen', 'pident', 'evalue']] = df[['qstart', 'qend', 'qlen', 'pident', 'evalue']].apply(pd.to_numeric, errors='coerce')

    #Calculate query coverage
    df['query_coverage'] = ((df['qend'] - df['qstart']) / df['qlen'] * 100).round(2)

    #Save all results before filtering (for --report-all case)
    if report_only_lowest_evalue is False:
        all_results_csv_path = os.path.join(results_output_dir, "all_results.csv")
        df.to_csv(all_results_csv_path, index=False)
        print(f" \033[92mAll results saved in {all_results_csv_path}\033[0m")

    # Sorting and filtering
    df = df.sort_values(by=['database', 'query_file_name', 'pident', 'query_coverage', 'evalue'], ascending=[True, True, False, False, True])
    filtered_df = df[(df['evalue'] <= evalue_threshold) & (df['pident'] >= perc_identity_threshold) & (df['query_coverage'] >= query_coverage_threshold)]

    # print("Sample of DataFrame after applying filters:")
    # print(filtered_df.head())

    if report_only_lowest_evalue:
        strongest_hits = filtered_df.groupby(['database', 'query_file_name']).first().reset_index()
        output_csv_path = os.path.join(results_output_dir, "filtered_results.csv")
        strongest_hits.to_csv(output_csv_path, index=False)
        print(f" \033[092mStrongest matches stored in {output_csv_path}\033[0m\n")
    else:
        #strongest_hits = filtered_df
        output_csv_path = os.path.join(results_output_dir, "all_filtered_results.csv")
        filtered_df.to_csv(output_csv_path, index=False)
        print(f"\n \033[92mAll filtered results stored in {output_csv_path}\033[0m\n")

    return filtered_df

def simplify_gene_names(gene_name):
    if gene_name is None:
        return None
    gene_name = gene_name.replace("'", "p")
    gene_name = gene_name.replace('"', "pp")
    gene_name = gene_name.replace(".", "_")
    gene_name = gene_name.replace("(", "")
    gene_name = gene_name.replace(")", "")
    gene_name = gene_name.replace("-", "_")
    gene_name = gene_name.replace("/", "_")
    return gene_name

def create_matrix(filtered_df, output_filename='curate_matrix.csv', metadata_filename='curate_metadata.csv', keep_gene_names=False, min_val=None, max_val=None):
    """
    Generates a presence/absence matrix and metadata file from filtered dataframe for Bayesian network analysis
    """

    try:

        output_directory = os.path.join(os.getcwd())

        filtered_df = filtered_df.copy()

        if not keep_gene_names:
            filtered_df['query_file_name'] = filtered_df['query_file_name'].apply(simplify_gene_names)

        filtered_df = filtered_df.drop_duplicates(subset=['database', 'query_file_name'])
        filtered_df.loc[:, 'presence'] = 1

        matrix_path = os.path.join(output_directory, output_filename)
        metadata_path = os.path.join(output_directory, metadata_filename)

        matrix_df = filtered_df.pivot_table(index='database', columns='query_file_name', values='presence', fill_value=0)

        if min_val is not None and max_val is not None and min_val > max_val:
                    raise ValueError("Minimum gene occurrence value cannot be greater than maximum value.")

        if min_val is not None:
            matrix_df = matrix_df.loc[:, (matrix_df.sum(axis=0) >= min_val)].astype(int)
        if max_val is not None:
            matrix_df = matrix_df.loc[:, (matrix_df.sum(axis=0) <= max_val)].astype(int)

        matrix_df.to_csv(matrix_path)

        metadata_df = filtered_df[['query_file_name']].drop_duplicates()
        metadata_df.rename(columns={'query_file_name': 'Gene'}, inplace=True)
        metadata_df['GeneClass'] = ''
        metadata_df.to_csv(metadata_path, index=False)

    except Exception as e:
        print(f" \033[91mAn error occurred: {e}\033[0m\n")
        return None
    
    print(f" \033[92mData matrix file saved to {matrix_path}\033[0m\n")
    print(f" \033[92mMedata data file saved to {metadata_path}\033[0m\n")

    return matrix_df

def run(args):
    start_time = time.time()

    report_only_lowest_evalue = not args.report_all

    #Create a temporary directory for BLAST databases
    with tempfile.TemporaryDirectory() as temp_dir:
        #Pass the temporary directory object to the function
        make_blast_db(args, temp_dir)
        filtered_df = run_multiblastp(args, temp_dir, report_only_lowest_evalue)
        if filtered_df is not None and not filtered_df.empty:
            create_matrix(filtered_df, keep_gene_names=args.keep_gene_names, min_val=args.min, max_val=args.max)
        else:
            print("\033[91mNo data to process for the presence/absence matrix\033[0m \n")

    #Performance
    end_time = time.time()
    total_time = end_time - start_time
    minute_time = float((end_time - start_time) / 60)

    if total_time < 60:
        print(f" Total runtime: {total_time:.2f} seconds")
    else:
        print(f" Total runtime: {minute_time:.2f} minutes")