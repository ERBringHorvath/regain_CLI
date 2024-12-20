import os
import glob
import csv
import pandas as pd

#Remove special characters for Bayesian network structure learning ReGAIN Module 2
def simplify_gene_names(gene_name):
    if gene_name is None:
        return None
    gene_name = gene_name.replace("'", "p")  # Replace single quote with p
    gene_name = gene_name.replace('"', "pp")  # Replace double quote with p
    gene_name = gene_name.replace(".", "_")  # Replace dot with underscore
    gene_name = gene_name.replace("(", "")  # Remove left parenthesis
    gene_name = gene_name.replace(")", "")  # Remove right parenthesis
    gene_name = gene_name.replace("-", "_")  # Replace hyphen with underscore
    gene_name = gene_name.replace("/", "_")  # Replace forward slash with underscore
    return gene_name

#Check tab or comma-separated input files
def detect_delimiter(file_path):
    with open(file_path, 'r', newline='') as file:
        first_line = file.readline()
        comma_count = first_line.count(',')
        tab_count = first_line.count('\t')
    return '\t' if tab_count > comma_count else ','

#Combine AMR output files without filtering
def combine_amrfinder_output_files(directory, output_file_path):
    all_filenames = glob.glob(os.path.join(directory, '*.csv'))

    if not all_filenames:
        raise FileNotFoundError(f"\n \033[91mNo CSV files found in the directory {directory}\033[0m \n")
    
    dfs = []
    for file in all_filenames:
        try:
            delimiter = detect_delimiter(file)
            df = pd.read_csv(file, sep=delimiter, na_filter=False)
            dfs.append(df)
        except Exception as e:
            print(f"\n \033[91mSkipping {file} due to an error: {e}\033[0m\n")
    
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv(output_file_path, index=False, sep='\t')
        print(f"\n \033[92mCombined results saved to {output_file_path}\033[0m\n")
    else:
        raise ValueError("\n \033[91mNo valid CSV files to combine\033[0m\n")

# Main function
def run(args):
    path = args.directory
    output_directory = os.path.join(os.getcwd(), 'ReGAIN_Dataset')
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    filtered_matrix_path = os.path.join(output_directory, 'filtered_matrix.csv')
    metadata_path = os.path.join(output_directory, 'metadata.csv')
    unfiltered_combined_path = os.path.join(output_directory, 'combined_AMR_results_unfiltered.csv')

    #Combine AMR output files without filtering and save unfiltered version
    combine_amrfinder_output_files(path, unfiltered_combined_path)
    combined_csv = pd.read_csv(unfiltered_combined_path, sep='\t')

    #Filter combined data based on gene type
    gene_type = args.gene_type.lower()
    if 'Subtype' in combined_csv.columns:
        if gene_type == 'resistance':
            combined_csv = combined_csv[(combined_csv['Subtype'].isin(['AMR', 'METAL', 'BIOCIDE', 'POINT']))]
        elif gene_type == 'virulence':
            combined_csv = combined_csv[(combined_csv['Subtype'].isin(['VIRULENCE', 'HEAT', 'ACID']))]
        elif gene_type == 'all':
            combined_csv = combined_csv[(combined_csv['Subtype'].isin(['AMR', 'METAL', 'BIOCIDE', 'POINT', 'VIRULENCE', 'HEAT', 'ACID']))]
    else:
        raise ValueError("\n \033[91mError with CSV file header\033[0m\n")
    
    combined_csv = combined_csv.drop_duplicates(subset=['Element symbol'])

    #Create metadata file
    metadata_df = combined_csv[['Element symbol', 'Class', 'Subclass']].copy()
    metadata_df.columns = ['Gene', 'GeneClass', 'GeneSubClass']
    metadata_df['GeneClass'] = metadata_df['GeneClass'].fillna('virulence').astype(str).str.title()
    metadata_df['GeneSubClass'] = metadata_df['GeneSubClass'].fillna('virulence').astype(str).str.title()

    #Create data matrix using unmodified gene names
    csv_files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.csv')]
    search_strings = metadata_df['Gene'].tolist()

    with open(filtered_matrix_path, 'w', newline='') as output_file:
        writer = csv.writer(output_file)
        headers = ['file'] + search_strings
        writer.writerow(headers)

        for f in csv_files:
            row = [os.path.basename(f).replace('.amrfinder.csv', '')]
            counts = [1 if search_str in open(f).read() else 0 for search_str in search_strings]
            row.extend(counts)
            writer.writerow(row)

    df_filtered = pd.read_csv(filtered_matrix_path)

    #Apply min/max thresholds to filter columns
    required_min = args.min
    required_max = args.max

    sums_df = pd.DataFrame(columns=['variable', 'sum'])
    sums_df['variable'] = df_filtered.columns[1:]

    for col in df_filtered.columns[1:]:
        sums_df.loc[sums_df['variable'] == col, 'sum'] = df_filtered[col].sum()

    cols_to_keep = list(sums_df[(sums_df['sum'] >= required_min) & (sums_df['sum'] <= required_max)]['variable'])
    cols_to_keep.insert(0, 'file')
    df_filtered = df_filtered[cols_to_keep]

    #Apply simplify_gene_names to metadata and matrix files unless --keep-gene-names is specified
    if not args.keep_gene_names:
        metadata_df['Gene'] = metadata_df['Gene'].apply(simplify_gene_names)
        df_filtered.columns = [simplify_gene_names(col) if col != 'file' else col for col in df_filtered.columns]
    
    metadata_df.to_csv(metadata_path, index=False)
    df_filtered.to_csv(filtered_matrix_path, index=False)

    #Generate unfiltered data matrix if --report-all flag is used
    if args.report_all:
        unfiltered_matrix_path = os.path.join(output_directory, 'unfiltered_matrix.csv')
        with open(unfiltered_matrix_path, 'w', newline='') as output_file:
            writer = csv.writer(output_file)
            writer.writerow(headers)
            for f in csv_files:
                row = [os.path.basename(f).replace('.amrfinder.csv', '')]
                counts = [open(f).read().count(search_str) for search_str in search_strings]
                row.extend(counts)
                writer.writerow(row)
        print(f" \033[92mUnfiltered matrix saved to {unfiltered_matrix_path}\033[0m\n")
    
    print(f" \033[92mFiltered matrix saved to {filtered_matrix_path}\033[0m\n")
    print(f" \033[92mMetadata file saved to {metadata_path}\033[0m\n")
