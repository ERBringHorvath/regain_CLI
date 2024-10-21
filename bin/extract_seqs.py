import os
import pandas as pd
import time
from Bio import SeqIO
import concurrent.futures

"""
ReGAIN Extract is designed to extract sequences based on gene hits identified using the ReGAIN Curate Module
As extracted sequences may not represent whole genes, the --translate flag should be used with care,
as sequence will be trimmed to the closest length to afford codons
"""

#Process input FASTA files
def find_fasta_file(basename, fasta_dir):
    
    extensions = {'fa', 'fas', 'fasta', 'fna', 'faa'}
    for file in os.listdir(fasta_dir):
        name_part, extension = os.path.splitext(file)
        extension = extension.lstrip('.').lower()
        if extension in extensions:
            if name_part == basename or file.startswith(basename + '_') or file.startswith(basename + '.'):
                return os.path.join(fasta_dir, file)
    return None

#Process seqs for extraction
def process_sequence_entry(row, fasta_dir, translate, min_evalue, min_perc, min_cov):
    if float(row['evalue']) > min_evalue or float(row['pident']) < min_perc or float(row['query_coverage']) < min_cov:
        print(f"\n \033[91mSkipping sequence ({row['database']} query: {row['query_file_name']} pident: {row['pident']} query coverage: {row['query_coverage']} evalue: {row['evalue']}) due to filtering thresholds\033[0m")
        return None

    original_fasta = find_fasta_file(row['database'], fasta_dir)
    if original_fasta is None:
        print(f"\n \033[91mNo matching FASTA file found for {row['database']}\033[0m")
        return None
    
    found = False
    for seq_record in SeqIO.parse(original_fasta, "fasta"):
        if seq_record.id == row['sseqid']:
            found = True
            sstart, send = int(row['sstart']), int(row['send'])
            strand = 1 if sstart < send else -1
            if strand == 1:
                sequence = seq_record.seq[sstart-1:send]  # Adjust for 0-based index
            else:
                sequence = seq_record.seq[send-1:sstart].reverse_complement()  # Adjust for 0-based index and reverse complement

            if translate:
                sequence = sequence[:len(sequence) - len(sequence) % 3].translate()

            header_id = f"{seq_record.id}_{row['database']}_{row['query_file_name']}_aligned_region"
            description = "translated" if translate else "aligned"

            return SeqIO.SeqRecord(sequence, id=header_id, description=description)

    if not found:
        print(f"\n \033[91mSequence ID {row['sseqid']} not found in {original_fasta}\033[0m")
    return None

#Execute in parallel
def extract_sequences_from_csv(csv_path, fasta_dir, output_fasta, translate=False, min_evalue=1e-5, min_perc=90.0, min_cov=75.0):
    df = pd.read_csv(csv_path)
    sequences = []
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_sequence_entry, row, fasta_dir, translate, min_evalue, min_perc, min_cov) for index, row in df.iterrows()]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                sequences.append(result)

    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"\n \033[92mExtracted {len(sequences)} sequences to {output_fasta}\033[0m\n")

def run(args):

    start_time = time.time()

    extract_sequences_from_csv(
        csv_path=args.csv_path,
        fasta_dir=args.fasta_directory,
        output_fasta=args.output_fasta,
        translate=args.translate,
        min_evalue=args.min_evalue,
        min_perc=args.min_perc,
        min_cov=args.min_cov
    )

    end_time = time.time()
    print(f" Total runtime: {end_time-start_time:.2f} seconds")