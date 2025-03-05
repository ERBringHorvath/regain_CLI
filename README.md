# <ins>**ReGAIN Installation and User guide**</ins> 

<img src="https://github.com/user-attachments/assets/149bad0c-fb6b-44c7-a8c6-a131eb979172" alt="image" width="363" height="418"/>

_________________________________________________________________________________

**Prerequisites**

Ensure that you have the following prerequisites installed on your system:

Python (version 3.8 or higher)

R (version 4 or higher)

NCBI AMRfinderPlus version 4.0.3 <br />
NCBI BLAST+ (Included in AMRfinderPlus installation)

[Install R](https://www.r-project.org/)

_________________________________________________________________________________


**We suggest that ReGAIN and all prerequisites are installed within a Conda environment**

Download [miniforge](https://github.com/conda-forge/miniforge/)

Create Conda environment and install [NCBI AMRfinderPlus](https://github.com/ncbi/amr/wiki/Install-with-bioconda)

`conda create -n regain python=3.10`

`conda activate regain`

Install AMRfinderPlus

`conda install -y -c conda-forge -c bioconda ncbi-amrfinderplus`

Check installation

`amrfinder -h`

Download ARMfinderPlus Database

`amrfinder -u`

Download ReGAIN to preferred directory

`git clone https://github.com/ERBringHorvath/regain_CLI`

Install Python dependencies

`pip install -r requirements.txt` or `pip3 install -r requirements.txt`

Add ReGAIN to your PATH

Add this line to the end of `.bash_profile`/`.bashrc` (Linux) or `.zshrc` (macOS):

`export PATH="$PATH:/path/to/regain_CLI/bin"`

Replace `/path/to/regain_CLI/bin` with the actual path to the directory containing the executable. <br />
Whatever the initial directory, this path should end with `/regain_CLI/bin`

Save the file and restart your terminal or run `source ~/.bash_profile` or `source ~/.zshrc`

Verify installation:

`regain --version`

use `-h`, `--help`, to bring up the help menu

`regain --help`

NOTE: ReGAIN utilizes shell scripts to execute some modules. You may need to modify your permissions <br />
to execute these scripts. If you run `regain --version` and see `permission denied: regain`, Navigate to <br />
`regain/bin`, then run both `chmod +x regain` and `chmod +x *.sh` and rerun `regain --version` <br />
and you should see something similar to: `regain v.1.5.0`

_________________________________________________________________________________

# <ins>**Programs and Example Usage**</ins> 

## **Resistance and Virulence Gene Identification** 

Module 1.0 `regain AMR`

`-d`, `--directory`, path to directory containing genome FASTA files to analyze <br />
`-O`, `--organism`, optional; specify what organism (if any) you want to analyze (optional flag) <br />
`-T`, `--threads`, number of cores to dedicate for parallel processing <br />
`-o`, `--output-dir`, output directory to store AMRfinder results

**Currently supported organisms and how they should be called:**

`Acinetobacter_baumannii` <br />
`Burkholderia_cepacia`<br />
`Burkholderia_pseudomallei`<br />
`Citrobacter_freundii`<br />
`Corynebacterium_diphtheriae`<br />
`Campylobacter`<br />
`Clostridioides_difficile`<br />
`Enterobacter_cloacae`<br />
`Enterobacter_asburiae`<br />
`Enterococcus_faecalis`<br />
`Enterococcus_faecium`<br />
`Escherichia`<br />
`Haemophilus_influenzae`<br />
`Klebsiella_oxytoca`<br />
`Klebsiella_pneumoniae`<br />
`Neisseria_gonorrhoeae`<br />
`Neisseria_meningitidis`<br />
`Pseudomonas_aeruginosa`<br />
`Salmonella`<br />
`Serratia_marcescens`<br />
`Staphylococcus_aureus`<br />
`Staphylococcus_pseudintermedius`<br />
`Streptococcus_agalactiae`<br />
`Streptococcus_pneumoniae`<br />
`Streptococcus_pyogenes`<br />
`Vibrio_cholerae`<br />
`Vibrio_vulnificus`<br />
`Vibrio_parahaemolyticus`<br />

**Module 1 example usage:**

Organism specific: <br />
`regain AMR -d path/to/FASTA/files -O Pseudomonas_aruginosa -T 8 -o path/to/output/directory`

Organism non-specific: <br />
`regain AMR -d path/to/FASTA/files -T 8 -o path/to/output/directory`

**Output files:**

One results file per submitted genome

_________________________________________________________________________________

## **Dataset Creation** 

**NOTE**: variable names <ins>**cannot**</ins> contain special characters–this transformation is automated during dataset creation <br />

Module 1.1 `regain matrix`
                                       
`-d`, `--directory`, path to AMRfinder results in CSV format <br />
`--gene-type`, searches for `resistance`, `virulence`, or `all` genes <br />
`--min`, minimum gene occurrence cutoff <br />
`--max`, maximum gene occurrence cutoff (should be less than number of genomes, see NOTE below) <br />
`--report-all`, optional; reports all genes identified, regardless of `--min`/`--max` threshold <br />
`--keep-gene-names`, optional; maintains special characters in variable names. Should **not** be used if proceeding to Module 2

**Module 1.1 example usage**

NOTE: Discrete Bayesian network anlyses requires all variables to exist in at least two states. For ReGAIN, these two states are 'present' and 'absent'. Ubiquitously occurring genes will break the analysis. 
Best practice is for *N* genomes, `--max` should MINIMALLY be defined as *N* - 1. Keep in mind that removing very low and very high abundance genes can reduce noise in the network.
                                            
`regain matrix -d path/to/AMRfinder/results --gene-type resistance --min 5 --max 475`

**NOTE: all results are saved in the 'ReGAIN_Dataset' folder, which will be generated within the directory <br /> 
defined by** `-d`/`--directory`

**Output files:**

`filtered_matrix.csv`: presence/absence matrix of genes <br />
`metadata.csv`: file containing genes identified in AMRfinderPlus analysis <br />
`combined_AMR_results_unfiltered.csv`: concatenated file of all AMRfinder/Plus results; this file contains contig and nucleotide location of all identified genes <br />
**If** `--report-all` **is used:** <br />
`unfiltered_matrix.csv`: presence/absence matrix of all genes identified, regardless of `--min`/`--max` thresholds

_________________________________________________________________________________

## **Bayesian Network Structure Learning** 

Module 2 `regain bnL` or `regain bnS`
                                            
`-i`, `--input`, input file in CSV format <br />
`-M`, `--metadata`, file containing gene names and descriptions <br />
`-o`, `--output_boot`, output bootstrap file <br />
`-T`, `--threads`, number of cores to dedicate for parallel processing <br />
`-n`, `--number_of_boostraps`, how many bootstraps to run (suggested 300-500) <br />
`-r`, `--number-of-resamples`, how many data resamples you want to use (suggested 100) <br />

**Module 2 example usage:**

**NOTE: We suggest using between 300 and 500 bootstraps and 100 resamples**

`bnS`, Bayesian network structure learning analysis for less than 100 genes <br />
`bnL`, Bayesian network structure learning analysis for 100 genes or greater

For less than 100 genes:

`regain bnS -i matrix_filtered.csv -M metadata.csv -o bootstrapped_network -T 8 -n 500 -r 100`
                                            
For 100 or more genes:

`regain bnL -i matrix_filtered.csv -M metadata.csv -o bootstrapped_network -T 8 -n 500 -r 100`

**Output files:**

`Results.csv`, results file of all conditional probability and relative risk values <br />
`post_hoc_analysis.csv`, results file of all bidirectional probability and fold change scores <br />
`Bayesian_Network.html`, interactive Bayesian network

_________________________________________________________________________________
# **ReGAIN Curate**

ReGAIN Curate is designed to allow users to generate a dataset for Bayesian network structure learning using <br />
a custom set of gene queries, independent of ReGAIN Module 1

`regain curate`

`-d`, `--directory`, path to genome FASTA files <br />
`-q`, `--query`, path to query files containing amino acid sequences in FASTA format <br />
`-T`, `--threads`, number of cores to dedicate for parallel processing <br />
`--min`, minimum gene occurrence threshold <br />
`--max`, maximum gene occurrence threshold (should be less than number of genomes, see NOTE below) <br />
`--nucleotide-query`, optional; use this to query nucleotide FASTA files <br />
`--report-all`, optional; use this to return all BLAST hits, regardless of internal identity thresholds <br />
`--perc`, optional; set a custom minimum percent identity threshold. Default = 90% <br />
`--cov`, optional; set a custom minimum query coverage threshold. Default = 75% <br />
`--min-seq-length`, optional; designate minimum allowed query sequence lenght. Use with caution <br />
`--keep-gene-names`, optional; maintains special characters in variable names. Should **not** be used if proceeding to Module 2 

**ReGAIN Curate example Usage:** <br />
`regain curate -d /path/to/genome/files -q /path/to/query/files -T 8 --min 5 --max 475`

**ReGAIN Curate output files:**

`filtered_results.csv`, all BLAST results meeting identity thresholds <br />
`curate_matrix.csv`, filtered data matrix <br />
`curate_metadata.csv`, metadata file for use in ReGAIN statistical modules <br />
**If** `--report-all` **is used:** <br />
`all_results.csv`, all BLAST results, regardless of identity thresholds

NOTE: Discrete Bayesian network anlyses requires all variables to exist in at least two states. For ReGAIN, these two states are 'present' and 'absent'. Ubiquitously occurring genes will break the analysis. 
Best practice is for *N* genomes, `--max` should MINIMALLY be defined as *N* - 1. Keep in mind that removing very low and very high abundance genes can reduce noise in the network.

## **ReGAIN Extract**

ReGAIN Extract is an optional module for use with ReGAIN Curate. This module extracts aligned sequences <br /> 
identified from `regain curate`. Offered as an additional quality control step for gene identification. <br />
Nucleotide sequences are extracted to a multi-FASTA file

`regain extract`

`-c`, `--csv-path`, path to ReGAIN Curate results file, such as `filtered_results.csv` <br />
`-f`, `--fasta-directory`, path to genome FASTA files used in ReGAIN curate <br />
`-T`, `--threads`, number of cores to dedicate for parallel processing <br />
`-o`, `--output-fasta`, multi-FASTA file output (`.fa`, `.fas`, `.fasta`, `.fna`, `.faa`) <br />
`--min-evalue`, optional; for use when `--report-all` flag is used. Sets minimum evalue threshold for sequence extraction <br />
`--min-perc`, optional; same guidelines as `--min-evalue` <br />
`--min-cov`, optional; same guidlines as `--min-evalue` and `--min-perc` <br />
`--translate`, optional; translates extracted nucleotide sequences (see NOTE below)

**ReGAIN Extract example usage:**

`regain extract -c /path/to/results/csv -f /path/to/genome/FASTA/files -T 8 -o sequences.fa`

NOTE: the `--translate` flag should be used with care. In the event an alignment returns an incomplete CDS, <br />
ReGAIN Extract will trim the sequence to the closest value divisible by 3 for codon prediction, which can result <br />
in frameshifts. `--translate` is only suggested for use if returned alignments represent full coding sequences, or <br />
manual validation of gene calls is performed

## **ReGAIN Combine**

ReGAIN Combine is an optional module for use in combination with the ReGAIN Curate and ReGAIN AMR modules. <br />
In the event users want to supplement the `regain AMR` results with a custom set of genes queried through <br />
`regain curate`, `regain combine` will merge both datasets into a single dataset for use in ReGAIN statistical modules

`regain combine`

`--matrix1`, path to ReGAIN AMR data matrix, `filtered_matrix.csv` <br />
`--matrix2`, path to ReGAIN Curate data matrix, `curate_matrix.csv` <br />
`--metadata1`, path to ReGAIN AMR metadata file, `metadata.csv` <br />
`--metadata2`, path to ReGAIN Curate metadata file, `curate_metadata.csv` <br />
`--delete-duplicates`, optional; automatically delete duplicate values from dataset

**ReGAIN Combine example usage:**

`regain combine --matrix1 /path/to/AMR/matrix/csv --matrix2 /path/to/curate/matrix/csv` <br />
`--metadata1 /path/to/AMR/metadata/csv --metadata2 /path/to/curate/metadata/csv` 

**ReGAIN Combine output files:**

`combined_matrix.csv`, combined presence/absence matrix <br />
`combined_metadata.csv`, combined metadata file

NOTE: in order for `regain combine` to function properly, do not modify values in column 1 (`file`) of the data matrix files

_________________________________________________________________________________
## **ReGAIN Accessary Modules**

**Stand Alone Network Visualization**

`regain network`

`-i`, `--input`, input RDS file generated from `bnS`/`bnL` analysis <br />
`-d`, `--data`, input filtered data matrix file <br />
`-M`, `--metadata`, input metadata file <br />
`-s`, `--statistics_results`, input 'Results.csv' file from `bnS`/`bnL` analysis

**Example usage:**

`regain network -i network.rds -d matrix_filtered.csv -M metadata.csv -s Results.csv`

This analysis is an integrated part of the standard `bnS`/`bnL` pipeline, but serves as a redundant measure in the event network visualization needs to be re-performed

**Output:**

`Bayesian_Network.html`, interactive Bayesian network <br />
_________________________________________________________________________________

**Multidimensional Analyses**

Optional Module 3 `regain MVA`

**Currently supported measures of distance:**

`manhattan`, `euclidean`, `canberra`, `clark`, `bray`, `kulczynski`, `jaccard`, `gower`, <br />
`horn`, `mountford`, `raup`, `binomial`, `chao`, `cao`, `mahalanobis`, `altGower`, `morisita`, <br />
`chisq`, `chord`, `hellinger`
                                           
`-i`, `--input`, input file in CSV format <br />
`-m`, `--method`, measure of distance method <br />
`-c`, `--centers`, number of centers (1-10) <br />
`-C`, `--confidence`, confidence interval for ellipses <br />
                                       
**Module 3 example usage:**

`regain MVA -i matrix.csv -m jaccard -c 3 -C 0.75`

**NOTE: the MVA analysis will generate 2 files: a PNG and a PDF of the plot**

_________________________________________________________________________________

## **Formatting External Data**

Bayesian network analysis requires both data matrix and metadata files. MVA analysis requires only a data matrix file <br />
Metadata file <ins>**MUST**<ins> have two column headers. Ideally, 'Gene' and 'GeneClass'. Second column may contain empty rows <br />
Data matrix <ins>**MUST**<ins> have headers for all columns

<img src="https://github.com/ERBringHorvath/regain_cl/assets/97261650/906de456-8368-4872-97c1-df3c9978d535" alt="image">

_________________________________________________________________________________

# Citing ReGAIN

Resistance Gene Association and Inference Network (ReGAIN): A Bioinformatics Pipeline for Assessing Probabilistic 
Co-Occurrence Between Resistance Genes in Bacterial Pathogens. <br />
Bring Horvath, E; Stein, M; Mulvey, MA; Hernandez, EJ; Winter, JM. <br />
*bioRxiv* 2024.02.26.582197; doi: https://doi.org/10.1101/2024.02.26.582197
