# ReGAIN: <ins>Re</ins>sistance <ins>G</ins>ene <ins>A</ins>ssociation and <ins>I</ins>nference <ins>N</ins>etwork

<img src="https://github.com/user-attachments/assets/149bad0c-fb6b-44c7-a8c6-a131eb979172" alt="image" width="363" height="418"/>

_________________________________________________________________________________
![Python](https://img.shields.io/badge/python-3.10%20-green)
![R](https://img.shields.io/badge/R-%E2%89%A5%204-blue)
[![License](https://img.shields.io/github/license/ERBringHorvath/regain_CLI)](./LICENSE)
![OS](https://img.shields.io/badge/platform-linux--64%20%7C%20osx--64%20%7C%20osx--arm64-red)
[![Issues](https://img.shields.io/github/issues/ERBringHorvath/regain_CLI)](https://github.com/OWNER/REPO/issues)
_________________________________________________________________________________

# <ins>**ReGAIN Installation and User guide**</ins> 


**Prerequisites**

Ensure that you have the following prerequisites installed on your system:

Python (version 3.10 or higher)

R (version 4 or higher)

NCBI AMRfinderPlus version 4.0.3 <br />
NCBI BLAST+

[Install R](https://www.r-project.org/)

_________________________________________________________________________________

## Install ReGAIN Dependencies

**We suggest that ReGAIN and all prerequisites are installed within a Conda environment**

Download [miniforge](https://github.com/conda-forge/miniforge/)

1. Create Conda environment and install [NCBI AMRfinderPlus](https://github.com/ncbi/amr/wiki/Install-with-bioconda)
* `conda create -n regain python=3.10`
* `source activate regain`

2. Install AMRfinderPlus<br/>
* `conda install -y -c conda-forge -c bioconda ncbi-amrfinderplus`

3. Check installation
* `amrfinder -h`

4. Download ARMfinderPlus Database
* `amrfinder -u`

5. Install NCBI BLAST+
* `conda install -y -c bioconda::blast`

## Install ReGAIN

6. Download ReGAIN to preferred directory
* `git clone https://github.com/ERBringHorvath/regain_CLI`

7. Install Python dependencies
* `cd /path/to/regain_CLI`
* `pip install -e .`

8. Add ReGAIN to your PATH
* Add this line to the end of `.bash_profile`/`.bashrc` (Linux) or `.zshrc` (macOS):

`export PATH="$PATH:/path/to/regain_CLI/src/regain"`

* Replace `/path/to/regain_CLI/bin` with the actual path to the directory containing the executable. <br />
Whatever the initial directory, this path should end with `/regain_CLI/src/regain`

9. Save the file and restart your terminal or run `source ~/.bash_profile` or `source ~/.zshrc`

10. Verify installation:
* `regain --version`
* use `-h`, `--help`, to bring up the help menu
* `regain --help`
* run `regain --module-health` to check status of ReGAIN modules
    * You should see:

    ```
    ReGAIN Module Status Report
    - AMR.run: Available
    - matrix.run: Available
    - curate.run: Available
    - extract.run: Available
    - combine.run: Available
    - bnS: Available
    - bnL: Available
    - MVA: Available
    - network: Available
    ```

NOTE: ReGAIN utilizes shell scripts to execute some modules. You may need to modify your permissions <br />
to execute these scripts. If you run `regain --version` and see `permission denied: regain`, Navigate to <br />
`regain/src/regain`, then run both `chmod +x regain` and `chmod +x *.sh` and rerun `regain --version` <br />

_________________________________________________________________________________

# <ins>**Programs and Example Usage**</ins> 

# **Resistance and Virulence Gene Identification** 

## Module 1.0 `regain AMR`

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

# **Dataset Creation** 

### NOTE: variable names <ins>cannot</ins> contain special characters; this transformation is automated during dataset creation <br />
Examples of special characters: '"/,().[]

## Module 1.1 `regain matrix`
                                       
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
`combined_AMR_results_unfiltered.csv`: concatenated file of all AMRfinder/Plus results; this file contains contig and nucleotide location of all identified genes

If `--report-all` is used: <br />
`unfiltered_matrix.csv`: presence/absence matrix of all genes identified, regardless of `--min`/`--max` thresholds

_________________________________________________________________________________

# **Bayesian Network Structure Learning** 

## Module 2 `regain bnL` or `regain bnS`
                                            
`-i`, `--input`, input file in CSV format <br />
`-M`, `--metadata`, file containing gene names and descriptions <br />
`-o`, `--output_boot`, output bootstrap file <br />
`-T`, `--threads`, number of cores to dedicate for parallel processing <br />
`-n`, `--number_of_boostraps`, how many bootstraps to run (suggested minimum of 300-500) <br />
`-r`, `--number-of-resamples`, how many data resamples you want to use (suggested minimum of 100) <br />
`--blacklist`, optional blacklist CSV (no header); 2 columns for variable 1 and variables 2 <br/>
`--iss`, imaginary sample size for BDe score (default = 10) <br/>
`--no-viz`, Skip HTML/PDF visualization (for use if specific options are wanted, see [`regain network`](#regain-accessory-modules))

**bnL only:**

`--cp-samples`, Monte Carlo samples for cpquery (default: 10000)

**Module 2 example usage:**

**NOTE: We suggest using between a <ins>minimum</ins> of 300 to 500 bootstraps and >100 resamples when possible**

`bnS`, Bayesian network structure learning analysis for less than 100 genes <br />
`bnL`, Bayesian network structure learning analysis for 100 genes or greater

For less than 100 genes:

`regain bnS -i matrix_filtered.csv -M metadata.csv -o bootstrapped_network -T 8 -n 500 -r 500 --blacklist list.csv`
                                            
For 100 or more genes:

`regain bnL -i matrix_filtered.csv -M metadata.csv -o bootstrapped_network -T 8 -n 500 -r 500 --blacklist list.csv`

**Output files:**

`<output>.rds`, bootstrapped network <br/>
`Query_Results.csv`, results file of all conditional probability and relative risk values <br />
`post_hoc_analysis.csv`, results file of all bidirectional probability and fold change scores <br />
`Bayesian_Network.html`, interactive Bayesian network
`Bayesian_Network.pdf`, static PDF of network using Fruchterman-Reingold force-directed layout algorithm

**NOTE:**<br/>
If using a blacklist, the CSV file should have no headers and only two columns. Column 1 should be your first variable in your variable pair to blacklist, and column 2 should be your second variable in the pair. Explicit bidirectional blacklisting must be passed, so if 'geneA' and 'geneB' are in the blacklist, they should be listed like:

col1, col2 <br/>
gene1, gene2 <br/>
gene2, gene1

Additionally, blacklists should be used with caution; ideally, only 'imposible' variables should be blacklisted. An example of this would be a gene with 2 different mutations at the same site (e.g., *gyrA_S83D* and *gyrA_S83E*). 

## Stand Alone Network Visualization (for use if `--no-viz` is passed)

`regain network`

`-i`, `--input`, input RDS file generated from `bnS`/`bnL` analysis <br />
`-d`, `--data`, input filtered data matrix file <br />
`-M`, `--metadata`, input metadata file <br />
`-s`, `--statistics_results`, input 'Results.csv' file from `bnS`/`bnL` analysis <br/>
`--threshold`, averaged network threshold (default: 0.5) <br/>
`--seed`, layout seed for Fruchterman-Reingold force-directed layout algorithm (pdf output) <br/>
`--html-out`, HTML output file name (default: Bayesian_Network.html) <br/>
`--pdf-out`, PDF output file name (default: Bayesian_Network.pdf) <br/>
`-b`, `--blacklist`, optional blacklist CSV (no header): from,to <br/>
`--width-metric`, edge-width metric selection (auto, abs_mean, abs_ci, cp_ci, cp_mean) (default: auto (abs_ci)) <br/>
`--rr-metric`, relative risk edge color threshold (default: 1.0, <1: red, >=1: black)

**Example usage:**

`regain network -i network.rds -d matrix_filtered.csv -M metadata.csv -s Results.csv`

This workflow is an integrated part of the standard `bnS`/`bnL` pipeline, but serves as a redundant measure in the event network visualization needs to be re-performed or specific parameters beyond `bnS`/`bnL` defaults are wanted.
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
# **ReGAIN Accessory Modules**

**Stand Alone Network Visualization**

`regain network`

`-i`, `--input`, input RDS file generated from `bnS`/`bnL` analysis <br />
`-d`, `--data`, input filtered data matrix file <br />
`-M`, `--metadata`, input metadata file <br />
`-s`, `--statistics-results`, input 'Results.csv' file from `bnS`/`bnL` analysis <br/>
`--threshold`, average network threshold (Default = 0.5) <br/>
`--seed`, set seed for Fruchterman-Reingold force-directed layout algorithm (PDF only, Default = 42) <br/>
`--html-out`, HTML output file name <br/>
`--pdf-out`, PDF output file name <br/>
`-b`, `--blacklist`, optional blacklist CSV (no header): from,to <br/>
`--width-metric`, edge width metric selection (auto, abs_mean, abs_ci, cp_ci, cp_mean) (Default = abs_ci) <br/>
`--rr-threshold`, relative risk edge color threshold (Default = 1, <1: red, >=1: black)

**Example usage:**

`regain network -i network.rds -d matrix_filtered.csv -M metadata.csv -s Results.csv`

This analysis is an integrated part of the standard `bnS`/`bnL` pipeline, but serves as a redundant measure in the event network visualization needs to be re-performed or specific options are wanted

**Default Output:**

`Bayesian_Network.html`, interactive Bayesian network <br />
`Bayesian_Network.pdf`, static network generated using Fruchterman-Reingold force-directed layout algorithm
_________________________________________________________________________________

**Multidimensional Analyses**

## Optional Module 3 `regain MVA`

**Currently supported measures:**

`manhattan`, `euclidean`, `canberra`, `clark`, `bray`, `kulczynski`, <br/>
`jaccard`, `gower`, `altGower`, `morisita`, `horn`, `mountford`, `raup`, <br/>
`binomial`, `chao`, `cao`, `chord`, `hellinger`, `aitchison`, `mahalanobis`

`regain MVA`
                                           
`-i`, `--input`, input file in CSV format <br />
`-m`, `--method`, options:  <br />
`--k`, k for k-means; 0 = auto (2..10) <br/>
`-C`, `--confidence`, Ellipse confidence (Default = 0.95) <br/>
`--seed`, set seed for reproducibility (default = 42) <br/>
`--label`, apply labels to data points [none, auto, all] (default = auto) <br/>
`--point-size`, size of data points (default = 3.5) <br/>
`--alpha`, opacity of data points (default = 0.75) <br/>
`--no-ellipses`, do not layer confidence ellipses <br/>
`--psuedocount`, for Aitchison method (default = 1e-6) <br/>
`--save-dist`, save distances to CSV <br/>
`--dist-out`, manually name output distance file <br/>
`--png-out`, manually name output PNG <br/>
`--pdf-out`, manually name output PDF <br/>
`--coords-out`, manually name output coordinate file (default = MVA_coordinates.csv) <br/>
`--pcoa-correction`, apply PCoA correction [auto, none, lingoes, cailliez] (default = auto)
                                       
**Module 3 example usage:**

`regain MVA -i matrix.csv -m jaccard --k 0 --pcoa-correction auto`

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
