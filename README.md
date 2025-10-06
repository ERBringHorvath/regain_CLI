# ReGAIN: <ins>Re</ins>sistance <ins>G</ins>ene <ins>A</ins>ssociation and <ins>I</ins>nference <ins>N</ins>etwork

<img src="https://github.com/user-attachments/assets/149bad0c-fb6b-44c7-a8c6-a131eb979172" alt="image" width="363" height="418"/>

_________________________________________________________________________________
[![bioconda](https://img.shields.io/conda/v/bioconda/regain-cli?color=yellow&logo=anaconda&label=bioconda)](https://anaconda.org/bioconda/regain-cli)
![Python](https://img.shields.io/badge/python-3.10%20-green)
![R](https://img.shields.io/badge/R-%E2%89%A5%204-blue)
[![License](https://img.shields.io/github/license/ERBringHorvath/regain_CLI?color=lavender)](./LICENSE)
![OS](https://img.shields.io/badge/platform-linux--64%20%7C%20osx--64%20%7C%20osx--arm64-red)
[![Issues](https://img.shields.io/github/issues/ERBringHorvath/regain_CLI)](https://github.com/OWNER/REPO/issues)
_________________________________________________________________________________
# ReGAIN Overview
Antimicrobial resistance is driven not only by single genes but by coordinated gene sets that co-occur and move on mobile elements, yet researchers lack a unified, reproducible way to quantify these patterns across genomic populations. The <ins>Re</ins>sistance <ins>G</ins>ene <ins>A</ins>ssociation and <ins>I</ins>nference <ins>N</ins>etwork (ReGAIN) is an open-source platform that utilizes Bayesian network structure learning to infer probabilistic co-occurrence among resistance and virulence genes in clinically important bacterial pathogens. ReGAIN delivers interpretable metrics (conditional probability, relative risk, absolute risk difference, confidence intervals) and post-hoc bidirectional probability scores, enabling detection of synergistic and mutually exclusive gene relationships. By standardizing analyses across species and studies, ReGAIN supports surveillance, hypothesis generation, and stewardship by highlighting gene constellations linked to multidrug resistance and potential co-selection of resistance determinants.

### Explanation of Probability Metrics Used By ReGAIN

1. Conditional Probability
    * Conditional probability is described as the probability of observing Variable A given the presence of Variable B within a given dataset, *P(**A**|**B**)*
    * This value is reported on a scale of 0–1 (conditional probability of 0.5 = 50%)
2. Absolute Risk Difference
    * Absolute risk difference is defined as the conditional probability of observing Variable A given the presence of Variable B <ins>minus</ins> the conditional probability of observing Variable A in the *absence* of Variable B, *P(**A**|**B**) - P(**A**|&not;B)*
    * This value is reported on a scale of -1–1; negative values can occurr if *P(**A**|&not;B)* >> *P(**A**|**B**)* and indicate a strong negative probabilistic relationship
3. Relative Risk
    * Relative risk is a ratio of the conditional probability of observing Variable A given the presence of Variable B to the conditional probability of observing Variable A given the *absence* of Variable B, *P(**A**|**B**)* / *P(**A**|&not;B)*
        * Relative risk < 1, negative relationship
        * Relative risk > 1, positive relationship
        * Relative risk = 1, neutral relationship (variable independence)
4. Baseline risk (`bnS` pipeline only)
    * Baseline risk of outcome: the % occurrence of Variable B in the dataset
    * This value will be identical for all "Variable B's" within a given dataset
    * Example: Gene 1 has a baseline risk of 0.12; this indicates that 12% of the genomic population encodes Gene 1. 

**Post-Hoc Analyses**

A relationship between two variables may be directionally asymmetric — that is, *P(A|B)* does not necessarily equal *P(B|A)* — because the conditional dependencies that determine A given B may differ from those that determine B given A (e.g., knowing that it’s raining (B) changes the probability that the grass is wet (A), but knowing that the grass is wet doesn’t change the probability of rain in exactly the same way — hence *P(A|B)* ≠ *P(B|A)*.) To better inform the direction of the relationship, ReGAIN calculates two post-hoc scores.

5. Bidirectional Probability Score (BDPS)
    * The ratio of the conditional probability *P(A|B)* to the conditional probability *P(B|A)*
        * BDPS > 1, *P(A|B)* > *P(B|A)*
        * BDPS < 1, *P(A|B)* < *P(B|A)*
        * BDPS = 1, equal bidirectional strength
6. Fold Change (FC)
    * The ratio of the relative risk (*P(A|B)* / *P(A|&not;B)* to the relative risk *P(B|A)* / *P(B|&not;A)*) / 2
    * Interpreted similarly to BPDS, with equal bidirectional strength being 0.5
_________________________________________________________________________________

# Table of Contents

* [Installation](#regain-installation)
    * [Bioconda Installation](#bioconda-installation)
    * [Standalone Installation](#standalone-installation)
* [Programs and Example Usage](#programs-and-example-usage)
    * [Basic Network Tutorial](#bayesian-network-structure-learning-basic-tutorial)
* [Module 1: Gene Identification](#resistance-and-virulence-gene-identification)
* [Module 1.1: Dataset Creation](#dataset-creation)
* [Module 2: Bayesian Network Structure Learning](#bayesian-network-structure-learning)
    * [Network Visualization](#standalone-network-visualization)
* [ReGAIN Curate](#regain-curate)
    * [Curate](#regain-curate)
    * [Sequence Extraction](#regain-extract)
    * [Combine Datasets](#regain-combine)
* [Multivariate Analysis](#regain-multivariate-analysis)
* [Formatting External Data](#formatting-external-data)
* [Citations](#citations)
_________________________________________________________________________________

# <ins>**ReGAIN Installation**</ins> 

## Bioconda Installation

1. If not done already, download [miniforge](https://github.com/conda-forge/miniforge/)
2. `conda create -y -n regain python=3.10 r-base=4.4`
    * As of patch 1.6.3, specifying r-base=4.4 during environment creation isn't necessary
3. `source activate regain`
4. `conda install -y bioconda::regain-cli`

One-liner:
`conda create -y -n regain -c bioconda python=3.10 r-base=4.4`

Test installation: <br/>
`regain -h`
`regain --module-health`

If these commands execute without error, the installation is successful.

For `regain --module-health`, you should see:
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

## Standalone Installation
### Install ReGAIN Dependencies

For the most cutting-edge ReGAIN releases, we suggest cloning from source

Create Conda environment and install [NCBI AMRFinderPlus](https://github.com/ncbi/amr/wiki/Install-with-bioconda) and BLAST+
1. `conda create -n regain python=3.10 r-base=4.4`
2. `source activate regain`
3. Install AMRfinderPlus<br/>
    * `conda install -y -c conda-forge -c bioconda ncbi-amrfinderplus`
4. Check AMRFinder installation
    * `amrfinder -h` 
5. Download ARMfinderPlus Database
    * `amrfinder -u`
6. BLAST+ is installed as an AMRFinderPlus dependency

One-liner: `conda create -y -n regain -c conda-forge -c bioconda python=3.10 r-base=4.4 ncbi-amrfinderplus`

## Install ReGAIN

6. Download ReGAIN to preferred directory
* `git clone https://github.com/ERBringHorvath/regain_CLI`

7. Install Python dependencies
* `cd /path/to/regain_CLI`
* `pip install .`

**Note:** R dependencies should be installed automatically when running ReGAIN for the first time

8. Add ReGAIN to your PATH
    * File might be `.bash_profile`, `.bashrc`, `.zshrc` depending on OS and Shell
    * Add this line to the end of the file:

`export PATH="$PATH:/path/to/regain_CLI/src/regain"`

* Replace `/path/to/regain_CLI/src/regain` with the actual path to the directory containing the executable
* Whatever the initial directory, this path should end with `/regain_CLI/src/regain`

9. Save the file and restart your terminal or run `source ~/.bashrc` or `source ~/.zshrc`
10. Verify installation:
* `regain --version`
* `regain --help`
11. Run `regain --module-health` to check status of ReGAIN modules
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
NOTE: ReGAIN utilizes shell scripts to execute some modules. You may need to modify your permissions 
to execute these scripts. If you run `regain --version` and see `permission denied: regain`, Navigate to `regain/src/regain`, then run both `chmod +x regain` and `chmod +x *.sh` and rerun `regain --module-health` <br />
_________________________________________________________________________________

# <ins>**Programs and Example Usage**</ins> 

## Bayesian Network Structure Learning Basic Tutorial

ReGAIN standalone installation comes with a subdirectory called `demo_dataset/` <br/>
For Bioconda installations, `demo_dataset/` files can be downloaded from this repository <br/>
* Input files:
    * Enterococcus_Demo_Dataset.csv
    * Enterococcus_Demo_Dataset_Metadata.csv

Basic `regain bnS` workflow: <br/>
`regain bnS -i Enterococcus_Demo_Dataset.csv \` <br/>
`-M Enterococcus_Demo_Dataset_Metadata.csv \` <br/>
`-o demo_network -n 500 -r 100 -T 4`

This will generate a Bayesian network using 500 bootstraps and 100 data resamples. Using 4 dedicated cores (`-T 4`), this analysis should take less than 30 minutes to run. If more cores are available, setting `-T` to `8` will decrease runtime. 

* Output files:
    * `demo_network.rds`
    * `Query_Results.csv`
    * `post_hoc_analysis.csv`
    * `Bayesian_Network.(html and pdf)`

The resulting network will have both red and black edges connecting nodes. By default, ReGAIN colors negatively associated variable edges red (relative risk < 1). Positively associated edges are colored black (relative risk ≥ 1).  

If specific network visualization parameters are desired, see [Module 2](#standalone-network-visualization) Standalone Network Visualization for available options.

## **Resistance and Virulence Gene Identification** 

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
`-n`, `--number_of_bootstraps`, how many bootstraps to run (suggested minimum of 300-500) <br />
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
`Bayesian_Network.html`, interactive Bayesian network <br/>
`Bayesian_Network.pdf`, static PDF of network using Fruchterman-Reingold force-directed layout algorithm

**NOTE:**<br/>
If using a blacklist, the CSV file should have no headers and only two columns. Column 1 should be your first variable in your variable pair to blacklist, and column 2 should be your second variable in the pair. Explicit bidirectional blacklisting must be passed, so if 'geneA' and 'geneB' are in the blacklist, they should be listed like:

col1, col2 <br/>
geneA, geneB <br/>
geneB, geneA

Additionally, blacklists should be used with caution; ideally, only 'imposible' variables should be blacklisted. An example of this would be a gene with 2 different mutations at the same site (e.g., *gyrA_S83D* and *gyrA_S83E*). 

## Standalone Network Visualization 
For use if `--no-viz` is passed or specific network parameters are wanted

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

## ReGAIN Multivariate Analysis

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
                                       
**Example usage:**

`regain MVA -i matrix.csv -m jaccard --k 0 --pcoa-correction auto`

**NOTE: the MVA analysis will generate 2 files: a PNG and a PDF of the plot**
_________________________________________________________________________________

## Formatting External Data

Bayesian network analysis requires both data matrix and metadata files. MVA analysis requires only a data matrix file <br />
Metadata file <ins>**MUST**<ins> have two column headers. Ideally, 'Gene' and 'GeneClass'. Second column may contain empty rows <br />
Data matrix <ins>**MUST**<ins> have headers for all columns

<img src="https://github.com/ERBringHorvath/regain_cl/assets/97261650/906de456-8368-4872-97c1-df3c9978d535" alt="image">

_________________________________________________________________________________

# Citations

**ReGAIN** <br/>
Bring Horvath E, Stein M, Mulvey MA, Hernandez EJ, Winter JM. <br />
Resistance Gene Association and Inference Network (ReGAIN): A Bioinformatics Pipeline for Assessing Probabilistic Co-Occurrence Between Resistance Genes in Bacterial Pathogens. <br />
*bioRxiv* 2024.02.26.582197; doi: https://doi.org/10.1101/2024.02.26.582197

**AMRFinder** <br/>
Feldgarden M, Brover V, Haft DH, Prasad AB, Slotta DJ, Tolstoy I, Tyson GH, Zhao S, Hsu CH, McDermott PF, Tadesse DA, Morales C, Simmons M, Tillman G, Wasilenko J, Folster JP, Klimke W. <br/>Validating the AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of Isolates. <br/>
Antimicrob Agents Chemother. 2019 Oct 22;63(11):e00483-19. doi: 10.1128/AAC.00483-19

**AMRFinderPlus** <br/>
Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. <br/>
AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <br/>
Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0