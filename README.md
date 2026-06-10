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
    * The ratio of the relative risk *P(A|B)* / *P(A|&not;B)* to the relative risk *P(B|A)* / *P(B|&not;A)*
    * Interpreted similarly to BPDS, with equal bidirectional strength being 1
_________________________________________________________________________________

# Table of Contents

* [Installation](#regain-installation)
  * [Bioconda Installation](#bioconda-installation)
  * [Standalone Installation](#standalone-installation)
* [Programs and Example Usage](#programs-and-example-usage)
  * [Basic Network Tutorial](#bayesian-network-structure-learning-basic-tutorial)
* [Module 1: Gene Identification](#resistance-and-virulence-gene-identification)
  * [Module 1.1: Dataset Creation](#dataset-creation)
  * [Module 1.2: Matrix Summary](#matrix-summary)
  * [Module 1.3: Collapse Features](#collapse-features)
* [Module 2: Bayesian Network Structure Learning](#bayesian-network-structure-learning)
  * [Network Visualization](#standalone-network-visualization)
  * [Module 2.1: Network Analysis](#network-analysis)
* [ReGAIN Curate](#regain-curate)
  * [Curate](#regain-curate)
  * [Sequence Extraction](#regain-extract)
  * [Combine Datasets](#regain-combine)
* [Optional Modules](#regain-accessory-modules)
  * [Multivariate Analysis](#regain-multivariate-analysis)
  * [Genome Similarity](#regain-genome-similarity)
* [Formatting External Data](#formatting-external-data)
* [Citations](#citations)
_________________________________________________________________________________

# <ins>**ReGAIN Installation**</ins> 

### **macOS Note (v1.7.1):** <br/>
Due to current limitations in the Conda R / gRain stack on macOS, the Bioconda recipe for ReGAIN v1.7.1 may fail to solve dependencies
on osx-64. Until this is resolved upstream, we recommend [this workaround](https://github.com/ERBringHorvath/regain_CLI/issues/3) 
using a minimal Python+R Conda environment and installing ReGAIN from source.

## Bioconda Installation

1. If not done already, download [miniforge](https://github.com/conda-forge/miniforge/)
2. `conda create -y -n regain-cli`
3. `source activate regain-conda`
4. `conda install -y -c conda-forge -c bioconda --strict-channel-priority regain-cli=1.8.0`
5. Install AMRFinderPlus database: `amrfinder -u`

Test ReGAIN installation: <br/>
`regain -h`
`regain --module-health`
`amrfinder -h` # Test AMRFinderPlus install

If these commands execute without error, the installation is successful.

For `regain --module-health`:
* You should see:
    ```
    ReGAIN Module Status Report
    - AMR.run: Available
    - matrix.run: Available
    - curate.run: Available
    - extract.run: Available
    - combine.run: Available
    - collapse_features.run: Available
    - matrix_summary.run: Available
    - genome_similarity.run: Available
    - network_analysis.run: Available
    - bnS: Available
    - bnL: Available
    - MVA: Available
    - network: Available
    ```

## Standalone Installation
### Install ReGAIN Dependencies

For the most up-to-date ReGAIN releases, we suggest cloning from source

Create Conda environment and install [NCBI AMRFinderPlus](https://github.com/ncbi/amr/wiki/Install-with-bioconda), BLAST+, FastANI, and Mash
1. `conda create -n regain python=3.10 r-base=4.4`
2. `source activate regain`
3. Install AMRfinderPlus<br/>
    * `conda install -y -c conda-forge -c bioconda ncbi-amrfinderplus fastani mash`
4. Check AMRFinder installation
    * `amrfinder -h` 
5. Download ARMfinderPlus Database
    * `amrfinder -u`
6. BLAST+ is installed as an AMRFinderPlus dependency

One-liner: `conda create -y -n regain -c conda-forge -c bioconda python=3.10 r-base=4.4 ncbi-amrfinderplus fastani mash`

## Install ReGAIN

6. Download ReGAIN to preferred directory (we suggest `$HOME` directory)
* `git clone https://github.com/ERBringHorvath/regain_CLI`

7. Install Python and R dependencies
* `cd /path/to/regain_CLI`
* `pip install .`
* `Rscript install_r_dependencies.R`

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
    - collapse_features.run: Available
    - matrix_summary.run: Available
    - genome_similarity.run: Available
    - network_analysis.run: Available
    - bnS: Available
    - bnL: Available
    - MVA: Available
    - network: Available
    ```
NOTE: ReGAIN utilizes shell scripts to execute some modules. You may need to modify your permissions 
to execute these scripts. If you run `regain --version` and see `permission denied: regain`, Navigate to `regain/src/regain`, then run both `chmod +x regain` and `chmod +x *.sh` and rerun `regain --module-health` <br/>
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
    * `Bayesian_Network.html`/`.pdf`

The resulting network will have both red and black edges connecting nodes. By default, ReGAIN colors negatively associated variable edges red (relative risk < 1). Positively associated edges are colored black (relative risk ≥ 1).  

If specific network visualization parameters are desired, see [Module 2](#standalone-network-visualization) Standalone Network Visualization for available options.

## **Resistance and Virulence Gene Identification** 

## Module 1.0 `regain AMR`

`-f`, `--fasta-directory`, path to directory containing genome FASTA files to analyze <br/>
`--mode`, input FASTA file format (nucleotide, protein, combined (required protein FASTA + GFF file); default = nucleotide) <br/>
`-O`, `--organism`, specify organism for AMRFinderPlus pipeline <br/>
`-o`, `--output-dir`, output directory path <br/>
`-T`, `--threads`, specify number of cores to dedicate

Additional options (AMRFinderPlus-specific parameters): <br/>
`--gff`, input GFF file for combined mode <br/>
`-D`, `--database`, alternate AMRFinderPlus database directory <br/>
`--no-plus`, run basic AMRFinder analysis (organism non-specific; cannot use with `--organism`) <br/>
`--name`, prepend a sample/run identifier to AMRFinderPlus output CSV file <br/>
`--quiet`, suppress AMRFinderPlus status messages <br/>
`--ident-min`, override minimum percent identity (0-1). Use with caution <br/>
`--print-node`, add heierarchy node column <br/>
`--mutation-all`, generate a per-sample point-mutation audit file <br/>
`--nucleotide-output`, write FASTA file of detected nucleotide regions <br/>
`--nucleotide-flank5-output`, write FASTA file of detected nucleotide regions + 5' flanking base pairs <br/>
`--nucleotide-flank5-size`, number of additional 5' flanking base pairs for flank output <br/>
`--protein-output`, write FASTA of detected proteins from input protein FASTA <br/>
`--organism-list`, print the built-in lilst of supported `--organism` options

**Currently supported organisms and how they should be called:**

`Acinetobacter_baumannii` <br/>
`Bordetella_pertussis` <br/>
`Burkholderia_cepacia`<br/>
`Burkholderia_pseudomallei`<br/>
`Burkholderia_mallei` <br/>
`Citrobacter_freundii`<br/>
`Corynebacterium_diphtheriae`<br/>
`Campylobacter`<br/>
`Clostridioides_difficile`<br/>
`Enterobacter_cloacae`<br/>
`Enterobacter_asburiae`<br/>
`Enterococcus_faecalis`<br/>
`Enterococcus_faecium`<br/>
`Escherichia`<br/>
`Haemophilus_influenzae`<br/>
`Klebsiella_oxytoca`<br/>
`Klebsiella_pneumoniae`<br/>
`Neisseria_gonorrhoeae`<br/>
`Neisseria_meningitidis`<br/>
`Pseudomonas_aeruginosa`<br/>
`Salmonella`<br/>
`Serratia_marcescens`<br/>
`Staphylococcus_aureus`<br/>
`Staphylococcus_pseudintermedius`<br/>
`Streptococcus_agalactiae`<br/>
`Streptococcus_pneumoniae`<br/>
`Streptococcus_pyogenes`<br/>
`Vibrio_cholerae`<br/>
`Vibrio_vulnificus`<br/>
`Vibrio_parahaemolyticus`<br/>

**Module 1 example usage:**

Organism specific:
```
regain AMR \
-d path/to/FASTA/files \
-O Pseudomonas_aruginosa \
-T 8 -o path/to/output/directory
```

Organism non-specific:
```
regain AMR \
-d path/to/FASTA/files \
-T 8 -o path/to/output/directory
```

**Output files:**

One result file per submitted genome

_________________________________________________________________________________

# **Dataset Creation** 

### NOTE: variable names <ins>cannot</ins> contain special characters; this transformation is automated during dataset creation <br/>
Examples of special characters: '"/,().[]

## Module 1.1 `regain matrix`
                                       
`-d`, `--directory`, path to AMRfinder results in CSV format <br/>
`--gene-type`, searches for `resistance`, `virulence`, or `all` genes <br/>
`--min`, minimum gene occurrence cutoff <br/>
`--max`, maximum gene occurrence cutoff (should be less than number of genomes, see NOTE below) <br/>
`--report-all`, optional; reports all genes identified, regardless of `--min`/`--max` threshold <br/>
`--keep-gene-names`, optional; maintains special characters in variable names. Should **not** be used if proceeding to Module 2

**Module 1.1 example usage:**

NOTE: Discrete Bayesian network anlyses requires all variables to exist in at least two states. For ReGAIN, these two states are 'present' and 'absent'. Ubiquitously occurring genes will break the analysis. 
Best practice is for *N* genomes, `--max` should MINIMALLY be defined as *N* - 1. Keep in mind that trimming very low and very high abundance genes can reduce noise in the network.
```                                           
regain matrix \
-d path/to/AMRfinder/results \
--gene-type resistance \
--min 5 --max 475
```

**NOTE: all results are saved in the 'ReGAIN_Dataset' folder, which will be generated within the directory <br/> 
defined by** `-d`/`--directory`

**Output files:**

`filtered_matrix.csv`: presence/absence matrix of genes <br/>
`metadata.csv`: file containing genes identified in AMRfinderPlus analysis <br/>
`combined_AMR_results_unfiltered.csv`: concatenated file of all AMRfinder/Plus results; this file contains contig and nucleotide location of all identified genes

If `--report-all` is used: <br/>
`unfiltered_matrix.csv`: presence/absence matrix of all genes identified, regardless of `--min`/`--max` thresholds

_________________________________________________________________________________

# **Matrix Summary**

## Module 1.2 `regain matrix-summary`

`-i`, `--input`, nput presence/absence matrix from `regain matrix` <br/>
`-o`, `--output-file`, output matrix summary file <br/>
`--id-col`, presence/absence matrix genome/sample ID column. Default is first column (header required) <br/>
`--missing-as-zero`, fill missing matrix values with 0 after writing missing-value report

NOTE: Only use `--missing-as-zero` if values are truly absent. Treating empirically missing values as absent <br/>
will affect Bayesian network structure learning results. <br/>
NOTE: value passed to `--output-file` will be used as the basename for all summary files. 

**Example Summary output:** 
| metric |	value |
| ------ | ------ |
| genome_id_column |	file |
| n_genomes	| 267 |
| n_features |	77 |
| duplicated_genome_ids |	0 |
| duplicated_feature_names |	0 |
| missing_values |	0 |
| non_binary_values |	0 |
| min_features_per_genome |	0 |
| max_features_per_genome |	36 |
| mean_features_per_genome |	8.91 |
| median_features_per_genome |	7 |
| sd_features_per_genome |	7.16 |
| n_features_present_in_all_genomes |	0 |
| n_features_absent_in_all_genomes |	0 |
| n_unique_profiles	| 217 |
| n_duplicate_profile_groups |	20 |
| n_genomes_in_duplicate_profiles |	70 |
| percent_genomes_in_duplicate_profiles |	26.22 |
| largest_profile_group_size |	12 |

`genome_id_column`: name of column containing genome/sample IDs <br>
`n_genomes`: number of genomes/samples in dataset <br/>
`n_features`: number of genes/variables in dataset <br/>
`duplicated_genome_ids`: checks if any genome/sample ID appears more than once <br/>
`duplicated_feature_ids`: check if any genes/variables appear more than once <br/>
`missing_values`: checks each row/column for missing data <br/>
`non_binary_values`: checks for non-0/1 data values <br/>
`min_`/`max_`/`mean_`/`median_`/`sd_features_per_genome`: basic stats of genes/variables across genomic population <br/>
`n_features_present`/`absent_in_all_genomes`: checks for genes/variables missing from or occuring in all genomes <br/>
`n_unique_profiles`: each genome is assigned a gene co-occurrence profile; this reports how many of those profiles are unique <br/>
`n_duplicate_profile_groups`: how many shared gene co-occurrence profiles are in the dataset <br/>
`n_genomes_in_duplicate_profiles`: how many genomes comprise duplicate profile groups <br/>
`percent_genomes_in_duplicate_profiles`: percent report of `n_genomes_in_duplicate_profiles` <br/>
`largest_profile_group_size`: how many genomes encode the most common gene co-occurrence profile <br/>

**Example feature_counts.csv output**
|feature	|count	|percent_present|
|---|---|---|
|blaTEM_1	|124	|46.44|
|sul2	|109	|40.82|
|dfrA1	|103	|38.58|
|aph3pp_Ib	|98	|36.7|
|aph6_Id	|98	|36.7|
|qacE	|96	|35.96|
|uhpT_E350Q	|90	|33.71|
|...|...|...|

**Example profile_groups.csv output**
|profile_id|profile_hash|features_present|n_genomes|genomes|n_features_present|status|
|---|---|---|---|---|---|---|
|profile_00018	|a584464eec587b49|	aph6_Id; aph3pp_Ib; sul2; blaTEM_1; gyrA_S83L; uhpT_E350Q; mphA; sul1; qacEdelta1; aadA5; dfrA17; tetA; aac3_IId; ptsI_V25I; parE_I529L; dfrA1; qacE; aac3_II|	2|	ECOLI_17_124_S187; ECOLI_13_126_S24	|18|	duplicate_profile|
|profile_00019	|bdac998e15c72026	|terD; terZ; terW	|2|	ECOLI_15_128_S109; ECOLI_14_134_S64	|3	|duplicate_profile |
|profile_00020	|c264b45c4256fbb5	|sul2; blaTEM_1; marR_S3N; blaEC_5; uhpT_E350Q	|2	|ECOLI_13_124_S22; ECOLI_15_141_S122	|5	|duplicate_profile|
|profile_00021	|e0c9696a80d71a39	|dfrA5; uhpT_E350Q; ptsI_V25I; parE_I529L; ble	|1	|ECOLI_16_117_S140	|5|	unique_profile|
|profile_00022	|26f1e1728c760faf	|tetA; dfrA1; aadA1; tetB; floR; aph3p_Ia; silE; blaCTX_M_55; qnrB19; aac3_IIe; dfrA14; qnrB; blaCTX_M; aac3_II	|1|	ECOLI_16_139_S162	|14|	unique_profile|

`profile_members.csv` is also generated. This is a simplifed version of `profile_groups.csv`
_________________________________________________________________________________

# **Collapse Features**

## Module 1.3 `regain collapse-features`

Regain Collapse Features is an optional module that allows binning of specific determinants into broader categories. Only features assigned a 'bin' will be collapse; other features will be left unchanged unless `--drop-unmapped` is passed. 

`-i`, `--input`, input presence/absence matrix generated from `regain matrix` <br/>
`-M`, `--metadata`, metadata file containing `Gene` and `bin` columns <br/>
`-o`, `--output-file`, output binned presence/absence matrix <br/>
Optional parameters: <br/>
`--id-col`, presence/absence matrix genome/sample ID column (default = first column) <br/>
`--gene-col`, column in metadata containing original feature name (default = Gene) <br/>
`--bin-col`, column in metaadata containing collapsed (binned) feature names (default = bin) <br/>
`--drop-unmapped`, drop matrix features not assigned a bin in metadata file <br/>
`--missing-as-zero`, fill missing matrix values with 0 after writing missing-value report <br/>
**NOTE:** `--missing-as-zero` is <ins>**not**</ins> appropriate for empirically missing data. Only use this parameter if rows were left empty to denote true absence

**Module 1.3 example usage:**

Example modified `regain-matrix` metadata file:
|Gene|bin|GeneClass|...|
|---|---|---|---|
|aac6_Ib_cr|aac6_Ib_cr_group|aminoglycoside/quinolone|
|aac6_Ib_cr5|aac6_Ib_cr_group|aminoglycoside/quinolone|
|aac6_Ib_cr6|aac6_Ib_cr_group|aminoglycoside/quinolone|
|tetA|tet_group|tetracycline|
|tetB|tet_group|tetracycline|
|...||

```
regain collapse-features \
-i filtered_matrix.csv \
-M metadata.csv \
-o collapsed_features_matrix.csv \
--bin-col bin \
--gene-col Gene
```
# **Bayesian Network Structure Learning** 

## Module 2.0 `regain bnL` or `regain bnS`
                                            
`-i`, `--input`, input file in CSV format <br/>
`-M`, `--metadata`, file containing gene names and descriptions <br/>
`-o`, `--output-boot`, output bootstrap file <br/>
`-T`, `--threads`, number of cores to dedicate for parallel processing <br/>
`-n`, `--bootstraps`, how many bootstraps to run (suggested minimum of 300-500) <br/>
`-r`, `--resamples`, how many data resamples you want to use (suggested minimum of 100) <br/>
`--blacklist`, optional blacklist CSV (no header); 2 columns for variable 1 and variables 2 <br/>
`--iss`, imaginary sample size for BDe score (default = 10) <br/>
`--no-viz`, Skip HTML/PDF visualization (for use if specific options are wanted, see [`regain network`](#regain-accessory-modules))

**bnL only:**

`--cp-samples`, Monte Carlo samples for cpquery (default: 10000)

**Module 2 example usage:**

**NOTE: We suggest using between a <ins>minimum</ins> of 300 to 500 bootstraps and >=100 resamples when possible**

`bnS`, Bayesian network structure learning analysis for less than 100 genes <br/>
`bnL`, Bayesian network structure learning analysis for 100 genes or greater

For less than 100 genes:
```
regain bnS \
-i matrix_filtered.csv \
-M metadata.csv \
-o bootstrapped_network \
-T 8 -n 500 -r 500 \
--blacklist list.csv
```                                     
For 100 or more genes:
```
regain bnL \
-i matrix_filtered.csv \
-M metadata.csv \
-o bootstrapped_network \
-T 8 -n 500 -r 500 \
--blacklist list.csv
```

**Output files:**

`<output>.rds`, bootstrapped network <br/>
`Query_Results.csv`, results file of all conditional probability and relative risk values <br/>
`post_hoc_analysis.csv`, results file of all bidirectional probability and fold change scores <br/>
`Bayesian_Network.html`, interactive Bayesian network <br/>
`Bayesian_Network.pdf`, static PDF of network using Fruchterman-Reingold force-directed layout algorithm

**NOTE:**<br/>
If using a blacklist, the CSV file should have no headers and only two columns. Column 1 should be your first variable in your variable pair to blacklist, and column 2 should be your second variable in the pair. Explicit bidirectional blacklisting must be passed, so if 'geneA' and 'geneB' are in the blacklist, they should be listed like:

|col1|col2|
|---|---|
|geneA|geneB|
|geneB|geneA|

Additionally, blacklists should be used with caution; ideally, only 'imposible' variables should be blacklisted. An example of this would be a gene with 2 different mutations at the same site (e.g., *gyrA_S83D* and *gyrA_S83E*). 

**NOTE:** It is possible for datasets of less than 100 genes to fail during `bnS` analysis due to insufficient resources needed for full graph inference. If this occurs, re-run analysis using `regain bnL`. 

## Standalone Network Visualization 
For use if `--no-viz` is passed or specific network parameters are wanted

`regain network <options>`

`-i`, `--input`, input RDS file generated from `bnS`/`bnL` analysis <br/>
`-d`, `--data`, input filtered data matrix file <br/>
`-M`, `--metadata`, input metadata file <br/>
`-s`, `--statistics-results`, input 'Results.csv' file from `bnS`/`bnL` analysis <br/>
`--threshold`, averaged network threshold (default: 0.5) <br/>
`--seed`, layout seed for Fruchterman-Reingold force-directed layout algorithm (pdf output) <br/>
`--html-out`, HTML output file name (default: Bayesian_Network.html) <br/>
`--pdf-out`, PDF output file name (default: Bayesian_Network.pdf) <br/>
`-b`, `--blacklist`, optional blacklist CSV (no header): from,to <br/>
`--width-metric`, edge-width metric selection (auto, abs_mean, abs_ci, cp_ci, cp_mean) (default: auto (abs_ci)) <br/>
`--rr-metric`, relative risk edge color threshold (default: 1.0, <1: red, >=1: black)

**Example usage:**
```
regain network \
-i network.rds \
-d matrix_filtered.csv \
-M metadata.csv \
-s Query_Results.csv
```

This workflow is an integrated part of the standard `bnS`/`bnL` pipeline, but serves as a redundant measure in the event network visualization needs to be re-performed or specific parameters beyond `bnS`/`bnL` defaults are wanted.
_________________________________________________________________________________

# Network Analysis

## Module 2.1 `regain network-analysis`
Summarizes gene pair-level differences between two `bnS`/`bnL`-generated `Query_Results.csv`

`--network1`, baseline probability results table <br/>
`--network2`, comparison results table <br/>
`-o`, `--output-dir`, Output directory for comparison results

Optional parameters for custom datasets: <br/>
`--gene-a-col`, column name for variable 1. Default: Gene_1 (standard ReGAIN output) <br/>
`--gene-b-col`, column name for variable 2. Default: Gene_2 (standard ReGAIN output) <br/>
`--cpr-mean-col`, conditional probability mean column. Default: Conditional_Probability_Mean <br/>
`--rr-mean-col`, relative risk mean column. Default: Relative_Risk_Mean <br/>
`--ard-mean-col`, absolute risk difference column. Default: Absolute_Risk_Mean <br/>
`--status-metric`, metric used to classify retained/weakened/strengthened. Default: ard <br/>
`--delta-threshold`, minimum mean delta magnitude used to define weakened/strengthened. Default: 0.1 <br/>

NOTE: `--delta-threshold` default value is just a suggestion. A change of 0.1 in absolute risk difference indicates a
10% variance from one network to the next. If relative risk is used as `--status-metric`, we suggest increasing the
`--delta-threshold` to accomodate the change in expected scale (ARD scale is 0-1, expressed as a decimal, RR scale is 0-&infin;)

**Example usage:** 
```
regain network-analysis \
--network1 Query_Results_n1.csv \
--network2 Query_Results_n2.csv \
-o /path/to/output/dir
```

**Example summary.csv:**
|metric	|value|
|---|---|
|n_pairs_network1|	2556|
|n_pairs_network2	|2926|
|n_shared_pairs	|2415|
|n_network1_only_pairs|	141|
|n_network2_only_pairs|	511|
|n_retained_pairs	|2223|
|n_weakened_pairs	|159|
|n_strengthened_pairs	|33|
|mean_delta_cpr_ab	|0.00039|
|mean_delta_cpr_ba	|-0.00229|
|mean_delta_rr_ab	|-8.96654|
|mean_delta_rr_ba	|-14.69088|
|mean_delta_ard_ab	|-0.01427|
|mean_delta_ard_ba	|-0.0129|

`n_pairs_network1`: total number of gene pairs seen in `--network1` <br/>
`n_pairs_network2`: total number of gene pairs seen in `--network2` <br/>
`n_shared_paris`: total number of gene pairs shared between both input results <br/>
`n_network1_only_pairs`: total number of gene pairs present only in `--network1` <br/>
`n_network2_only_pairs`: total number of gene pairs present only in `--network2` <br/>
`n_retained_pairs`: number of gene pairs showing variance < `--delta-threshold` <br/>
`n_weakened_pairs`: number of gene pairs showing decreased probability of co-occurrence > `--delta-threshold` <br/>
`n_strengthened_pairs`: number of gene pairs showing increased probability of co-occurece > `--delta-threshold` <br/>
`mean_delta_*`: average change in conditional probability (cpr), relative risk (rr), and absolute risk difference (ard) <br/>
`*_ab`: direction of relationship, *P(A|B)* <br/>
`*_ba`: direction of relationship, *P(B|A)* <br/>

**NOTE:** `n_weakened_pairs`/`n_strengthened_pairs` reporting is relative to `--status-metric` and `--delta-threshold` <br/>
Example using the table above: <br/>
If `--status-metric` is `ard` and `--delta-threshold` is `0.1`, we can interpret `n_weakened_pairs` as 159 / 2415 (6.6%) of gene pairs showed
a >10% decreased probability of co-occurrence based on ARD, while 33 / 2415 (1.4%) of gene pairs showed a >10% increased
probability of co-occurrence based on ARD. 
_________________________________________________________________________________
# **ReGAIN Curate**

ReGAIN Curate is designed to allow users to generate a dataset for Bayesian network structure learning using <br/>
a custom set of gene queries, independent of ReGAIN Module 1

`regain curate <options>`

`-d`, `--directory`, path to genome FASTA files <br/>
`-q`, `--query`, path to query files containing amino acid sequences in FASTA format <br/>
`-T`, `--threads`, number of cores to dedicate for parallel processing <br/>
`--min`, minimum gene occurrence threshold <br/>
`--max`, maximum gene occurrence threshold (should be less than number of genomes, see NOTE below) <br/>
`--nucleotide-query`, optional; use this to query nucleotide FASTA files <br/>
`--report-all`, optional; use this to return all BLAST hits, regardless of internal identity thresholds <br/>
`--evalue`, optional; set a custom maximum e-value threshold. Default = 1e-5 <br/>
`--perc`, optional; set a custom minimum percent identity threshold. Default = 90% <br/>
`--cov`, optional; set a custom minimum query coverage threshold. Default = 75% <br/>
`--min-seq-length`, optional; designate minimum allowed query sequence lenght. Use with caution <br/>
`--keep-gene-names`, optional; maintains special characters in variable names. Should **not** be used if proceeding to Module 2 

**ReGAIN Curate example Usage:** <br/>
```
regain curate \
-d /path/to/genome/files \
-q /path/to/query/files \
-T 8 --min 5 --max 475
```

**ReGAIN Curate output files:**

`filtered_results.csv`, all BLAST results meeting identity thresholds <br/>
`curate_matrix.csv`, filtered data matrix <br/>
`curate_metadata.csv`, metadata file for use in ReGAIN statistical modules <br/>
**If** `--report-all` **is used:** <br/>
`all_results.csv`, all BLAST results, regardless of identity thresholds

NOTE: Discrete Bayesian network anlyses requires all variables to exist in at least two states. For ReGAIN, these two states are 'present' and 'absent'. Ubiquitously occurring genes will break the analysis. 
Best practice is for *N* genomes, `--max` should MINIMALLY be defined as *N* - 1. Keep in mind that removing very low and very high abundance genes can reduce noise in the network.

## **ReGAIN Extract**

ReGAIN Extract is an optional module for use with ReGAIN Curate. This module extracts aligned sequences <br/> 
identified from `regain curate`. Offered as an additional quality control step for gene identification. <br/>
Nucleotide sequences are extracted to a multi-FASTA file

`regain extract <options>`

`-c`, `--csv-path`, path to ReGAIN Curate results file, such as `filtered_results.csv` <br/>
`-f`, `--fasta-directory`, path to genome FASTA files used in ReGAIN curate <br/>
`-T`, `--threads`, number of cores to dedicate for parallel processing <br/>
`-o`, `--output-fasta`, multi-FASTA file output (`.fa`, `.fas`, `.fasta`, `.fna`, `.faa`) <br/>
`--evalue`, optional; for use when `--report-all` flag is used. Sets maximum evalue threshold for sequence extraction <br/>
`--perc`, optional; same guidelines as `--evalue` <br/>
`--cov`, optional; same guidlines as `--evalue` and `--min-perc` <br/>
`--translate`, optional; translates extracted nucleotide sequences (see NOTE below)

**ReGAIN Extract example usage:**
```
regain extract \
-c /path/to/results/csv \
-f /path/to/genome/FASTA/files \
-T 8 -o sequences.fa
```

NOTE: the `--translate` flag should be used with care. In the event an alignment returns an incomplete CDS, <br/>
ReGAIN Extract will trim the sequence to the closest value divisible by 3 for codon prediction, which can result <br/>
in frameshifts. `--translate` is only suggested for use if returned alignments represent full coding sequences, or <br/>
manual validation of gene calls is performed

## **ReGAIN Combine**

ReGAIN Combine is an optional module for use in combination with the ReGAIN Curate and ReGAIN AMR modules. <br/>
In the event users want to supplement the `regain AMR` results with a custom set of genes queried through <br/>
`regain curate`, `regain combine` will merge both datasets into a single dataset for use in ReGAIN statistical modules

`regain combine <options>`

`--matrix1`, path to ReGAIN AMR data matrix, `filtered_matrix.csv` <br/>
`--matrix2`, path to ReGAIN Curate data matrix, `curate_matrix.csv` <br/>
`--metadata1`, path to ReGAIN AMR metadata file, `metadata.csv` <br/>
`--metadata2`, path to ReGAIN Curate metadata file, `curate_metadata.csv` <br/>
`--delete-duplicates`, optional; automatically delete duplicate values from dataset

**ReGAIN Combine example usage:**
```
regain combine \
--matrix1 /path/to/AMR/matrix/csv \
--matrix2 /path/to/curate/matrix/csv \
--metadata1 /path/to/AMR/metadata/csv \
--metadata2 /path/to/curate/metadata/csv 
```

**ReGAIN Combine output files:**

`combined_matrix.csv`, combined presence/absence matrix <br/>
`combined_metadata.csv`, combined metadata file

NOTE: in order for `regain combine` to function properly, do not modify values in column 1 (`file`) of the data matrix files

_________________________________________________________________________________
# **ReGAIN Accessory Modules**

## ReGAIN Multivariate Analysis

**Currently supported measures:**

`manhattan`, `euclidean`, `canberra`, `clark`, `bray`, `kulczynski`, <br/>
`jaccard`, `gower`, `altGower`, `morisita`, `horn`, `mountford`, `raup`, <br/>
`binomial`, `chao`, `cao`, `chord`, `hellinger`, `aitchison`, `mahalanobis`

`regain MVA <options>`
                                           
`-i`, `--input`, input file in CSV format <br/>
`-m`, `--method`, options:  <br/>
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
```
regain MVA \
-i matrix.csv \
-m jaccard \
--k 0 --pcoa-correction auto
```

**NOTE: the MVA analysis will generate 2 files: a PNG and a PDF of the plot**
_________________________________________________________________________________

## ReGAIN Genome Similarity

ReGAIN `genome-similarity` is an optional module to assess population similarity using either average nucleotide identity (ANI) or Mash distance. This module is not a part of the default ReGAIN pipeline and is primarily useful for custom datasets where there is concern that highly-similar (e.g., clonal) lineages may impact Bayesian network structure learning results. 

**NOTE:** This module does not automatically trim genomes from an input population. Instead, it provides similarity results for user review to direct population adjustments, if needed. 

`regain genome-similarity <options>`

`-f`, `--fasta-dir`, directory of genomes in FASTA format for comparison <br/>
`-o`, `--output-dir`, output directory for genome reports <br/>
`--method`, genome similarity method {fastani, mash}. Default = fastani <br/>
`--threshold`, one or more space-separated similarity thresholds. For FastANI, thresholds are ANI percent values (default = 99.9 99.5 99.0). For Mash, thresholds are Mash distance values (default = 0.001 0.005 0.01) <br/>
`--genome-list`, optional file containing one genome filename or relative path per row. Entries are resolved relative to `--fasta-dir` <br/>
`-e`, `--extension`, FASTA extensions (space-separated) to include. Default: .fna, .fa, .fas, .ffn, .fasta <br/>
`--recursive`, search `--fasta-dir` recursively for FASTA files (one or more levels deep) <br/>
`-T`, `--threads`, number of cpus to pass to selected method (default = 4) <br/>
`--sketch-size`, Mash sketch size. Used only with `--method mash`. Default: 10000 <br/>
`--kmer-size`, Mash k-mer size. Used only with `--method mash`. Default: 21

**Example usage:**
```
regain genome-similarity \
-f /path/to/genome/FASTA/files \
-o /path/to/results/dir \
--method fastani \
-T 4 
```
_________________________________________________________________________________

## Formatting External Data

Bayesian network analysis requires both data matrix and metadata files. MVA analysis requires only a data matrix file <br/>
Metadata file <ins>**MUST**<ins> have two column headers. Ideally, 'Gene' and 'GeneClass'. Second column may contain empty rows <br/>
Data matrix <ins>**MUST**<ins> have headers for all columns

<img src="https://github.com/ERBringHorvath/regain_cl/assets/97261650/906de456-8368-4872-97c1-df3c9978d535" alt="image">

_________________________________________________________________________________

# Citations

**ReGAIN** <br/>
Bring Horvath E, Stein M, Mulvey MA, Hernandez EJ, Winter JM. <br/>
ReGAIN: A Bioinformatics Pipeline for Assessing Probabilistic Co-Occurrence Between Resistance Genes in Bacterial Pathogens. <br/>
*bioRxiv* 2024.02.26.582197; doi: https://doi.org/10.1101/2024.02.26.582197

**AMRFinder** <br/>
Feldgarden M, Brover V, Haft DH, Prasad AB, Slotta DJ, Tolstoy I, Tyson GH, Zhao S, Hsu CH, McDermott PF, Tadesse DA, Morales C, Simmons M, Tillman G, Wasilenko J, Folster JP, Klimke W. <br/>Validating the AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of Isolates. <br/>
Antimicrob Agents Chemother. 2019 Oct 22;63(11):e00483-19. doi: 10.1128/AAC.00483-19

**AMRFinderPlus** <br/>
Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. <br/>
AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <br/>
Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0

**FastANI** <br/>
Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. <br/>
High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. <br/>
Nat Commun. 2018;9(1):5114. doi:10.1038/s41467-018-07641-9

**Mash** <br/>
Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. <br/>
Mash: fast genome and metagenome distance estimation using MinHash. <br/>
Genome Biol. 2016 Jun 20;17(1):132. doi: 10.1186/s13059-016-0997-x