This demo dataset was created using 502 *Enterococcus faecium* genomes, resulting in 39 resistance genes <br/>
This dataset is meant to be used starting with Module 2, Bayesian network structure learning.

Example Module 2 command:

`regain bnS -i Enterococcus_Demo_Dataset.csv -M Enterococcus_Demo_Dataset_Metadata.csv -o TEST -T 4 -n 500 -r 100` <br />
or <br />
`regain bnL -i Enterococcus_Demo_Dataset.csv -M Enterococcus_Demo_Dataset_Metadata.csv -o TEST -T 4 -n 500 -r 100`

Example Multivariate Analysis command:

`regain MVA -i Enterococcus_Demo_Dataset.csv -m jaccard -c 3 -C 0.95`

This specifies three clusters with 95% confidence using the Jaccard measure of similarity/dissimilarity
