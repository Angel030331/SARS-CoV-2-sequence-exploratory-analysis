# SARS-CoV-2-sequence-exploratory-analysis
This repository serves as a documentation on SARS-CoV-2 analysis trials. Different bioinformatic tools and experimental set up may be involved. This documentation will contribute to a transdisciplinary research project regarding SARS-CoV-2 evolution.

## Environmnet set up
### 1. Download sequences from NCBI 
ref sequence: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta

sample sequence: https://www.ncbi.nlm.nih.gov/nuccore/OP599898.1

### 2. Download MAFFT via conda
Assume that the conda environment is set up in the terminal, type the lines below to download the required tool, MAFFT.

conda install bioconda::mafft

mafft --help

The analysis can be done on local computer or on ssh server.

### 3. Biopython for sequence analysis (comparison)



