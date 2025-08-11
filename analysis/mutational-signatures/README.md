# Mutational Signature Analysis

This repository contains scripts used for SBS mutaional signature analysis used in the paper (Fig. 4)
- **`run_SigProfiler_Assignment_sampleBasis_git.py`** – Describes the run for SBS mutational signatures. This file requires a input folder where each cluster from pyclone is saved a a .vcf.
- Uses COSMIC V3.2 signature and a custom signature database (excludes UV signatures and non platinum chemotherapy signatures)
- **`UWTAN_Mutation_Signatures.Rmd`** – R markdown script contains the code used for plotting and analysis of mutational signatues such as APOBEC, Clock-like, and platinum chemotherapy. 
- **`Hartwig_Cis_FA.Rmd`** – R markdown script contains the code for plotting cisplatin signatures from Hartwig medical foundation downloaded data.

## Usage
1. python3 run_SigProfiler_Assignment_sampleBasis_git.py
2. To run from the command line:
R -e "rmarkdown::render('UWTAN_Mutation_Signatures.Rmd')"
R -e "rmarkdown::render('Hartwig_Cis_FA.Rmd')"
or run the Rmd scripts in a RStudio

