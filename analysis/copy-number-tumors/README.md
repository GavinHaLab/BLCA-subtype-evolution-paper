# CNA matrix generation and other analysis in Fig. 3 and Supplements

This repository contains scripts for generating the CNA matrix from a copy number caller and downstream analysis.
## Contents
- **`makeMatrixFromTITAN-ICHOR_segBased_git.R`** –R script to generate a matrix for copy number profiles.

- **`CNA_analysis.Rmd** – CNA analysis and Supplementary Fig. 4 CNA related analysis, processing, cleaning and visualization script.

  ## Usage
1.  Rscript makeMatrixFromTITAN-ICHOR_segBased_hg38.R /path/to/TitanSegFiles/ GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt sampleList.txt common Gene 1000 output_directory.
2.  Where sampleList.txt is a list of sample names that should match with TITAN output file names, commom, Gene and 1000 are other flags for the script.
3. To run from the command line for the Rmd script: R -e "rmarkdown::render('CNA_analysis.Rmd')" or run the Rmd scripts in a RStudio.

