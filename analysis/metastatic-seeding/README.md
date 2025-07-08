# Metastatic BLCA Seeding Analysis

This directory contains custom scripts for analyzing metastatic seeding patterns in bladder cancer. One script is used to run the seeding tool, and another processes the output for downstream analysis. Usage details are included within each script.

## Overview
- Seeding analysis generates inputs converting and formatting LiCHeE outputs to run the MACHINA tool. 
- Outputs include clone-level phylogenies and summary annotations per patient inferred using migrations from MACHINA and CPs from Pyclone-VI.

## Key Scripts
- `clone_tree_generation_for_machina_BLCA.R`: Runs the seeding analysis per patient using phylogenetic input data from LiCHEE. Generates three tsv outputs - colorfile, labelling and tree files. 
- `Metastatic_seeding_analysis.Rmd`: Parses and summarizes seeding results for integration with clinical and genomic data.

