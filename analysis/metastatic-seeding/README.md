# Metastatic BLCA Seeding Analysis

This directory contains custom scripts for analyzing metastatic seeding patterns in bladder cancer. One script is used to run the seeding tool, and another processes the output for downstream analysis. Usage details are included within each script.

## Overview
- Seeding analysis generates inputs converting and formatting LiCHeE outputs to run the MACHINA tool. 
- Outputs include clone-level phylogenies and summary annotations per patient inferred using migrations from MACHINA and CPs from Pyclone-VI.

## Key Scripts
- `clone_tree_generation_for_machina_BLCA.R`: Runs the seeding analysis per patient using phylogenetic input data from LiCHEE. Generates three tsv outputs - colorfile, labelling and tree files. 
- `Metastatic_seeding_analysis.Rmd`: Parses and summarizes seeding results for integration with clinical and genomic data.

##Resources
- Added an example folder 19-001 with the resulting files from using clone_tree_generation_for_machina_BLCA.R and the command to run MACHINA along with the expected results file is added in folder 19-001_input_example
--  Pmh
-- pmh -p P001PriFFPE -c 19_001/19_001_colorfile.tsv 19_001/19_001_tree.tsv 19_001/19_001_labeling.tsv -o 19_001/pmh/ > 19_001/pmh/result.txt
-- Pmh_tr
-- pmh_tr -p P001PriFFPE -c19_001/19_001_colorfile.tsv 19_001/19_001_tree.tsv 19_001/19_001_labeling.tsv -o 19_001/pmh_tr/ > 19_001/pmh_tr/result.txt


