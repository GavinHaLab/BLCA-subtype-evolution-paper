# Metastatic BLCA Seeding Analysis

This directory contains custom scripts for analyzing metastatic seeding patterns in bladder cancer. One script is used to run the seeding tool, and another processes the output for downstream analysis. Usage details are included within each script.

## Overview

- Seeding analysis generates inputs to run the MACHINA tool. 
- Outputs include clone-level phylogenies and summary annotations per patient.

## Key Scripts

- `run_seeding_tool_git.R`: Runs the seeding analysis per patient using phylogenetic input data from LiCHEE. 
- `process_seeding_output_git.R`: Parses and summarizes seeding results for integration with clinical and genomic data.

## Resources

- Clone assignments and phylogenies are assumed to be generated using LICHeE or PyClone-VI.
- Output from this analysis feeds into survival and treatment correlation analyses described in .
