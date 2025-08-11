# Structural Variant Clonality from SVClone Analysis  Demo

This repository contains a small demo for merging structural variants (SVs) across multiple tumor samples from the same patient.

## Contents
- **`TwoOverlap_SVClone_output/`** – Example SVClone output for 4 samples.
- **`SVCloneFix_git.R`** – Generates per-sample SV summary files.
- **`MergeSVs_500bp.R`** – Merges SVs across tumors within a patient using a 500 bp breakpoint overlap.

## Usage
1. Generate summaries:  
   ```bash
   Rscript SVClR <sample_folder>
