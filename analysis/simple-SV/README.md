# Metastatic BLCA  SV Analysis

This repository contains custom scripts used for generating consensus structural variants (SVs) and BEDPE files for bladder cancer evolution analysis.
Scripts also include breakpoint annotation using GENCODE Release 44 (GRCh38.p14).
Usage instructions are included within each script.

## Contents

- `CombiningSVCallers_TwoOverlap_git.R`: Merges and formats SVs across 3 callers for each sample
- `annotate_sv_genes.py`: Annotates SV breakpoints with intersecting genes
