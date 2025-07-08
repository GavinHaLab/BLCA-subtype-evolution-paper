# Metastatic BLCA Simple SV Analysis

This directory contains custom scripts for generating consensus structural variants (SVs) and BEDPE files for metastatic bladder cancer evolution analysis. Breakpoint annotation is performed using GENCODE Release 44 (GRCh38.p14). Usage details are included within each script.

## Overview

- SVs are retained if detected by ≥2 of 3 callers (SVABA, MANTA, GRIDSS), >1000 bp in size, and outside ENCODE blacklist regions.
- Gene annotation is based on GENCODE v44; blacklist filtering uses `encodeBlacklistV2_20240927.txt`.

## Key Scripts

- `CombiningSVCallers_TwoOverlap_git.R`: Merges SVs across callers per sample.
- `annotate_SV_geneInfo_perSample_git.py`: Annotates breakpoints with intersecting genes.
- `annotate_SV_geneInfo_perSample_wrapper_git.py`: Wrapper to batch-annotate SVs across samples.

## Resources

- GENCODE v44 GTF: [GENCODE Human Release 44](https://www.gencodegenes.org/human/release_44.html). Choosing BASIC gene annottaion
- Driver gene list: `analysis/comut_mutations_cna_sv/drivergene_manual_annotations_finalBLCAPaperJune2025.csv`
- ENCODE blacklist: `encodeBlacklistV2_20240927.txt`
