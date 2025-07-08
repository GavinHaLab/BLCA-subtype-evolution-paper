# Metastatic BLCA  Simple SV Analysis

This repository contains custom scripts used for generating consensus structural variants (SVs) and BEDPE files for bladder cancer evolution analysis.
Scripts also include breakpoint annotation using GENCODE Release 44 (GRCh38.p14).
Usage instructions are included within each script.

## Contents
- Each SV is retained if it is detetected by atleast 2 out of 3 callers, has size > 1000 base pairs, and is not part of the ENCODE blacklisted regions. 
- `CombiningSVCallers_TwoOverlap_git.R`: Merges and formats SVs across 3 callers for each sample. The callers included in this analyis are SVABA, MANTA and GRIDDS. 
- `annotate_SV_geneInfo_perSample_git.py`: Annotates SV breakpoints with intersecting genes but has to be used with the wrapper annotate_SV_geneInfo_perSample_wrapper_git,py
- `annotate_SV_geneInfo_perSample_wrapper_git.py`: The wrapper that calls the annotate_SV_geneInfo_perSample_git to annotate the breakpoints.
- The GENCODE v44 file is not inlcuded here due to space limitations but can be found here https://www.gencodegenes.org/human/release_44.html
- The driver genes are located here analysis/comut_mutations_cna_sv/drivergene_manual_annotations_finalBLCAPaperJune2025.csv
- The encode BLACKLIST used is uploaded as encodeBlacklistV2_20240927.txt
