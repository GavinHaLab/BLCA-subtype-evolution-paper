# Comut and other analysis in Fig. 3

This repository contains scripts for making the landscape comut/oncoprint which has pathogenic Mutations, CNAs and SVs

## Contents
- **`comut_timing_withCNA_git.py`** – Comut to make a landscape plot.
- **`drivergene_manual_annotations_finalBLCAPaperJune2025.csv** – Driver genes with annotation used for analysis. Refer to the Category columb to identify tumor suppressors(TSP) and oncogenes, Other and None are the genes that have no defined roles for TSP or Oncogenes.
- **`fig2_Boxplots_Hist_git.py`** – This script uses different data files generated after curation, cleaning and identifying pathogenic alterations. Different functions in this script refer to various plots used in Fig.3 such as TMB/FGA/SVs by histology.
sampleNameRename_UpsetComut.py
- **`sampleNameRename_UpsetComut.py`** – This script uses a json kind of object to rename the sample IDs in a better or more simpler format. This was used for test purposes and wasn't used in the final paper. 

