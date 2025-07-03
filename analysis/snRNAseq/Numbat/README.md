# Numbat Pipeline

##  User Guide  
[https://kharchenkolab.github.io/numbat/](https://kharchenkolab.github.io/numbat/)

---

##  Steps to Run Numbat

### 1. Pileup and Phasing  
**`numbat_pileup_example.sh`**  
Example script for running the first step: pileup and phasing.

### 2. CNA Calling  
**`Numbat_Step2.R`**  
Run Numbat CNA calling using the output from the pileup step.

---

##  Numbat Analysis Scripts

- **`Numbat_analysis.Rmd`** – Graph tumor vs. normal clones  
- **`Numbat_GeneSpecific_makeDFs.Rmd`** – Create dataframes of gene-level clonality status  
- **`Numbat_GeneSpecific.Rmd`** – Visualize and analyze gene-specific clonality using the above dataframes

---

##  Preparing TITAN Bulk CNA Calls for Numbat Input

### 1. Convert TITAN Output  
**`Numbat_TITAN_CNV_to_Numbat_input.Rmd`**  
Convert TITAN `segs.txt` file to Numbat input format.

### 2. Combine CNA Profiles  
**`Numbat_combine_segs.Rmd`**  
Create a "union" CNA profile from multiple tumors after conversion.
