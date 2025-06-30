# BLCA-subtype-evolution-paper
## Keywords

`Bladder cancer` · `Tumor heterogeneity` · `Clonal evolution` · `Tumor immune microenvironment` `cell-free DNA` · `WGS` · `Phylogenetics` · `Genomic profiling`· `Transcriptomics`

---

## 🔍 Summary

This repository contains code, data, and results for analyzing clonal and genomic heterogeneity in metastatic bladder cancer. We integrate multi-tumors DNA(WGS/WES) sequencing and cell-free DNA to reconstruct subclonal architectures and infer patterns of tumor evolution across primary and metastatic sites.

---

## 📁 Repository Structure

```bash
├── data/                # Raw and processed sequencing data includes source data files
├── scripts/             # Analysis and visualization code for all downstream analysis. Each type if analysis has its own folder with a name that describes the analysis. 
├── results/             # Output figures and tables. Along with the source data files used to generate figures in main figures.
├── pipelines/           # Used to share pipelines used for the project such as mutation callers and a shell scripts if a custom pipleline was utilized. 
├── LICENSE
└── README.md
