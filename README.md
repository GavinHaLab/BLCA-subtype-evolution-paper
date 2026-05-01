# BLCA-subtype-evolution-paper
## Keywords

`Metastatic Bladder cancer` · `Tumor heterogeneity` · `Clonal evolution` · `Tumor immune microenvironment` `cell-free DNA` · `WGS` · `Phylogenetics` · `Transcriptomics`. `Genomic profiling`·

---

## 🔍 Summary

This repository contains code, data, and results associated with a manuscript investigating clonal and genomic heterogeneity in metastatic bladder cancer. We integrate multi-tumor DNA (WGS/WES) sequencing and cell-free DNA to reconstruct subclonal architectures and infer tumor evolutionary patterns across primary and metastatic sites.

The evolution and heterogeneity of lethal metastatic bladder cancer subtypes

Pushpa Itagi1,2*, Samantha L. Schuster1,2*, Sonali Arora2, Thomas W. Persse1, Jennifer A. Waters2, Michael Yang1, Alan Min1, Pooja Chandra1, Mohamed Adil1,2, Patricia C. Galipeau1,2, Dmytro Rudoy2, Allie S. Kreitman3, Yixin Lin4, Minjeong Ko1, Erolcan Sayar2, Robert D. Patton1,2, Lori Kollath5, Abby Meis5, Samuel Lindergren5, Nathan Ji5, Khursheed Ali5, Hrishi Venkatesh6, A. Patrick McDeed1,2, Claire B. Mills2,7, Manasvita Vashisth1,2, Cynthia L. Wladyka2, Rosa Nadal8,9, Jessica Hawley8,9, Todd Yezefski8, Sarah P. Psutka5, John L. Gore5, Daniel W. Lin1,5, Peter S. Nelson1,2,3,5,9, Heather H. Cheng8,9, Michael T. Schweizer8,9, Lawrence Fong6, John K. Lee2,10, Evan Yu8,9, Eva Corey5, Colm Morrissey5, Petros Grivas8,9, 

---

## 📁 Repository Structure

```bash
├── data/                # Raw and processed sequencing data includes source data files
├── scripts/             # Analysis and visualization code for all downstream analysis. Each type if analysis has its own folder with a name that describes the analysis. 
├── results/             # Output figures and tables. Along with the source data files used to generate figures in main figures.
├── pipelines/           # Used to share pipelines used for the project such as mutation callers and a shell scripts if a custom pipleline was utilized. 
├── LICENSE
└── README.md

## Resources
The GRCh38 version used for this project is downloaded from the google bucket gs://genomics-public-data/resources/broad/hg38/v0/  
