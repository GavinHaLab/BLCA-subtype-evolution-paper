# BLCA-subtype-evolution-paper
## Keywords

`Metastatic Bladder cancer` · `Tumor heterogeneity` · `Clonal evolution` · `Tumor immune microenvironment` `cell-free DNA` · `WGS` · `Phylogenetics` · `Transcriptomics`. `Genomic profiling`·

---

## 🔍 Summary

This repository contains code, data, and results associated with a manuscript investigating clonal and genomic heterogeneity in metastatic bladder cancer. We integrate multi-tumor DNA (WGS/WES) sequencing and cell-free DNA to reconstruct subclonal architectures and infer tumor evolutionary patterns across primary and metastatic sites.

**<u>The evolution and heterogeneity of lethal metastatic bladder cancer subtypes</u>**

Pushpa Itagi1*, Samantha L. Schuster*, Sonali Arora, Thomas W. Persse, Jennifer A. Waters, Michael Yang, Alan Min, Pooja Chandra, Mohamed Adil1, Patricia C. Galipeau1, Dmytro Rudoy, Allie S. Kreitman, Yixin Lin, Minjeong Ko, Erolcan Sayar, Robert D. Patton, Lori Kollath, Abby Meis, Samuel Lindergren, Nathan Ji, Khursheed Ali, Hrishi Venkatesh, A. Patrick McDeed, Claire B. Mills, Manasvita Vashisth, Cynthia L. Wladyka, Rosa Nadal, Jessica Hawley, Todd Yezefski, Sarah P. Psutka, John L. Gore, Daniel W. Lin1, Peter S. Nelson, Heather H. Cheng, Michael T. Schweizer, Lawrence Fong6, John K. Lee, Evan Yu, Eva Corey, Colm Morrissey, Petros Grivas, Robert B. Montgomery, Jonathan L. Wright, Michael C. Haffner, Funda Vakar-Lopez, Omar Y. Mian†, Hung-Ming Lam†, Andrew C. Hsieh†, Gavin Ha†
###



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
