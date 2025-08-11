# Variant Calling Pipelines
All pipelines in the paper use GRCh38 version.  
This folder contains custom pipelines used for variant calling in this project. Pipelines already published and available on GitHub are referenced in the manuscript Methods section with direct links.

## Contents
- **ichorCNA/**  
   Contains the Snakemake workflow and configuration files used to run ichorCNA. The orginal repo can be found at https://github.com/GavinHaLab/ichorCNA

- **mutect2/**  
  The orginal repo can be found at https://github.com/GavinHaLab/ichorCNA](https://github.com/GavinHaLab/mutect2_snakemake. The same pipeline was used with new Panel of Normals(PONs). These were made for WES and WGS samples from the matched normals of the cohort to remove sequencing artifacts.  
  
- **Strelka2/**  
   Contains the Snakemake workflow and configuration files used to run Strelka2.

- **MuSE.sh**, **Varscan2.sh**  
  Shell scripts to run MuSE and Varscan2 for a single tumor-normal pair. These were executed via wrappers to parallelize submissions on an HPC cluster.

- **TITAN-SVABA/**  
  The orginal repo can be found at https://github.com/GavinHaLab/TitanCNA_SV_WGS. The same pipeline was used with adjusting configs for bins (WES = 50kb, and WGS = 10kb) and was run with the corresponding matched normal and UCSC style for chromosomes. 

## Notes

- Each script is self-contained with usage instructions in the header.
- Paths and parameters may need adjustment based on your environment or cluster setup.
