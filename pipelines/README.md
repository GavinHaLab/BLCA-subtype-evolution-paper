# Variant Calling Pipelines

This folder contains custom pipelines used for variant calling in this project. Pipelines already published and available on GitHub are referenced in the manuscript Methods section with direct links.

## Contents

- **Strelka2/**  
  Contains the Snakemake workflow and configuration files used to run Strelka2.

- **MuSE.sh**, **Varscan2.sh**  
  Shell scripts to run MuSE and Varscan2 for a single tumor-normal pair. These were executed via wrappers to parallelize submissions on an HPC cluster.

## Notes

- Each script is self-contained with usage instructions in the header.
- Paths and parameters may need adjustment based on your environment or cluster setup.
