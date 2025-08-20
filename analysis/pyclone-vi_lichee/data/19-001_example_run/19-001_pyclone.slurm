#!/bin/bash
#SBATCH --job-name=pyclonevi_19-001
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --output=19-001/logs/19-001.20250820_154946.out
#SBATCH --error=/19-001.20250820_154946.err

# Load all the modules and test on local login nodes before runnin
ml purge
ml Python
conda activate pyclone-vi # or use ml pyclone -vi

# Make the directory for logs
mkdir -p analysis/pyclone-vi/data/19-001/logs
# Fit method
pyclone-vi fit -i analysis/pyclone-vi/data/19-001/19-001_pyclone_Input.tsv -o analysis/pyclone-vi/data/19-001/19-001.hd5 -c 40 -d beta-binomial -r 10
# Write results to a file
pyclone-vi write-results-file -i analysis/pyclone-vi/data/19-001/19-001.hd5 -o analysis/pyclone-vi/data/19-001/19-001_output.tsv
