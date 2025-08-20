# Pyclone-VI and LICHeE Analysis

This directory contains custom scripts for generating, running pyclone-VI and LICHeE tools. 
Original tool github : https://github.com/Roth-Lab/pyclone-vi and https://github.com/viq854/lichee/tree/master

## Overview
- Steps to run Pyclone-Vi and LICHeE
- Scripts used in the analysis which are avaiable for reference/use

## Key Scripts
- `masterScript_Pyclone_git.py`: Used for generating the inputs for pylcone-VI.
- `pycloneToLichee.py`: Converts the pyclone-VI to LICHeE require format. 


## Overview of steps
1. For each sample/tumor get a consensus list of SNVs and INDELs (high confidence of mutations).
2. Next, make a Union list of all mutations across multiple-samples from a patient (if there is multi samples data).
3. Run collect allelic counts to get the reads for mutations, run each sample with this Union list made in step2. See Commands below for further details. 
4. Place the results from Step3 and curated copy number files (we used TITAN) in two folders.
5. Use the masterScript_Pyclone_git.py script and place the absoulte paths to above folders in the variables allelicCounts_path and CNA_dir (with segment/bin-level copy number files and tumor fraction which are in the .params file).
6. Get all the imputs reads such as samples file with a list of sample and patient IDs, Indel folders, output directory, patient clusters (default Pyclone-VI takes 40) if there are any, and a Yes or No to run LICHeE (if we want only Pyclone-Vi to be run)
7. Run the masterScript_Pyclone_git.py. See Commands below for further details. This should produce a folder and a pycloneInput.tsv file. The demo has sucessful runs example.
8. The above script directly submits a job to an HPC cluster if ran successfully and also provides a generates run_lichee.sh script and outputs in the patient folder.
9. The LICHeE requires another script pycloneToLichee.py which generates a sSNV.txt file that is required for running LICHeE. Set the right sizes and error thresholds as needed, see Methods.
10. Use FishPlot or pairtee if needed to use as an alternate for LICHeE.

   
##Resources
## Example: Running of Pyclone-Vi and LICHeE on Sample `19-001`
This repository includes an example folder, `19-001/`, containing output files generated using the script `masterScript_Pyclone_git.py`.
Input files and expected results are provided in the `19-001_example/` directory.

### Workflow for clonality Analysis

**Commands and usage of scripts**
```bash
**GATK function to get read counts**
gatk CollectAllelicCounts -I sample.bam -L patient.intervals -o sample_allelicCounts_output.tsv
**Script to generate and run Pyclone-Vi, this also includes the sbatch script to submit the slurm job to HPC**
python3 masterScript_Pyclone_git.py --samples_file  sampleList_demo.txt --patient_clusters 40 --indels "No" --indels_Dir "" --runLichee "Yes" --output_dir ./data/sample_name/
**LICHeE**
python3 pycloneToLichee.py /path/to/Pyclone_output/PatientID/PatientIDoutput.tsv ./data/PatientID/ PatientID
tools/lichee/LICHeE/release/lichee -build -i ./data/PatientID/PatientID_lichee_sSNV.txt -cp -sampleProfile -minRobustNodeSupport 2 -minClusterSize 2 -maxClusterDist 0.2 -minPrivateClusterSize 1 -e 0.1 -o ../data/PatientID/PatientID_lichee_trees.txt -dotFile  ./data/PatientID/PatientID_lichee_tree.dot -color -dot
**DOT to PDF Conversion**
This output gives a pdf for the best scoring trees, if you need all trees or wish  to remove certain nodes/clusters from Pyclobe-VI then choose the GUI mode to get the trees.  https://github.com/viq854/lichee/tree/master and use the -showTree Flag. 
dot -Tpdf ./data/PatientID/PatientID_lichee_tree.dot -O

