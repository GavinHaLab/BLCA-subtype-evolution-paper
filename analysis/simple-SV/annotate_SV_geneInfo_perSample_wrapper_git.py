#!/usr/bin/python
# ===============================================================
# Script       : annotate_SV_geneInfo_perSample_wrapper
# Purpose      : Wrapper script to generate annotated SV files with gene information for each sample
# Usage        : Run on terminal as python3 annotate_SV_geneInfo_perSample_wrapper.py --input_sv_file sampleName.mergedTitan.sv.bedpe --sample_id sampleName
# Dependencies : pandas, glob
# Author       : Pushpa Itagi
# Note:
#   - This script uses a input parameters that can be modified:
#       * input_files    – path to input data where mutation data per sample as a .vcf file is stored where each column has details such as chromosome, start, end, reference, alternate tab separated values
#       * output_files   – a file annotate_SV_TwoOverlap_geneInfo_perSample.sh that contains command to run the script for each sample annotate_SV_geneInfo_perSample_git.py

import glob
import pandas as pd

# Specify the directory where the .bedpe files are located post using the CombiningSVCallers_TwoOverlap_git.R script
path = "../../results/TwoOverlap/1kbgap_annotated_SV/" # ******* CHANGE AS NEEDED *******

# Get a list of all .vcf files in the directory and its subdirectories
files = glob.glob(path + '**/*.mergedTitan.sv.bedpe', recursive=True) #match with the required pattern ******* CHANGE AS NEEDED *******

# Main function to iterate over the files
if __name__ == '__main__':
     f1 = open("annotate_SV_TwoOverlap_geneInfo_perSample.sh", "w") # output file name ******* CHANGE AS NEEDED *******
     script_name = "annotate_SV_geneInfo_perSample_git.py" # code used to annotate SV
     # For each file in the list of files
     for file in files:
          # Print the file path
          print(file)
          sample_id = file.split("/")[-1].split(".")[0] # extracts the sample id
          string_out = " python " + script_name + " --input_sv_file " + file + " --sample_id " + sample_id + "\n"
          print(string_out)
          f1.write(string_out)

     f1.close()
