#!/usr/bin/env python3
''' This is a important wrapper and preprocessing script runs Pyclone-VI after doing preprocessing steps and until obtaining the clusters file for inputs
Step1: Filter allelic Counts - allelic Counts should be ready before this script runs, this step makes unique allelic counts ready for each sample/group of related samples.
Step2: Extract Copy Number Information - This step extracts the copy number information from TITAN or your chosen copy number caller segment files and combines it with allelic counts. and extracts details such as tumor fraction.
Step3: Run Pyclone-VI - This step runs Pyclone-VI and generates the output files. Sbatch the files. 
Example Run: python3 masterScript_Pyclone_git.py --samples_file  sampleList_demo.txt --patient_clusters 40 --indels "No" --indels_Dir "" --runLichee "Yes" --output_dir ./data/sample_name/
Author : Pushpa Itagi
'''
####awk -F'\t' 'NR>1 {print $4":"$5"-"$5"\t"$3}' OFS="\t" C1_C2_822_union_final_2.tsv 

# Load Modules
import argparse
import os
import sys
from pathlib import Path
import pandas as pd
import subprocess
import pyranges as pr


# Get the list of all files and directories that have allelic counts per sample
allelicCounts_path = "/path/to/allelicCounts/"  # Change this to the path where allelic counts files are stored
# Get the list of all Copy number solutions from TITAN (Curated) one for each sample in the sample list
CNA_dir = "/path/to/CNAsolutions/"

import os, shlex, tempfile, subprocess
from datetime import datetime


def sbatchPyclone(pyclone_input_file,output_dir,patient,patient_clusters,runLichee,indels_folder_path,with_indels,time_limit="08:00:00",mem_gb=4000, cpus=1,partition="campus-new", account=None):

    in_file = os.path.abspath(pyclone_input_file)
    print("allegCounts_path:", output_dir)
    out_dir = os.path.abspath(output_dir)
    sample_out = out_dir+"/"+patient+".hd5"
    sample_out2 = out_dir+"/"+patient+"_output.tsv"
    os.makedirs((out_dir+"/"), exist_ok=True)
    print("patient clusters:", patient_clusters)
    log_dir = os.path.join((out_dir+"/"), "logs")
    os.makedirs(log_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    stdout_log = os.path.join(log_dir, f"{patient}.{ts}.out")
    stderr_log = os.path.join(log_dir, f"{patient}.{ts}.err")

    pyclone_cmd = f"pyclone-vi fit -i {shlex.quote(in_file)} -o {shlex.quote(sample_out)} -c {patient_clusters} -d beta-binomial -r 10".strip()
    pyclone_write = f"pyclone-vi write-results-file -i {shlex.quote(sample_out)} -o {shlex.quote(sample_out2)}".strip()

    sbatch_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=pyclonevi_{patient}",
        f"#SBATCH --time=08:00:00",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --nodes=1",
        f"#SBATCH --mem={mem_gb}",
        f"#SBATCH --output={shlex.quote(stdout_log)}",
        f"#SBATCH --error={shlex.quote(stderr_log)}",
    ]

    sbatch_lines += [
    "ml purge", # clear the modules
    "ml Python", # clear the modules
    "conda activate pyclone-vi", # Load the PyClone-VI module, changed this as cluster runs were throwing errors for openMPI. So, test and run accordingly, 
    f"mkdir -p {shlex.quote(log_dir)}",] + [pyclone_cmd] + [pyclone_write]
    script_path = out_dir +"/"+ patient+ "_pyclone.slurm"
    with open(script_path, "w") as fh:
        fh.write("\n".join(sbatch_lines) + "\n")

    print(f"Generated sbatch script: {script_path}")
    # Use os module to submit the sbatch script
    submit_command = f"sbatch {shlex.quote(script_path)}"
    os.system(submit_command)
    print(f"Submitted sbatch script: {submit_command}")
    return sample_out2;


# Function to extract allelic counts
def extractAllelicCounts(allelicCounts_path,patient_list_file):
    print("Extracting Allelic Counts from ", allelicCounts_path)

    # Step 1: Read the patient list file
    patient_list = pd.read_csv(patient_list_file, sep='\t', header=None, names=['Patient', 'Sample'])
    print("Patient List:\n", patient_list.head())
    
    # Step 2: Group samples by patient
    patient_samples = patient_list.groupby('Sample')['Patient'].apply(list).to_dict()
    print("Grouped Samples by Patient:\n", patient_samples)
    combined_data = []
    # Step 3: Loop through patients and their samples
    for patient, samples in patient_samples.items():
        print(f"Processing Patient: {patient}")
        for sample in samples:
            file_name = f"{sample}_allelicCounts.tsv"
            file_path = os.path.join(allelicCounts_path, file_name)
            print(f"Looking for file: {file_path}")
            # Create a new DataFrame for the current patient
            df_new = pd.DataFrame(columns=["mutation_id", "desc", "Sample"])
            if os.path.exists(file_path):
                print(f"Reading file: {file_path}")
                df = pd.read_csv(file_path, sep='\t', comment='@')
                
                 # Create new columns for each row in the DataFrame
                df_new["mutation_id"] = df.apply(lambda row: f"{row['CONTIG']}:{row['POSITION']}-{row['POSITION']}", axis=1)
                df_new["desc"] = df.apply(lambda row: f"{row['CONTIG']}:{row['POSITION']}_{row['REF_NUCLEOTIDE']}_{row['ALT_NUCLEOTIDE']}", axis=1)
                df_new["Sample"] = sample
                df_new["REF_COUNT"] = df["REF_COUNT"]
                df_new["ALT_COUNT"] = df["ALT_COUNT"]
                combined_data.append(df_new)
            else:
                print(f"File not found: {file_path}")
    
        # Step 4: Combine all data into a single DataFrame
        if combined_data:
            final_df = pd.concat(combined_data, ignore_index=True)
            out_file = os.path.join(allelicCounts_path, patient + '_combined_allelicCounts.tsv')
            final_df.to_csv(out_file, sep='\t', index=False)
            print(f"Combined allelic counts saved to: {out_file}")
            
        else:
            print("No data to combine.")
        # clear combined_data
        combined_data.clear()

    return allelicCounts_path;

# Write a function to extract copy number information and add it to the allelic counts file
def copyNumberAddition(allelicCounts_path,patient_list_file,patient_clusters,output_dir):
    print("Extracting Copy Number Information from ", CNA_dir)
    # read all the files ending with _combined_allelicCounts.tsv in allelicCounts_path
    allelic_counts_files = [f for f in os.listdir(allelicCounts_path) if f.endswith('_combined_allelicCounts.tsv')]
    print("Allelic Counts Files:\n", allelic_counts_files)
    
    for file in allelic_counts_files:
        file_path = os.path.join(allelicCounts_path, file)
        print(f"Processing file: {file_path}")
        patient = os.path.basename(file_path).split('/')[-1].split('_')[0]
        print(f"Patient ID: {patient}")
        # Read the allelic counts file
        df_allelic_counts_all = pd.read_csv(file_path, sep='\t')
        samples_list = df_allelic_counts_all['Sample'].unique()
        print("Samples List:\n", samples_list)
        patient_df = pd.DataFrame()
        for sample in samples_list:
            print(f"Processing Sample: {sample}")
            # open the file in CNA_dir that matches the sample name
            df_allelic_counts = df_allelic_counts_all[df_allelic_counts_all['Sample'] == sample]
            print("Allelic Counts DataFrame for sample:\n", df_allelic_counts.head())
           
            # find the params.txt file for the sample
            params_filename = next((os.path.join(CNA_dir, fn) for fn in os.listdir(CNA_dir) if fn.startswith(sample) and fn.endswith("params.txt")),None)
            print("test params_filename:", params_filename)
            if params_filename:
                normal = subprocess.check_output(
                    f"awk -F'\t' '/Normal contamination estimate:/ {{print $2}}' {params_filename}",
                    shell=True, text=True
                ).strip()
                tumor_fraction = round(1 - float(normal), 3)
                print(f"Normal contamination estimate for sample {sample}: {normal}, Tumor fraction: {tumor_fraction}")

            # Next handle the seg file
            seg_filename = next(
                (os.path.join(CNA_dir, fn) for fn in os.listdir(CNA_dir) if sample in fn and fn.endswith(".titan.ichor.seg.noSNPs.txt")),
                None )
            ###print(f"Segment file for sample {sample}: {seg_filename}")
            if seg_filename:
                df_segments = pd.read_csv(seg_filename, sep='\t')
                # Read the segment file
                df_segments = pd.read_csv(seg_filename, sep='\t')
                ###print("Alelic counts columns:\n", df_allelic_counts.head())
                df_allelic_counts["Chr"] = df_allelic_counts["mutation_id"].apply(lambda x: x.split(":")[0])
                df_allelic_counts["Start"] = df_allelic_counts["mutation_id"].apply(lambda x: int(x.split(":")[1].split("-")[0]))
                df_allelic_counts["End"] = df_allelic_counts["Start"]
                ###print("Alelic counts columns:\n", df_allelic_counts.head())
                ###print("Segment DataFrame:\n", df_segments.head())  
                
                print("Updated Allelic Counts DataFrame with CN_Major and CN_Minor:\n", df_allelic_counts.head())
                # Convert both to PyRanges
                gr_segments = pr.PyRanges(df_segments.rename(columns={"Chromosome": "Chromosome", "Start": "Start", "End": "End"}))
                gr_counts = pr.PyRanges(df_allelic_counts.rename(columns={"Chr": "Chromosome", "Start": "Start", "End": "End"}))  # treat Start as 1-bp interval
                # Join on overlaps
                joined = gr_counts.join(gr_segments)
                # Back to pandas
                df_joined = joined.df
                print("Joined DataFrame:\n", df_joined.head())
                
                # Add a column called Normal with a value 2 except for ChrX where it is 1
                df_joined['normal_cn'] = 2
                df_joined.loc[df_joined['Chromosome'] == 'ChrX', 'normal_cn'] = 1
                df_joined["tumour_content"] = tumor_fraction
                # retain only these columms
                df_joined = df_joined[['mutation_id', 'desc', 'Sample', 'normal_cn', 'tumour_content','Corrected_MajorCN', 'Corrected_MinorCN','REF_COUNT', 'ALT_COUNT']]
                # arrange the columns in this order
                ordered_columns = ['mutation_id', 'desc', 'Sample', 'REF_COUNT', 'ALT_COUNT','normal_cn',  'Corrected_MajorCN', 'Corrected_MinorCN','tumour_content']
                df_joined = df_joined[ordered_columns]
                # rename the columns to match Pyclone-VI input
                df_joined.rename(columns={'mutation_id': 'mutation_id',
                                            'desc': 'desc',
                                            'Sample': 'sample_id',
                                            'REF_COUNT': 'ref_counts',
                                            'ALT_COUNT': 'alt_counts',
                                            'normal_cn': 'normal_cn',
                                            'Corrected_MajorCN': 'major_cn',
                                            'Corrected_MinorCN': 'minor_cn',
                                            'tumour_content': 'tumour_content'}, inplace=True)
                #concatentate the sample df to a major DF
                patient_df = pd.concat([patient_df, df_joined], ignore_index=True)
                print("Patient DataFrame:\n", patient_df.head())
                ###tumor_fraction = 0.0  # reset tumor_fraction for the next sample
        # write the joined DataFrame to a new file
        pyclone_input_file = os.path.join(output_dir, patient + '_pyclone_Input.tsv')
        
        # drop rows where major_cn is 0
        patient_df = patient_df[patient_df['major_cn'] != 0]
        # fill All NAs with 0
        patient_df.fillna(0, inplace=True)
        # Convert norma_cn, major_cn and minor_cn to integers
        cols = ["normal_cn", "major_cn","minor_cn"]  # columns you want
        patient_df[cols] = patient_df[cols].astype(int) 
        patient_df.to_csv(pyclone_input_file, sep='\t', index=False)
        print(f"PyClone-VI input file saved to: {pyclone_input_file}")
        # Call a function to sbatch Pyclone-VI
        sbatchPyclone(pyclone_input_file,output_dir, patient,patient_clusters,output_dir,runLichee,indels_folder_path,with_indels)
    return 0;   


# Function to run LiCHEE for visualization of trees and generate sh script  
def runLicheSbatch(pyclone_output, output_dir, runLichee,patient_list_file):
    if runLichee.lower() == "yes":
        # Step 1: Read the patient list file
        patient_list = pd.read_csv(patient_list_file, sep='\t', header=None, names=['Patient', 'Sample'])
        print("Patient List:\n", patient_list.head())
        
        # Step 2: Group samples by patient
        patient_samples = patient_list.groupby('Sample')['Patient'].apply(list).to_dict()
        print("Grouped Samples by Patient:\n", patient_samples)
        combined_data = []
        # Step 3: Loop through patients and their samples
        for patient, samples in patient_samples.items():
                print("Running LiCHEE for visualization of trees")
                # Call the function to run LiCHEE
                # make a shell script 
                lichee_script = os.path.join(output_dir, "run_lichee.sh")
                with open(lichee_script, "w") as f_out:
                    f_out.write("ml Java/1.8\n")
                    f_out.write("python3 pycloneToLichee.py " + pyclone_output +" " + output_dir + " "+patient+"\n")
                    # old method that uses clusters file and is not required but optional
                    ####second_string = "LICHeE/release/lichee -build -i " + output_dir +"/" + patient +"_lichee_sSNV.txt" + " -clustersFile " + output_dir +"/" + patient +"_lichee_cluster.txt" +" -cp -sampleProfile -minRobustNodeSupport 2 -minClusterSize 2 -maxClusterDist 0.2 -minPrivateClusterSize 1 -e 0.1 -o " + output_dir +"/" + patient +"_lichee_trees.txt" + " -dotFile  " + output_dir +"/" + patient +"_lichee_tree.dot" + " -color -dot"
                    second_string = "LICHeE/release/lichee -build -i " + output_dir +"/" + patient +"_lichee_sSNV.txt" +" -cp -sampleProfile -minRobustNodeSupport 2 -minClusterSize 2 -maxClusterDist 0.2 -minPrivateClusterSize 1 -e 0.1 -o " + output_dir +"/" + patient +"_lichee_trees.txt" + " -dotFile  " + output_dir +"/" + patient +"_lichee_tree.dot" + " -color -dot"
                    f_out.write(second_string + "\n")
                    third_string = "dot -Tpdf  " + output_dir +"/" + patient +"_lichee_tree.dot" + " -O"
                    f_out.write(third_string + "\n")
                f_out.close()
                print(f"Generated LiCHEE script: {lichee_script}")
                # Make the script executable
                os.chmod(lichee_script, 0o755)
                print(f"LiCHEE script is ready to run: {lichee_script}")
    return 0;


# Main function to run the script
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--samples_file',help='please enter a file with a list of Sample and Patient IDs specified in a file',required=True) # look for sampleList_demo.txt
    parser.add_argument('--patient_clusters',type=int,help='please enter a limit to clusters for each Patient ID 7',required=True) # enter the number of clusters you ant Pyclone -VI to use, default is 40
    parser.add_argument('--indels', type=str, help='please enter True or False to include indels for this analysis or not, ex: True',required=True) # enter a Yes or No to include indels in the analysis
    parser.add_argument('--indels_Dir', type=str, help='please enter the directory to where Indel files are located [ if you give  --indels as True',required = True) # get the indels directory path such as /path/to/Indels
    parser.add_argument('--runLichee', type=str, help='please enter Yes or No to inlcude making trees from Lichee this analysis or not, ex: True',required=True) # for the visualization of tree use lichee and generate a .sh script to do that
    parser.add_argument('--output_dir', type=str, help='please enter the output directory to save the results',required=True) # specify the output directory where you want to save the results
    
    args = parser.parse_args()

    patient_list_file = args.samples_file
    patient_clusters = args.patient_clusters
    runLichee = args.runLichee
    indels_folder_path = args.indels_Dir
    with_indels = str(args.indels)
    output_dir = args.output_dir
    print("\n Begin Processing Pyclone-VI and LiChee ~~~")
    print("Patients list ",patient_list_file," indels include :",with_indels)

    # a function to extract the read counts for each sample
    allelicCounts_path=extractAllelicCounts(allelicCounts_path,patient_list_file)
    #allelicCounts_path = "./data/"
    # Pass the path for the merged allelc counts file, did a seperate method in case the allelic counts are not ready
    pyclone_output = copyNumberAddition(allelicCounts_path,patient_list_file,patient_clusters,output_dir)
    pyclone_output = "data/19-001/19-001_output.tsv"
    runLicheSbatch(pyclone_output, output_dir, runLichee,patient_list_file) 
   
