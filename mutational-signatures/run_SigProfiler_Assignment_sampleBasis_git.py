
#!/usr/bin/python
# ===============================================================
# Script       : run_SigProfiler_Assignment_sampleBasis_git.py
# Purpose      : Script to generate SBS and DBS based mutational signatures for samples, each sample is provided as a vcf
# Usage        : Run on terminal as python3 run_SigProfiler_extractorAssignment_sampleBasis.py --input_dir ../data/all_samples/ --project_name WGS_BLCA_TAN --exome False
# Dependencies : SigProfilerMatrixGenerator, SigProfilerExtractor, SigProfilerAssignment,pandas
# Author       : Pushpa Itagi
# Note:
#   - This script uses a input parameters that can be modified:
#       * input_file directories    – path to input data where mutation data per sample as a .vcf file is stored where each column has details such as chromosome, start, end, reference, alternate tab separated values
#       * output_file directories   – path to output data where results will be saved 



print("Loading required modules.....")
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
#import SigProfilerExtractor
from SigProfilerExtractor import sigpro as sig
import pandas as pd
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
import sys
import os
print(sys.path)
import argparse


# Function to call SBS Mutational Signatures
def sigProfilerAssignmentSBS(input_file_dir,output_dir,project_name,exome_parameter):
    #set directories and paths to signatures and samples
    print(" In sig profiler sigProfilerAssignment SBS Function",exome_parameter)

    dir_inp     = './data/'
    samples_file     = input_file_dir+"output/SBS/"+project_name+".SBS96.all"
    samples_file = input_file_dir
    output_dir = input_file_dir+"/SBS/"
    ###output_dir      =  dir_inp+"19-022/"
    #signatures  = spa.__path__[0]+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
    sigs        = "COSMIC_v3_SBS_GRCh37_noSBS84-85.txt" #Custom Signature Database, IF NEEDED ^^^^^
    #Analysis of SP Assignment
    Analyze.cosmic_fit( samples_file,
                    output_dir,
                    genome_build="GRCh38",
                    cosmic_version=cosmic_version,
                    verbose=True,
                    context_type="96",
                    input_type="vcf",
                    collapse_to_SBS96=True,
                    signature_database="../data/Assignment_Solution_Signatures_UVart_odd.txt",
                    collapse_to_SBS96=True,
                    make_plots=True,
                    sample_reconstruction_plots=True,
                    exclude_signature_subgroups=['UV_signatures'],
                    exome=exome_parameter,export_probabilities=True)

# Function to call DBS Mutational Signatures
def sigProfilerAssignmentDBS(input_file_dir,output_dir,project_name,exome_parameter):
    #set directories and paths to signatures and samples
    print(" In sig profiler sigProfilerAssignment DBS  Function")
    ###dir_inp     = '../data/'
    ###samples_file     = input_file_dir+"output/DBS/"+project_name+".DBS78.all"
    ###output_dir      =  dir_inp+"19-037/DBS/"
    samples_file = input_file_dir
    output_dir = input_file_dir+"/DBS/"
    #Analysis of SP Assignment
    Analyze.cosmic_fit( samples_file,
                    output_dir,
                    genome_build="GRCh38",
                    cosmic_version=cosmic_version,
                    verbose=True,
                    context_type="DINUC",
                    input_type="vcf",
                    make_plots=True,
                    collapse_to_SBS96 = False,
                    exclude_signature_subgroups=None,
                    sample_reconstruction_plots=True,
                    exome=exome_parameter,export_probabilities=True)


if __name__=="__main__":
    cosmic_version = 3.2
    # Create command-line argument parser
    parser = argparse.ArgumentParser()
    # Add positional argument
    parser.add_argument('--input_dir',required=True,help="name of the input directory containing file(s)")
    parser.add_argument('--project_name',required=True,help="name of the project")
    parser.add_argument('--exome',required=True,help="True if the data is from exomes else Fals")
    # Parse arguments from terminal
    args = parser.parse_args()
    # Access the arguments
    input_file_dir = args.input_dir
    project_name = args.project_name
    exome_parameter = args.exome
    if(exome_parameter == "True"):
        exome_parameter = True
    else:
        exome_parameter = False
    print("processing files in : ",input_file_dir," ",project_name)

    sigProfilerAssignmentSBS(input_file_dir,input_file_dir,project_name,exome_parameter)
    ####sigProfilerAssignmentDBS(input_file_dir,input_file_dir,project_name,exome_parameter)
