#!/usr/bin/python
__author__="Pushpa Itagi"
__purpose__=" Run annovar on these files"
_comments_="22May2023 script to process the strelka2 vcf output and annotate them by ANNOVAR *******"

#Load modules
import os
import sys
import argparse
import pandas as pd

#please load the module annovar if using this script directly by typing : ml annovar on the shell
#******* Change this if needed *******
#Refers to the file extensions that are created after running mutation callers
file_pattern1 = ".vcf.gz"
input_path_to_humandb = "/fh/fast/ha_g/grp/reference/annovar/humandb/HG38/" #databases for annovar files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process files.')
    parser.add_argument('--input_vcf_file_path',help='please enter a path to the vcf file',required=True)
    args = parser.parse_args()
    input_vcf_file_path = args.input_vcf_file_path
    out_file = input_vcf_file_path.split(file_pattern1)[0]
    
    #Add the GT tag to the annotated files so that annovar works on them
    out_file_GT  = out_file + "_GT.vcf"
    df = pd.read_csv(input_vcf_file_path,sep="\t",comment="#",compression='gzip',names =["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR"]) #format reading the DF such that we skip additional filter lines
    #add the required GT tag for annovar to work
    df["FORMAT"] = "GT:"+df["FORMAT"]
    df["NORMAL"] = "0/0:"+df["NORMAL"]
    df["TUMOR"] = "0/1:"+df["TUMOR"]
    print(df.head()," ",out_file_GT)
    df.to_csv(out_file_GT,sep="\t",index=None)

    #******* Change this if needed *******
    command1 = "table_annovar.pl "+out_file_GT+" "+input_path_to_humandb+" -buildver hg38 -out "+out_file+" -remove -protocol refGene,cytoBand,esp6500siv2_all,avsnp144,avsnp150,ALL.sites.2015_08,gnomad_genome,gnomad_exome,exac03,clinvar_20190305,dbnsfp33a,gnomad312_genome,cosmic70 -operation g,r,f,f,f,f,f,f,f,f,f,f,f  -nastring . -polish --vcfinput"
    print(command1)
    print(" Annotated file for input file -",input_vcf_file_path," is here: ",out_file)
    print("\n")
    os.system(command1) #execute the command
