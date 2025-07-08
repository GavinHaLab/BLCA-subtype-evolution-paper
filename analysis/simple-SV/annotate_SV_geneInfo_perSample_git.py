#!/usr/bin/python
# ===============================================================
# Script       : annotate_SV_geneInfo_perSample_git.py
# Purpose      : Script to gannotate gene information for the breakpoints based on GENCODE v44 basic annotation file
# Usage        : Run on terminal as 
# Dependencies : pandas, argparse
# Author       : Pushpa Itagi
# Note:
#   - This script uses a input parameters that can be modified:
#       * input_files and directories    – path to BLCA driver genes and gencode basic annotation
#       * output_file directories   – two files for each sample with gene information for the breakpoints and excluding non-protein coding genes



import pandas as pd
import sys
import argparse

# Specify the directory

driver_genes_allCancers = "analysis/comut_mutations_cna_sv/drivergene_manual_annotations_finalBLCAPaperJune2025.csv"
df1 = pd.read_csv(driver_genes_allCancers, sep='\t')

def add_gene_info(df, gene_info_file,sample_name,input_sv_file):

    # Read the gene information file into a DataFrame
    gene_info = pd.read_csv(gene_info_file, sep='\t',header=None, comment='#')
    # Rename columns for clarity
    gene_info.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    # Extract gene name from attribute column
    gene_info['gene_name'] = gene_info['attribute'].str.extract('gene_name "([^"]+)"')
    # Extract gene type from attribute column
    gene_info['gene_type'] = gene_info['attribute'].str.extract('gene_type "([^"]+)"')
    # Convert 'start' and 'end' columns to integers
    gene_info = gene_info[gene_info['start'] != '.']
    gene_info['start'] = gene_info['start'].astype(int)
    print(gene_info['end'])
    # Drop rows where 'end' is '.'
    gene_info = gene_info[gene_info['end'] != '.']
    # Convert 'end' column to integers
    gene_info['end'] = gene_info['end'].astype(int)
    print(gene_info.head(2))
    # Print column 6
    print(gene_info['gene_name'])

    # Initialize gene type columns in df
    df['gene_1_type'] = ''
    df['gene_1_type'] = ''

    # Iterate over each row in df
    for i, row in df.iterrows():
        # Convert row[1] and row[2] to integers
        print("Variant Num.",i)
        start = int(row[1])
        end = int(row[2])
        # Find matching rows in gene_info
        matches_a = gene_info[(gene_info['chr'] == row[0]) & (gene_info['start'] <= start) & (gene_info['end'] >= end)]
        matches_b = gene_info[(gene_info['chr'] == row[3]) & (gene_info['start'] <= row[4]) & (gene_info['end'] >= row[5])]
        #print("matches",matches.head(2))
        # If there are matches, add the gene names to the 'gene' column in df
        if not matches_a.empty:
             df.at[i, 'gene_1'] = ', '.join(matches_a['gene_name'].unique())
             df.at[i, 'gene_1_type'] = ', '.join(matches_a['gene_type'].unique())
             #print("matches a",df.at[i, 'gene_a'])
        else:
             df.at[i, 'gene_1'] = 'NA'
             df.at[i, 'gene_1_type'] = 'NA'
             #print("matches a",df.at[i, 'gene_a'])
        if not matches_b.empty:
             df.at[i, 'gene_2'] = ', '.join(matches_b['gene_name'].unique())
             df.at[i, 'gene_2_type'] = ', '.join(matches_b['gene_type'].unique())
             #print("matches b",df.at[i, 'gene_b'])
        else:
             df.at[i, 'gene_2'] = 'NA'
             df.at[i, 'gene_2_type'] = 'NA'
             #print("matches b",df.at[i, 'gene_b'])
        #print("annotated lines..",df.head())
    output_file2 = input_sv_file.replace('.bedpe', '.bedpe_gene.txt')
    output_file3 = input_sv_file.replace('.bedpe', '.bedpe_proteinCodinggene.txt')
    # Add 'sample_name' column to df
    df['sample_name'] = sample_name
    print(df.head(2))
    # Create a new column 'cancer_driver' in 'output2' and 'output3'
    df['cancer_driver'] = df['gene_1'].isin(df1['SYMBOL']) | df['gene_2'].isin(df1['SYMBOL'])
    df['cancer_driver'] = df['cancer_driver'].replace({True: 'yes', False: 'no'})
    df.to_csv(output_file2, sep='\t', index=False)
    # New line to filter rows containing Protein coding genes
    df = df[(df['gene_1_type'].str.contains('protein_coding')) | (df['gene_2_type'].str.contains('protein_coding'))]
    df.to_csv(output_file3, sep='\t', index=False)
    print("output file written here ....",output_file2)
   
    return df;


# Main function to iterate over the files
if __name__ == '__main__':
     # Create an ArgumentParser instance
     parser = argparse.ArgumentParser(description='Process some files.')
     # Add two arguments
     parser.add_argument('--sample_id', help='The sample id to process.',required=True)
     parser.add_argument('--input_sv_file', help='The SV file to process.',required=True)
     # Parse the arguments
     args = parser.parse_args()
     sample_id = args.sample_id
     input_sv_file = args.input_sv_file

     print("Reading file",input_sv_file)
     output_file = input_sv_file.replace('.vcf', '.bedpe')
     df = pd.read_csv(output_file, sep='\t')
     # Add 'chr' prefix to the 1st and 4th columns
     df['chrom1'] =  df['chrom1'].astype(str)
     df['chrom2'] =  df['chrom2'].astype(str)
     sample_name = sample_id
     print("sample name",sample_name)

     # Add gene information, for HG 38
     gene_info_file = 'gencode.v44.basic.annotation.gtf' # download from the GENOCODE website it is > 25 MB to upload
     # Read the bedpe file into a pandas DataFrame, skipping the first row
     df_return = add_gene_info(df, gene_info_file,sample_name,input_sv_file)
