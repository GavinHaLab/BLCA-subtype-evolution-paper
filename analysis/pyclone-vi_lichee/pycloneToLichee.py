#!/usr/bin/python
'''
Purpose: "Convert Pyclone output file to lichee required inputs"
Usage:
        python3 pycloneToLichee.py /path/toDir/Pyclone_output.tsv /path/toOutputDir/

Author:
        Pushpa Itagi 
'''
#Import modules
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, '/fh/fast/ha_g/projects/Ha_Hsieh_BladderCancer/analysis/important_files_References/') #Refer to the file name with sample IDs which is the analysis folder
from sampleNameRename import *

#Generating header function
def create_header(sample_list):
    full_list = ""
    print("Sample list: ", sample_list)
    # if you do not need to rename files or samples then you can use a function here to rename samples
    for name in sample_list:
        full_list = full_list+"\t"+str(name)
    final_header = "#chr\tposition\tdescription\tprofile\tNormal"+str(full_list)
    return final_header;

#Generating SNV file call from Pyclone Output
def licheeInputGenerationSNV(f,f_out):
    previous=""
    profile="";cp_final="0.0"
    for line in f:
        #print(line)
        str1=line.split("\t")
        if((previous!=str1[0] and len(previous) >0) and (line != last_line1)):
            #New Mutation line
            #cp_final=cp_final+cp+"\t"
            #str_out=str(chrom[0])+"\t"+str(position)+"\t"+"desc1"+"\t"+profile+"\t"+str(cluster)+"\t"+"0.0"+"\t"+cp_final
            str_out=str(chrom[0])+"\t"+str(position)+"\t"+"desc1"+"\t"+"0"+profile+"\t"+cp_final
            f_out.write(str_out)
            f_out.write("\n")
            profile="";cp_final="0.0"

        chrom=str1[0].split(":");mutation_id=str1[0].split("-")
        position=mutation_id[-1]
        sample=str1[1]
        cluster=str1[2]
        cp=str1[3]
        #^^^Certain clusters with a value lower than Cellular Prevalence thrshold are not included
        if(float(cp) >= CP_cutoff): #*******Change this number to remove clusters with lower CP of 0.01, as these can be low read and not true mutations but are counted due to their presence in other samples
            profile=profile+"1"
        else:
            profile=profile+"0"

        cp_final=cp_final+"\t"+cp
        if(line == last_line1): #Appends the very last time correctly
            str_out=str(chrom[0])+"\t"+str(position)+"\t"+"desc1"+"\t"+"0"+profile+"\t"+cp_final
            f_out.write(str_out)
            f_out.write("\n")
        previous=str1[0]
        #cp_final=cp

    f.close()
    return 0;

#Generating Cluster Information from Pyclone
def licheeInputGenerationCluster(f1,f2):
    print("Inside licheeInputGenerationCluster")
    f_out=open(f2,"w")
    df1 = pd.read_csv(f1,sep="\t")
    print(df1)
    num_fo_samples = df.sample_id.nunique()
    ##print(num_fo_samples)
    df_all = pd.DataFrame()
    mutations_cluster = pd.DataFrame()
    df2 = pd.DataFrame()
    for sample in range(0,(num_fo_samples)):
        df2 = df1.groupby("cluster_id").nth(sample) #first group by cluster
        data = [df_all, df2]
        df_all= pd.concat(data,ignore_index=False, sort=False)
    df_all['cluster_num'] = df_all.index #Add index which are basically cluster_ids here into a new column
    df3 = pd.DataFrame()
    ####df3 = df1.groupby(['cluster_id','sample_id'],as_index=False).size() #get the size of each cluster
    df3 = df1.groupby(["cluster_id","sample_id"]).size().reset_index(name="count")
    print("df3",df3)
    ##df3 = df1.groupby(['cluster_id','sample_id'],as_index=False).size().reset_index(name='size') #get the size of each cluster changed for 02/15/2024 as we encountered sample_id errors if using wes+cfdna serum
    df3 = df3.reset_index()
    df3.index.name = 'name'
    ###print(df3.columns)
    one_sample = df3['sample_id'][0] #We need one sample ID to pick up unqiue pairs of clusters-number of mutations
    print("One sample",one_sample)
    mutations_cluster['cluster_id'] = df3.loc[df3['sample_id'] == one_sample, ('cluster_id')] #extract the cluster number/ID and the size of it
    print("Mutations cluster",mutations_cluster)
    #####mutations_cluster['size'] = df3.loc[df3['sample_id'] == one_sample, ('size')] #extract the size of it
    mutations_cluster["size"] = df3.loc[df3["sample_id"] == one_sample, "count"].values #extract the size of it
    ####mutations_cluster['size'] = df3.loc[df3['sample_id'] == one_sample].get('size', default_value)
    #mutations_cluster = df3['size'].unique()
    ####mutations_cluster = df3['size']
    ####mutations_cluster = df3.groupby(['size', 'cluster_id']).ngroups
    ##print("Length.......",(mutations_cluster))

    prev_cluster_count = 0;df_cluster = pd.DataFrame();
    string1 = string2 = string3 = out1 = out2 = "";
    for cluster in range(0,len(mutations_cluster)):
        per_cluster_size = mutations_cluster['size'].iloc[cluster]
        cluster = mutations_cluster['cluster_id'].iloc[cluster]
        print("Cluster ",cluster,"Cluster size ",per_cluster_size,"DF",df_all)
        #index = mutations_cluster.index(str(cluster))
        df_cluster = (df_all.loc[df_all['cluster_num'] == int(cluster)])
        #df_cluster = (df_all.loc[df_all['cluster_num'] == int(index)])
        #df_cluster['cellular_prevalence'] = df_cluster['cellular_prevalence'].astype(float)
        #df_new = df_cluster[df_cluster.cellular_prevalence > CP_cutoff] #Cut offs for CPs
        string1 = (np.where(df_cluster['cellular_prevalence'] > CP_cutoff, df_cluster['cellular_prevalence'], '0.0'))
        string2 = (np.where(df_cluster['cellular_prevalence'] > CP_cutoff, '1', '0'))
        string3 = ""
        ###print(string1, "..",string2)
        out1 = "\t".join([str(item) for item in string1]);out1 = "0.0"+"\t"+out1 #CPs description
        out2 = "".join([str(item) for item in string2]);out2 = "0"+out2 #profile description for lichee
        for mutation_id in range (1,per_cluster_size+1):
            string3 = string3+str(prev_cluster_count+(mutation_id))+","
        out3 = string3[:-1] #remove the last comma appended
        f_out.write(out2+"\t"+out1) #All this is a format required for the cluster file, they need to be tab separated, including the mutations id list
        f_out.write("\t")
        f_out.write(out3)
        f_out.write("\n")
        prev_cluster_count = prev_cluster_count+per_cluster_size
    return 0;

if __name__ == "__main__":
    #Gather all input and output files needed for this script
    pyclone_outputfile=sys.argv[1]
    outputdir_path=sys.argv[2]
    pid=sys.argv[3]
    print("Pyclone output file: ",pyclone_outputfile,"Output directory: ",outputdir_path,"Patient ID: ",pid)
    output1=outputdir_path+"/"+pid+"_lichee_sSNV.txt"
    output2=outputdir_path+"/"+pid+"_lichee_cluster.txt"

    #^^^Part one for LiCHEE input to generate the CP information # ******* CHANGE CP_cutoff if needed *******
    CP_cutoff = 0.03 #clusters below this cellular prevalence(CP)threshold will be set to 0, this is to remove low confidence read counts and hence CP for mutations in samples which appear because we consider a union of all SNVs between the METs

    f=open(pyclone_outputfile,"r")
    next(f)
    last_line1 = f.readlines()[-1]
    f.close()
    f=open(pyclone_outputfile,"r")
    next(f)
    df = pd.read_csv(pyclone_outputfile, sep="\t")
    samples = df['sample_id'].unique()
    header_list = list(samples)
    header = create_header(header_list) #Call the function to finish the preprocessing
    #final_header = "#chr\tposition\tdescription\tprofile\tNormal"+str(header)
    f_out=open(output1,"w")
    f_out.write(header)
    f_out.write("\n")
    #Function for SNV information generation
    licheeInputGenerationSNV(f,f_out)

    #Function for Cluster information generation
    ####licheeInputGenerationCluster(pyclone_outputfile,output2)
