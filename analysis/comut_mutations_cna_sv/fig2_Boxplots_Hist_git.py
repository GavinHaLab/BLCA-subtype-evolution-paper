
#/usr/bin/python3
#script used to plot Histology based plots for Fig.3 such as TMB, FGA etc and any box plots based on Histology
__Author__ = "Pushpa Itagi"
__Purpose__= "This script is used to generate plots for the Fig.2 of the manuscript"
#Importing required libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import sys
import matplotlib.patches as mpatches
import matplotlib as mpl
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
from scipy.stats import ttest_ind
# Set Type 42 font embedding
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# Add the histology information for patient IDs
UC=["19-001","18-065","19-044","15-108","16-079","17-020"]
PUC=["18-016","17-071","19-007","18-120","15-109"]
UCSD =["17-030","19-022","18-101","16-097","17-026","16-070"]
NE=["17-047","19-037"]
UC_Sarc = ["19-004"]


def tmbExonicPlot(in_file):
    print("Processing the file ", in_file)
    print("Make TMB plot per histology")
    # Load the data
    data_input = pd.read_csv(in_file)
    data_input['patient_id'] = data_input['sample'].str[:6] # uses the first 6 characters of the sample name as the patient id
    # Drop the hypermutated Patient as it skews the data
    data_input = data_input[~data_input['patient_id'].isin(['17-047'])]
    samples_to_drop = ["17-047pD10","18-120pA14","18-101M3","19-044pB11"] # these did not pass QC for mutation analysis
    data_input = data_input[~data_input['sample'].isin(samples_to_drop)]
    # drop samples containing a letter p which denotes the FFPEs
    data_input = data_input[~data_input['sample'].str.contains('p')]
    # Add histology based on the patient id
    global UC, PUC, UCSD, NE, UC_Sarc
    conditions = [ data_input["patient_id"].isin(UC),data_input["patient_id"].isin(PUC), data_input["patient_id"].isin(UCSD),data_input["patient_id"].isin(NE),data_input["patient_id"].isin(UC_Sarc)]
    choices = [ "UC","PUC", "UCSD", "NE", "UC_Sarc"]
    data_input["Histology"] = np.select(conditions, choices, default="Other")
    print(data_input["Histology"])
    plt.figure(figsize=(4, 4))
    histology_order = ["PUC","UC","UCSD","NE", "UC_Sarc"]
    data_input['Histology'] = pd.Categorical(data_input['Histology'], categories=histology_order, ordered=True)
    histology_colors = {'UC':'#815CA6','PUC':'#84C2EB','NE':'#E24D74','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D','UCSD':'#43823D','UCSD_focal':'#ACD7A0','UCSD_focal':'#ACD7A0','UC_Sarc':"#d2bea5"} # August 2024
    patients = {
        '15-109': '#e6194b', '16-070': '#3cb44b', '16-097': '#ffe119', '17-026': '#fabed4', '17-030': '#4363d8', '17-071': '#f58231', '18-101': '#42d4f4', '18-120': '#911eb4', '19-007': '#f032e6', '19-044': '#bfef45', '15-108': '#469990', '16-079': '#dcbeff', 
        '17-020': '#9a6324', '17-047': '#fffac8', '18-016': '#800000', '18-065': '#aaffc3', '19-001': '#808000', '19-004': '#ffd8b1', '19-022': '#000075','19-037': '#a9a9a9'}
    sns.boxplot(data=data_input,y="SNV_INDELS_coding_all_row_count",x="Histology",hue="Histology",notch=False,palette=histology_colors,width=0.9,showfliers=False,showcaps=False,medianprops={"color": "#000000", "linewidth": 3},gap=0.1,linecolor="#000000",whis=4.0,linewidth=1.7,whiskerprops={"linewidth": 1.5})
    sns.stripplot(data=data_input,y="SNV_INDELS_coding_all_row_count",x="Histology",hue="patient_id",palette=patients,dodge=False,size=8,jitter=0.25)
    ##sns.stripplot(data=data_input,y="SNV_INDELS_coding_all_row_count",x="Histology",color="#000000",dodge=False,size=7,jitter=0.2)
    ##plt.title('Tumor Mutational Burden')
    plt.ylabel('exonic region TMB',fontsize=14)
    plt.legend([],[], frameon=False)
    # remove side spines
    sns.despine()
    # set line width to 1.5
    ax = plt.gca()
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    # tick size to 18
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # Adjust y-axis limits to add spacing before the first box
    min_value = data_input["SNV_INDELS_coding_all_row_count"].min()
    max_value = data_input["SNV_INDELS_coding_all_row_count"].max()
    spacing = 0.5  # Adjust as needed
    plt.ylim(min_value - spacing, max_value + spacing)  # Add spacing
    plt.tight_layout()

    #Statistical test 
    # do a MWU test between UC and PUC and UC and UCSD and UCSD and PUC
    df_UC = data_input[data_input["Histology"] == "UC"]
    df_PUC = data_input[data_input["Histology"] == "PUC"]
    df_UCSD = data_input[data_input["Histology"] == "UCSD"]
    df_NE = data_input[data_input["Histology"] == "NE"]
    df_UC_Sarc = data_input[data_input["Histology"] == "UC_Sarc"]
    # Perform the Mann-Whitney U test
    stat_UC_PUC, p_UC_PUC = mannwhitneyu(df_UC['SNV_INDELS_coding_all_row_count'], df_PUC['SNV_INDELS_coding_all_row_count'])
    print('UC vs PUC MWU p-value=%.3f, p=%.4f' % (stat_UC_PUC, p_UC_PUC))
    stat_UC_UCSD, p_UC_UCSD = mannwhitneyu(df_UC['SNV_INDELS_coding_all_row_count'], df_UCSD['SNV_INDELS_coding_all_row_count'])
    print('UC vs UCSD MWU p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))
    stat_UCSD_PUC, p_UCSD_PUC = mannwhitneyu(df_UCSD['SNV_INDELS_coding_all_row_count'], df_PUC['SNV_INDELS_coding_all_row_count'])
    print('UCSD vs PUC MWU p-value=%.3f, p=%.4f' % (stat_UCSD_PUC, p_UCSD_PUC))
    # add the p-values to the plot
    plt.text(4, 9, 'UC-PUC p=%.3f' % p_UC_PUC, ha='center', va='center',fontsize=8)
    plt.text(4, 9.5, 'UC-UCSD p=%.3f' % p_UC_UCSD, ha='center', va='center',fontsize=8)
    plt.text(4, 10, 'UCSD-PUC p=%.3f' % p_UCSD_PUC, ha='center', va='center',fontsize=8)
    plt.savefig('./images/boxplot_SNV_INDEL_coding_excludes17047_TMB_by_Histology_Feb2025_colored_PUCFirst.pdf',bbox_inches=None)

    # LME model of PUC vs UC+UCSD+NE+SARC
    df_PUC['Group'] = 'PUC'
    df_combined = pd.concat([df_UC, df_UCSD, df_NE, df_UC_Sarc])
    df_combined['Group'] = 'UC_UCSD_NE_SARC'
    # Combine the two datasets
    df_lme = pd.concat([df_PUC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    lme_model = smf.mixedlm("SNV_INDELS_coding_all_row_count ~ Group", df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_UCSD_NE_SARC]'] 
    print('LME model  for PUC vs UC_UCSD_NE_SARC p-value=%.3f' % (p_value_lme))
      # LME model of PUC vs UC+UCSD
    df_PUC['Group'] = 'PUC'
    df_combined = pd.concat([df_UC, df_UCSD])
    df_combined['Group'] = 'UC_UCSD'
    # Combine the two datasets
    df_lme = pd.concat([df_PUC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    lme_model = smf.mixedlm("SNV_INDELS_coding_all_row_count ~ Group", df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_UCSD]'] 
    print('LME model  for PUC vs UC+UCSD p-value=%.3f' % (p_value_lme))
    return 0;

def averageMeanCPClonality(in_file):
    data_input = pd.read_csv(in_file)
    print("Inside the function to plot the average mean clonality")
    print(data_input)
    # rename column hist_patient to Histology
    data_input.rename(columns={"hist_patient":"Histology"},inplace=True)
    plt.figure(figsize=(4, 4))
    histology_order = ["UC", "UCSD", "PUC", "NE", "SARC"]
    data_input['Histology'] = pd.Categorical(data_input['Histology'], categories=histology_order, ordered=True)
    histology_colors = {'UC':'#815CA6','PUC':'#84C2EB','NE':'#E24D74','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D','UCSD':'#43823D','UCSD_focal':'#ACD7A0','UCSD_focal':'#ACD7A0','UC_Sarc':"#d2bea5","SARC":"#d2bea5"} # August 2024
    sns.boxplot(data=data_input,y="avg_CP",x="Histology",hue="Histology",notch=False,palette=histology_colors,width=0.9,showfliers=False,showcaps=False,medianprops={"color": "#000000", "linewidth": 3},gap=0.1,linecolor="#000000",whis=4.0,linewidth=1.7,whiskerprops={"linewidth": 1.5})
    sns.stripplot(data=data_input,y="avg_CP",x="Histology",color="#000000",dodge=False,size=7,jitter=0.2)
    ##plt.title('Tumor Mutational Burden')
    plt.ylabel('Mean Non-Founder CP',fontsize=14)
    plt.legend([],[], frameon=False)
    # remove side spines
    sns.despine()
    # set line width to 1.5
    ax = plt.gca()
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    # tick size to 18
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # Adjust y-axis limits to add spacing before the first box
    min_value = data_input["Histology"].min()
    max_value = data_input["Histology"].max()
    spacing = 0.1  # Adjust as needed
    ##plt.ylim(min_value - spacing, max_value + spacing)  # Add spacing
    plt.tight_layout()

    #Statistical test 
    # do a MWU test between UC and PUC and UC and UCSD and UCSD and PUC
    df_UC = data_input[data_input["Histology"] == "UC"]
    df_PUC = data_input[data_input["Histology"] == "PUC"]
    df_UCSD = data_input[data_input["Histology"] == "UCSD"]
    df_NE = data_input[data_input["Histology"] == "NE"]
    df_UC_Sarc = data_input[data_input["Histology"] == "UC_Sarc"]
    # Perform the Mann-Whitney U test
    stat_UC_PUC, p_UC_PUC = mannwhitneyu(df_UC['avg_CP'], df_PUC['avg_CP'])
    print('UC vs PUC MWU p-value=%.3f, p=%.4f' % (stat_UC_PUC, p_UC_PUC))
    stat_UC_UCSD, p_UC_UCSD = mannwhitneyu(df_UC['avg_CP'], df_UCSD['avg_CP'])
    print('UC vs UCSD MWU p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))
    stat_UCSD_PUC, p_UCSD_PUC = mannwhitneyu(df_UCSD['avg_CP'], df_PUC['avg_CP'])
    print('UCSD vs PUC MWU p-value=%.3f, p=%.4f' % (stat_UCSD_PUC, p_UCSD_PUC))
    # add the p-values to the plot
    plt.text(4, 0.9, 'UC-PUC p=%.3f' % p_UC_PUC, ha='center', va='center',fontsize=4)
    plt.text(4, 0.95, 'UC-UCSD p=%.3f' % p_UC_UCSD, ha='center', va='center',fontsize=4)
    plt.text(4, 1.0, 'UCSD-PUC p=%.3f' % p_UCSD_PUC, ha='center', va='center',fontsize=4)
    plt.savefig('./images/Mutations_meanCP_nonFounder_eachTumor_Feb2025.pdf',bbox_inches=None)
    return 0;


def fgaHistology(in_file):
    data_input = pd.read_csv(in_file)
    print("Inside the function to plot the FGA.....")
    print(data_input)
    # rename column hist_patient to Histology
    data_input.rename(columns={"hist_patient":"Histology"},inplace=True)
    plt.figure(figsize=(4, 4))
    histology_order = ["UCSD", "UC", "PUC", "NE", "SARC"]
    data_input['Histology'] = pd.Categorical(data_input['Histology'], categories=histology_order, ordered=True)
    histology_colors = {'UC':'#815CA6','PUC':'#84C2EB','NE':'#E24D74','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D','UCSD':'#43823D','UCSD_focal':'#ACD7A0','UCSD_focal':'#ACD7A0','UC_Sarc':"#d2bea5","SARC":"#d2bea5"} # August 2024
    sns.boxplot(data=data_input,y="FGA",x="Histology",hue="Histology",notch=False,palette=histology_colors,width=0.9,showfliers=False,showcaps=False,medianprops={"color": "#000000", "linewidth": 3},gap=0.1,linecolor="#000000",whis=4.0,linewidth=1.7,whiskerprops={"linewidth": 1.5})
    sns.stripplot(data=data_input,y="FGA",x="Histology",color="#000000",dodge=False,size=7,jitter=0.2)
    ##plt.title('Tumor Mutational Burden')
    plt.ylabel('Fraction genome Altered per Mb',fontsize=14)
 
    plt.legend([],[], frameon=False)
    # remove side spines
    sns.despine()
    # set line width to 1.5
    ax = plt.gca()
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    # tick size to 18
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # Adjust y-axis limits to add spacing before the first box
    min_value = data_input["Histology"].min()
    max_value = data_input["Histology"].max()
    spacing = 0.1  # Adjust as needed
    ##plt.ylim(min_value - spacing, max_value + spacing)  # Add spacing
    plt.tight_layout()

    #Statistical test 
    # do a MWU test between UC and PUC and UC and UCSD and UCSD and PUC
    df_UC = data_input[data_input["Histology"] == "UC"]
    df_PUC = data_input[data_input["Histology"] == "PUC"]
    df_UCSD = data_input[data_input["Histology"] == "UCSD"]
    df_NE = data_input[data_input["Histology"] == "NE"]
    df_UC_Sarc = data_input[data_input["Histology"] == "UC_Sarc"]
     # MWU is not performed rather a one-sided MWU is performed as we know the direction of the difference
    # do a one sided MWU between UCSD and UC+PUC+NE
    combined_FGA = pd.concat([df_UC['FGA'], df_PUC['FGA'], df_NE['FGA']])
    stat_UCSD_UC, p_UCSD_UC = mannwhitneyu(df_UCSD['FGA'], combined_FGA,alternative='greater')
    print('UCSD vs UC+PUC+NE 1 sided MWU p-value=%.3f, p=%.4f' % (stat_UCSD_UC, p_UCSD_UC))
    plt.text(2, -0.1, '1 sided MWU UCSD-UC+PUC+NE p=%.3f' % p_UCSD_UC, ha='center', va='center',fontsize=6)
    # two sided MWU between UC and UCSD+PUC+NE
    combined_FGA = pd.concat([df_UC['FGA'], df_PUC['FGA'], df_NE['FGA'], df_UC_Sarc['FGA']])
    stat_UC_UCSD, p_UC_UCSD = mannwhitneyu(df_UCSD['FGA'], combined_FGA)
    print('UC vs UCSD+PUC+NE MWU p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))
    # T-test for UCSD vs UC+PUC+NE
    combined_FGA = pd.concat([df_UC['FGA'], df_PUC['FGA']])
    stat_UCSD_UC, p_UCSD_UC = ttest_ind(df_UCSD['FGA'], combined_FGA)
    print('UCSD vs UC+PUC+NE  sided t-test p-value=%.3f, p=%.4f' % (stat_UCSD_UC, p_UCSD_UC))
    '''
    # one sided t-test between UC and UCSD
    from scipy.stats import ttest_ind
    stat_UC_UCSD, p_UC_UCSD = ttest_ind(df_UC['FGA'], df_UCSD['FGA'],alternative='greater')
    print('UC vs UCSD 1 sided t-test p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))
    plt.text(4, 1.1, 'UC-UCSD p=%.3f' % p_UC_UCSD, ha='center', va='center',fontsize=4)
    plt.savefig('./images/CN_FGA_byHistology_Feb2025.pdf',bbox_inches=None)
    '''
    return 0;


def cnNumbatMeanCP(in_file):
    data_input = pd.read_csv(in_file)
    print("Inside the function to plot the cnNumbatMeanCP .....")
    print(data_input)
    # rename column hist_patient to Histology
    data_input.rename(columns={"histology":"Histology"},inplace=True)
    print(data_input)
    plt.figure(figsize=(3, 3))
    histology_order = ["UC", "UCSD", "PUC"]
    data_input['Histology'] = pd.Categorical(data_input['Histology'], categories=histology_order, ordered=True)
    histology_colors = {'UC':'#815CA6','PUC':'#84C2EB','NE':'#E24D74','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D','UCSD':'#43823D','UCSD_focal':'#ACD7A0','UCSD_focal':'#ACD7A0','UC_Sarc':"#d2bea5","SARC":"#d2bea5"} # August 2024
    sns.boxplot(data=data_input,y="avg_CP",x="Histology",hue="Histology",notch=False,palette=histology_colors,width=0.9,gap=0.4,showfliers=False,showcaps=False,medianprops={"color": "#000000", "linewidth": 3},linecolor="#000000",whis=4.0,linewidth=1.7,whiskerprops={"linewidth": 1.5})
    sns.stripplot(data=data_input,y="avg_CP",x="Histology",color="#000000",dodge=False,size=9,jitter=0.25)
    ##plt.title('Tumor Mutational Burden')
    plt.ylabel('Mean CP copy number',fontsize=14)
    plt.legend([],[], frameon=False)
    # remove side spines
    sns.despine()
    # set line width to 1.5
    ax = plt.gca()
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    # tick size to 18
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # do a MWU test between UC and PUC and UC and UCSD and UCSD and PUC
    df_UC = data_input[data_input["Histology"] == "UC"]
    df_PUC = data_input[data_input["Histology"] == "PUC"]
    df_UCSD = data_input[data_input["Histology"] == "UCSD"]
    # Perform the Mann-Whitney U test UC vs UCSD, UC vs PUC, UCSD vs PUC
    stat_UC_UCSD, p_UC_UCSD = mannwhitneyu(df_UC['avg_CP'], df_UCSD['avg_CP'])
    print('UC vs UCSD MWU p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))
    stat_UC_PUC, p_UC_PUC = mannwhitneyu(df_UC['avg_CP'], df_PUC['avg_CP'])
    print('UC vs PUC MWU p-value=%.3f, p=%.4f' % (stat_UC_PUC, p_UC_PUC))
    stat_UCSD_PUC, p_UCSD_PUC = mannwhitneyu(df_UCSD['avg_CP'], df_PUC['avg_CP'])
    print('UCSD vs PUC MWU p-value=%.3f, p=%.4f' % (stat_UCSD_PUC, p_UCSD_PUC))
    # add the p-values to the plot
    #plt.text(4, 0.9, 'UC-UCSD p=%.3f' % p_UC_UCSD, ha='center', va='center',fontsize=4)
    #plt.text(4, 0.95, 'UC-PUC p=%.3f' % p_UC_PUC, ha='center', va='center',fontsize=4)
    #plt.text(4, 1.0, 'UCSD-PUC p=%.3f' % p_UCSD_PUC, ha='center', va='center',fontsize=4)
    plt.tight_layout()
    plt.savefig('./images/CN_Numbat_meanCP_perPatient.pdf',bbox_inches=None)
    return 0;


# This is to count the events for the paper based on the comut
def countEvents(in_file):
    print("Inside the function to count the events.....")
    df_in = pd.read_csv(in_file)
    # first three letters of sample are patient id
    df_in["patient_id"] = df_in["sample"].apply(lambda x: x[:4])
    # get all the unique values in raw_value column
    unique_values = df_in["raw_value"].unique()
    print("Unique values in raw_value column", unique_values)
    
    # If is_mutation is in raw_value then it is a mutation
    df_mutations = df_in[df_in["raw_value"].str.contains("is_mutation")]
    mutations_interest = ["Founder","Shared","Private"]
    # If mutations_interest is in value then retain those rows
    df_mutations = df_mutations[df_mutations["value"].isin(mutations_interest)]
    print("Mutations of interest", mutations_interest)
    # for each patient count the number of unique categories in the column value
    df_mutations_counts = df_mutations.groupby("patient_id")["category"].nunique().reset_index()
    print("Mutations counts", df_mutations_counts)
    # Calculate the sum of the last column
    sum_counts = df_mutations_counts.iloc[:, -1].sum()
    # Print the result
    print(f"Sum of counts in the last column for mutations: {sum_counts}")

      # If is_SV is in raw_value then it is a mutation
    df_sv = df_in[df_in["raw_value"].str.contains("is_SV")]
    sv_interest = ["Founder_SV","Shared_SV","Private_SV"]
    # If mutations_interest is in value then retain those rows
    df_sv = df_sv[df_sv["value"].isin(sv_interest)]
    print("SVs of interest", sv_interest)
    # drop the rows where the category is NO SV
    df_sv = df_sv[df_sv["category"] != "NO SV"]
    # for each patient count the number of unique categories in the column value
    sv_interest_counts = df_sv.groupby("patient_id")["category"].nunique().reset_index()
    ##sv_interest_counts = df_sv.groupby(["patient_id", "category"]).size().reset_index(name="unique_counts")

    print("SVs counts", sv_interest_counts)
    # Calculate the sum of the last column
    sum_counts = sv_interest_counts.iloc[:, -1].sum()
    # Print the result
    print(f"Sum of counts in the last column for SVs: {sum_counts}")

    # If its a CNA thencount the number of unique categories in the column value
    raw_values_cna = ["is_mutation","is_SV","NoCoverage","Indeterminate"]
    df_cna = df_in[~df_in["raw_value"].isin(raw_values_cna)]
    cna_interest = ["Founder_Amp","Shared_Amp","Private_Amp","Founder_Del","Shared_Del","Private_Del"]
    df_cna = df_cna[df_cna["value"].isin(cna_interest)]
    print("CNAs of interest", cna_interest)
    # for each patient count the number of unique categories in the column value
    df_cna_counts = df_cna.groupby("patient_id")["category"].nunique().reset_index()
    print("CNAs counts", df_cna_counts)
    # Calculate the sum of the last column
    sum_counts = df_cna_counts.iloc[:, -1].sum()
    # Print the result
    print(f"Sum of counts in the last column for CNAs: {sum_counts}")

    return 0;

# function to read per sample 1MB bin file and do a box plot histology
def fgaHistologyperSample(in_file):
    print("Inside the function to plot the FGA per sample.....")
    df_in = pd.read_csv(in_file,sep="\t")
    # add a patient_id which is the first six characters of the sample
    df_in["patient_id"] = df_in["Sample"].apply(lambda x: x[:6])
    print("Processing the file ", df_in)
    # assign histology based on the patient_id
    global UC, PUC, UCSD, NE, UC_Sarc
    conditions = [ df_in["patient_id"].isin(UC),df_in["patient_id"].isin(PUC), df_in["patient_id"].isin(UCSD),df_in["patient_id"].isin(NE),df_in["patient_id"].isin(UC_Sarc)]
    choices = [ "UC","PUC", "UCSD", "NE", "UC_Sarc"]
    df_in["Histology"] = np.select(conditions, choices, default="Other")
    print(df_in["Histology"])
    print("NUmber of columns in each row are", df_in.shape[1])
    df_in["length"] = df_in.apply(lambda row: len(row.dropna()), axis=1)
    # for each row/sample calcuate the number of columns that are do not contain the word Neutral
    df_in['non_neutral_count'] = df_in.apply(lambda row: sum('Neutral' not in str(value) for value in row), axis=1)
    # calculate  FGA by dividing the non_neutral_count by the length of the row
    df_in["FGA"] = df_in["non_neutral_count"]/df_in["length"]
    print("FGA", df_in["FGA"])
    data_input = df_in
    plt.figure(figsize=(4, 4))
    histology_order = ["UCSD", "UC", "PUC", "NE", "UC_Sarc"]
    data_input['Histology'] = pd.Categorical(data_input['Histology'], categories=histology_order, ordered=True)
    patients = {
    '15-109': '#e6194b', '16-070': '#3cb44b', '16-097': '#ffe119', '17-026': '#fabed4', '17-030': '#4363d8', '17-071': '#f58231', '18-101': '#42d4f4', '18-120': '#911eb4', '19-007': '#f032e6', '19-044': '#bfef45', '15-108': '#469990', '16-079': '#dcbeff', 
    '17-020': '#9a6324', '17-047': '#fffac8', '18-016': '#800000', '18-065': '#aaffc3', '19-001': '#808000', '19-004': '#ffd8b1', '19-022': '#000075','19-037': '#a9a9a9'}
    histology_colors = {'UC':'#815CA6','PUC':'#84C2EB','NE':'#E24D74','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D','UCSD':'#43823D','UCSD_focal':'#ACD7A0','UCSD_focal':'#ACD7A0','UC_Sarc':"#d2bea5","SARC":"#d2bea5"} # August 2024
    sns.boxplot(data=data_input,y="FGA",x="Histology",notch=False,palette=histology_colors,width=0.9,showfliers=False,showcaps=False,medianprops={"color": "#000000", "linewidth": 3},gap=0.1,linecolor="#000000",whis=4.0,linewidth=1.7,whiskerprops={"linewidth": 1.5})
    sns.stripplot(data=data_input,y="FGA",x="Histology",color="#000000",dodge=False,size=7,jitter=0.25,hue='patient_id', palette=patients)
    ##plt.title('Tumor Mutational Burden')
    plt.ylabel('Fraction genome Altered',fontsize=14)
 
    plt.legend([],[], frameon=False)
    # remove side spines
    sns.despine()
    # set line width to 1.5
    ax = plt.gca()
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    # tick size to 18
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # Adjust y-axis limits to add spacing before the first box
    min_value = data_input["Histology"].min()
    max_value = data_input["Histology"].max()
    spacing = 0.1  # Adjust as needed
    ##plt.ylim(min_value - spacing, max_value + spacing)  # Add spacing
    plt.tight_layout()
    plt.savefig('./images/CN_FGA_byHistology_perSample_Feb2025.pdf',bbox_inches=None)
    # do a MWU test on UCSD vs UC+PUC+NE
    df_UC = data_input[data_input["Histology"] == "UC"]
    df_PUC = data_input[data_input["Histology"] == "PUC"]
    df_UCSD = data_input[data_input["Histology"] == "UCSD"]
    df_NE = data_input[data_input["Histology"] == "NE"]
    df_SARC = data_input[data_input["Histology"] == "UC_Sarc"]
    df_combined = pd.concat([df_UC, df_PUC, df_NE, df_SARC])
    # two sided MWU between UCSD and UC+PUC+NE
    stat_UC_UCSD, p_UC_UCSD = mannwhitneyu(df_UCSD['FGA'], df_combined['FGA'])
    print('UCSD vs UC+PUC+NE_UC_SARC MWU p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))

    # LME model of UCSD vs UC+PUC+NE+UC_Sarc
    df_UCSD['Group'] = 'UCSD'
    df_combined['Group'] = 'UC_PUC_NE_SARC'
    # Combine the two datasets
    df_lme = pd.concat([df_UCSD, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    lme_model = smf.mixedlm("FGA ~ Group", df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_PUC_NE_SARC]'] 
    print('LME model  for UCSD vs UC+PUC+NE p-value=%.3f' % (p_value_lme))
    # write thd df_combined to a file
    df_combined.to_csv("../data/combined_FGA_UCSD_UC_PUC_NE.csv",index=False)
    #write the df_UCSD to a file
    df_UCSD.to_csv("../data/UCSD_FGA.csv",index=False)

     # LME model of UCSD vs UC+PUC
    df_combined = pd.concat([df_UC, df_PUC])
    df_UCSD['Group'] = 'UCSD'
    df_combined['Group'] = 'UC_PUC'
    # Combine the two datasets
    df_lme = pd.concat([df_UCSD, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    lme_model = smf.mixedlm("FGA ~ Group", df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_PUC]'] 
    print('LME model  for UCSD vs UC+PUC p-value=%.3f' % (p_value_lme))


    # LME model of PUC vs others
    df_combined = pd.concat([df_UC, df_UCSD, df_NE, df_SARC])
    # LME model of PUC vs UC+UCSD+NE+UC_Sarc
    df_PUC['Group'] = 'PUC'
    df_combined['Group'] = 'UC_UCSD_NE_SARC'
    # Combine the two datasets
    df_lme = pd.concat([df_PUC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    lme_model = smf.mixedlm("FGA ~ Group", df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_UCSD_NE_SARC]'] 
    print('LME model  for PUC vs UC+UCSD+PUC+NE p-value=%.3f' % (p_value_lme))

    # LME model of PUC vs UC+UCSD
    df_PUC['Group'] = 'PUC'
    df_combined['Group'] = 'UC_UCSD'
    # Combine the two datasets
    df_lme = pd.concat([df_PUC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    lme_model = smf.mixedlm("FGA ~ Group", df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_UCSD]'] 
    print('LME model  for PUC vs UC_UCSD p-value=%.3f' % (p_value_lme))
    return 0;



# Count the number of SV's per patient
def svPerTumor(in_file):
    print("Inside the function to plot the SV counts per tumor.....")
    data_input = pd.read_csv(in_file,sep="\t")
    print(data_input)
    ##data_input["patient_id"] = data_input["patient"]
    # Add histology based on the patient id
    global UC, PUC, UCSD, NE, UC_Sarc
    conditions = [ data_input["patient_id"].isin(UC),data_input["patient_id"].isin(PUC), data_input["patient_id"].isin(UCSD),data_input["patient_id"].isin(NE),data_input["patient_id"].isin(UC_Sarc)]
    choices = [ "UC","PUC", "UCSD", "NE", "UC_Sarc"]
    data_input["Histology"] = np.select(conditions, choices, default="Other")
    print(data_input["Histology"])
    plt.figure(figsize=(4, 4))
    ##histology_order = ["UC", "UCSD", "PUC", "NE", "UC_Sarc"]
    histology_order = ["UC", "UCSD", "PUC", "NE", "UC_Sarc"]
    patients = {
        '15-109': '#e6194b', '16-070': '#3cb44b', '16-097': '#ffe119', '17-026': '#fabed4', '17-030': '#4363d8', '17-071': '#f58231', '18-101': '#42d4f4', '18-120': '#911eb4', '19-007': '#f032e6', '19-044': '#bfef45', '15-108': '#469990', '16-079': '#dcbeff', 
        '17-020': '#9a6324', '17-047': '#fffac8', '18-016': '#800000', '18-065': '#aaffc3', '19-001': '#808000', '19-004': '#ffd8b1', '19-022': '#000075','19-037': '#a9a9a9'}
    data_input['Histology'] = pd.Categorical(data_input['Histology'], categories=histology_order, ordered=True)
    histology_colors = {'UC':'#815CA6','PUC':'#84C2EB','NE':'#E24D74','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D','UCSD':'#43823D','UCSD_focal':'#ACD7A0','UCSD_focal':'#ACD7A0','UC_Sarc':"#d2bea5","SARC":"#d2bea5"} # August 2024
    ##column_of_interest = "counts"
    column_of_interest = "protein_coding_genes"
    sns.boxplot(data=data_input,y=column_of_interest,x="Histology",notch=False,palette=histology_colors,width=0.9,showfliers=False,showcaps=False,medianprops={"color": "#000000", "linewidth": 3},gap=0.1,linecolor="#000000",whis=4.0,linewidth=1.7,whiskerprops={"linewidth": 1.5})
    sns.stripplot(data=data_input,y=column_of_interest,x="Histology",hue="patient_id",palette=patients,dodge=False,size=7,jitter=0.2)
    ##plt.title('Tumor Mutational Burden')
    plt.ylabel('SV counts per tumor',fontsize=14)
    plt.legend([],[], frameon=False)
    # remove side spines
    sns.despine()
    # set line width to 1.5
    ax = plt.gca()
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    # tick size to 18
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # Adjust y-axis limits to add spacing before the first box
    min_value = data_input["Histology"].min()
    max_value = data_input["Histology"].max()
    spacing = 0.1  # Adjust as needed
    ##plt.ylim(min_value - spacing, max_value + spacing)  # Add spacing
    plt.tight_layout()
    out_file = "../images/SV_counts_consensus_TwoOverlap_protein_coding_genes.pdf" #~~~~~~~ Change the output file name as needed
    plt.savefig(out_file,bbox_inches=None)
    # do a MWU test on UCSD vs UC+PUC+NE
    df_UC = data_input[data_input["Histology"] == "UC"]
    df_PUC = data_input[data_input["Histology"] == "PUC"]
    df_UCSD = data_input[data_input["Histology"] == "UCSD"]
    df_NE = data_input[data_input["Histology"] == "NE"]
    df_SARC = data_input[data_input["Histology"] == "UC_Sarc"]
    df_combined = pd.concat([df_UCSD, df_PUC, df_NE, df_SARC])
    # two sided MWU between UCSD and UC+PUC+NE
    stat_UC_UCSD, p_UC_UCSD = mannwhitneyu(df_UC[column_of_interest], df_combined[column_of_interest])
    print('UC vs UCSD+PUC+NE_UC_SARC MWU p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))
    # stat test between UC and PUC
    stat_UC_PUC, p_UC_PUC = mannwhitneyu(df_UC[column_of_interest], df_PUC[column_of_interest])
    print('UC vs PUC MWU p-value=%.3f, p=%.4f' % (stat_UC_PUC, p_UC_PUC))
    # stat test between UC and UCSD
    stat_UC_UCSD, p_UC_UCSD = mannwhitneyu(df_UC[column_of_interest], df_UCSD[column_of_interest])
    print('UC vs UCSD MWU p-value=%.3f, p=%.4f' % (stat_UC_UCSD, p_UC_UCSD))
    df_combined = pd.concat([df_UCSD, df_UC, df_NE, df_SARC])
    # stat test between PUC and others
    stat_PUC_UCSD, p_PUC_UCSD = mannwhitneyu(df_PUC[column_of_interest], df_combined[column_of_interest])
    print('PUC vs UCSD+UC+NE_UC_SARC MWU p-value=%.3f, p=%.4f' % (stat_PUC_UCSD, p_PUC_UCSD))
    # LME model of PUC vs UC+PUC+NE+UC_Sarc
    df_PUC['Group'] = 'PUC'
    df_combined = pd.concat([df_UC, df_UCSD, df_NE, df_SARC])
    df_combined['Group'] = 'UC_UCSD_NE_SARC'
    # Combine the two datasets
    df_lme = pd.concat([df_PUC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    ####lme_model = smf.mixedlm("column_of_interest ~ Group", df_lme, groups=df_lme["patient_id"])
    formula = f"{column_of_interest} ~ Group"
    lme_model = smf.mixedlm(formula, df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_UCSD_NE_SARC]'] 
    print('LME model  for PUC vs UC+UCSD+PUC+NE p-value=%.3f' % (p_value_lme))

     # LME model of PUC vs UC+UCSD
    df_PUC['Group'] = 'PUC'
  
    df_combined = pd.concat([df_UCSD, df_PUC])
    df_combined['Group'] = 'UC_UCSD'
    # Combine the two datasets
    df_lme = pd.concat([df_PUC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    ####lme_model = smf.mixedlm("column_of_interest ~ Group", df_lme, groups=df_lme["patient_id"])
    formula = f"{column_of_interest} ~ Group"
    lme_model = smf.mixedlm(formula, df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UC_UCSD]'] 
    print('LME model  for PUC vs UC_UCSDp-value=%.3f' % (p_value_lme))

     # LME model of UC vs PUC+UCSD
    df_UC['Group'] = 'UC'
    df_combined = pd.concat([df_UCSD, df_PUC])
    df_combined['Group'] = 'UCSD_PUC'
    print(df_combined)
    print(df_UC)
    # Combine the two datasets
    df_lme = pd.concat([df_UC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    ####lme_model = smf.mixedlm("column_of_interest ~ Group", df_lme, groups=df_lme["patient_id"])
    formula = f"{column_of_interest} ~ Group"
    lme_model = smf.mixedlm(formula, df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    print(lme_result.pvalues.filter(like='Group'))
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.UCSD_PUC]'] 
    print('LME model  for UC vs UCSD_PUC p-value=%.3f' % (p_value_lme))

    # LME model of UC vs UC+PUC+NE+UC_Sarc
    df_UC['Group'] = 'UC'
    df_combined = pd.concat([df_UCSD, df_PUC, df_NE, df_SARC])
    df_combined['Group'] = 'PUC_UCSD_NE_SARC'
    # Combine the two datasets
    df_lme = pd.concat([df_PUC, df_combined])
    print(df_lme.head())
    # Fit the LME model
    # Assuming 'FGA' is the dependent variable and 'Group' is the fixed effect
    ####lme_model = smf.mixedlm("column_of_interest ~ Group", df_lme, groups=df_lme["patient_id"])
    formula = f"{column_of_interest} ~ Group"
    lme_model = smf.mixedlm(formula, df_lme, groups=df_lme["patient_id"])
    lme_result = lme_model.fit()
    # Print the summary of the LME model
    print(lme_result.summary())
    # P-value from the LME model
    p_value_lme = lme_result.pvalues['Group[T.PUC_UCSD_NE_SARC]'] 
    print('LME model  for UC vs PUC_UCSD_NE_SARC p-value=%.3f' % (p_value_lme))
    return 0;


# main function
if __name__ == "__main__":
    print("This script is used to generate plots for the Fig.2 of the manuscript")
    tmbExonicPlot("./data/SNV_coding_excludeSynonymous_INDEL_coding_all_Final.csv")
    averageMeanCPClonality("./data/Mutations_meanCP_nonFounder_eachTumor.csv")
    fgaHistology("./data/CN_FGA_byHistology.csv")
    cnNumbatMeanCP("./data/CN_Numbat_meanCP_perPatient.csv")

    #SV analysis
    # per tumor SV counts
    svPerTumor("./data/SV_counts_per_patient_consensus_twoOverlap.txt") # Using consensus SV numbers
