
#!/usr/bin/python
# ===============================================================
# Script       : comut_timing_withCNA_git.py
# Purpose      : Script to generate a comut/oncoprint plot using mutations, CNAS and SVs. Make a landscape styled Figure 3A comut plot
# Usage        : Run on terminal as python3 comut_timing_withCNA_git.py
# Dependencies : comput, pandas
# Author       : Pushpa Itagi
# Note:
#   - This script uses a input parameters that can be modified:
#       * input_file directories    – path to input data where we have annotation data such as ./../data/mutations_cna_sv_in_allDrivergenes_noMixed_filtered_singleOccurrence.csv, 
#       * input_file directories    – snv_coding_nonsynonmous_splicing.csv, WESWGS_variant.csv, WESWGS_tissue.csv, ../../data/snv_coding_nonsynonmous_splicing.csv, ../../data/tumorFraction_ploidy.csv
#       * output_file directories   – path to output data where files are stored

# Load required modules
import sys
sys.path.insert(0, '/path/To/Comut/package') # Adjust the path to where comut is installed *******
import comut
#from comut import comut
import fileparsers
import palettable
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patches as patches
import os
sys.path.insert(0, os.getcwd()) #Refer to the file name with sample IDs which is the analysis folder this is to use different Sample IDs like a shorter version if needed.
from sampleNameRename_UpsetComut import *
import matplotlib as mpl

# Set Type 42 font embedding
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
driver_genes_blca = "./drivergene_manual_annotations_finalBLCAPaperJune2025.csv"


# Function to generate side bar data for the comut plot
def generateSideBarData(mutation_data):
    print(" Inside side bar generation")
    # get all the unique values in the category column
    categories = mutation_data['category'].unique()
    values = mutation_data['value'].unique()
    # remove no mutation from the values
    values = values[(values != 'No mutation') & (values != 'No power')]
     #  add a column for each value in the value column
    columns =  list(values)
    print("columns",columns)
    side_bar_data = pd.DataFrame(columns=columns)
    for category in categories:
        # get the unique values in the value column
        for value in values:
            count = mutation_data[(mutation_data['category'] == category) & (mutation_data['value'] == value)].shape[0]
            side_bar_data.loc[category, value] = count
    print("side bar data",side_bar_data.head(2))
    # convert the  index to a column
    side_bar_data = side_bar_data.reset_index().rename(columns={'index': 'category'})
    return side_bar_data;

# Function to calculate proportions of mutation types
def calculate_proportions(mutation_data):
    # Filter relevant mutation types
    mutation_data_filtered = mutation_data[mutation_data['value'].isin(['Founder', 'Shared', 'Private', 'Founder_Amp', 'Founder_Del', 'Shared_Amp', 'Shared_Del', 'Private_Amp', 'Private_Del', 'Founder_SV', 'Shared_SV', 'Private_SV'])]
    
    # Group by sample and calculate value counts
    mutation_data_filtered_bar = mutation_data_filtered.groupby('sample')['value'].value_counts().unstack().fillna(0).reset_index()
    # Calculate total mutations
    mutation_data_filtered_bar["Total"] = mutation_data_filtered_bar.sum(axis=1)
    # Calculate proportions
    for category in ['Founder', 'Shared', 'Private']:
        mutation_data_filtered_bar[f"{category}_sum"] = mutation_data_filtered_bar.filter(like=category).sum(axis=1)
        mutation_data_filtered_bar[f"{category}_proportion"] = mutation_data_filtered_bar[f"{category}_sum"] / mutation_data_filtered_bar["Total"]
    # Retain only sample and proportions
    result = mutation_data_filtered_bar[['sample', 'Founder_proportion', 'Shared_proportion', 'Private_proportion']]
    return result;

# Function to rename the sample names in the mutation data
def sampleNameRename_UpsetComut(mutation_data,matched_column,new_column,sample_column):
    print("Inside the function to rename the samples",mutation_data.head(2),".....",sample_column)
    new_df = pd.DataFrame()
    # Iterate over all samples in mutation_data
    for sample in mutation_data[sample_column]:
        sample_name = sample.split(".txt")[0]
        if sample_name in all_samples:
            mutation_data['sample'] = mutation_data[sample_column].replace(sample, all_samples[sample_name])
    # rename second column to value
    mutation_data = mutation_data.rename(columns={matched_column: 'value'})
    
    samples_to_drop =['P004Lym2', 'P120Pri', 'P047Pri', 'P101Pri', 'P044Pri']
    mutation_data = mutation_data[~mutation_data[sample_column].isin(samples_to_drop)]
    # convert the values to log scale only for the TMB
    if new_column == 'log_tmb':
        mutation_data['value'] = np.log(mutation_data['value'])
        mutation_data["category"] = "TMB"
    else:
        mutation_data["category"] = "TFx"
    
    # rename value to log_tmb
    mutation_data = mutation_data.rename(columns={'value': new_column})
    # make sure sample is the first column followed bycategory and value and then the rest of the columns
    mutation_data = mutation_data[['sample',new_column,'category']]
    return mutation_data;

custom_rcParams = {'font.family': 'Times New Roman','font.size': 8}
# update rcParams
rcParams.update(custom_rcParams)
plt.rcParams["font.family"] = "Arial" # if you don't have arial it will fall back to Deja Vu
font = {'family': 'serif','color':  'darkred', 'weight': 'normal','size': 8}
fig = plt.figure(figsize=(25, 5), tight_layout=True)
ax = fig.add_subplot(111)

object_comut = comut.CoMut()
figsize = (20,26) # the size of the figure - changes the shape of the squares in the comut
dpi = 100 # change the output resolution
extension = '.pdf' # extension for saving - can otherwise be .pdf, .png, .jpg, etc
bar_mapping = {'Met Counts': '#275056'}
gender_colors={'male':'#759383','female':'#53465c'}
mapping_colors_binary = {'Yes':'#5034a4','No':'#c6ccc2'}
smoking_colors = {'Never':'#c6ccc2','past':'#737770','current':'#372a24'}
value_order = ['Never', 'past', 'current']
indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1.5, 'markersize': 3.0}
cat_mapping = {'Absent': {'facecolor': '#000000'}}
side_mapping = {'counts': 'darkgrey'}
bar_kwargs = {'width': 0.7}

# set the order of the values in the mutation data
value_order_mutations = ['Founder', 'Shared', 'Private','No mutation','No power','Below Threshold','Founder_SV','Shared_SV','Private_SV','Founder_Amp','Founder_Del','Shared_Amp','Shared_Del','Private_Amp','Private_Del','NoCoverage','No sCNV','Indeterminate']
priority = ['Founder', 'Shared', 'Private','Founder_Amp','Founder_Del','Shared_Amp','Shared_Del','Private_Amp','Private_Del','Founder_SV','Shared_SV','Private_SV']

# import categorical data and continuous data (purity)
indicator_data = pd.read_csv('../../data/sample_indicators_histology.csv', sep = ',') # prev: histology_set2.csv
mutation_data = pd.read_csv("../../data/mutations_cna_sv_in_allDrivergenes_noMixed_filtered_singleOccurrence.csv",sep=",")
print("Number of unique Genes in this comut ",mutation_data['category'].nunique())
print("\n\n\nThe unique values are ",mutation_data['value'].unique())

# Drop genes where all samples have no mutation
####mutation_data = mutation_data.groupby('category').filter(lambda x: (x['value'] != 'No mutation').any())
mutation_data = mutation_data.groupby('category').filter(lambda x: (x['value'] != 'No mutation').any() and (x['value'] != 'No sCNV').any())
variant_data = pd.read_csv('../../data/WESWGS_variant.csv', sep = ',')#
#category_order = sorted(mutation_data['category'].unique())
tissue_data = pd.read_csv('../../data/WESWGS_tissue.csv', sep = ',')
side_bar_data = generateSideBarData(mutation_data)
print("returned side bar data",side_bar_data.head(2))
tmb_coding_data = pd.read_csv("../../data/snv_coding_nonsynonmous_splicing.csv", sep = ',')
tmb_coding_data = sampleNameRename_UpsetComut(tmb_coding_data,'SNV_coding_non_synonymous_row_count','log_tmb','sample')
tumor_Fraction_data = pd.read_csv("../../data/tumorFraction_ploidy.csv", sep = ',')
# rename sample_id to sample
tumor_Fraction_data = tumor_Fraction_data.rename(columns={'sample_id': 'sample'})
tumor_Fraction_data = sampleNameRename_UpsetComut(tumor_Fraction_data,'tumor_purity','tfx','sample')

#drop the genes from the side bar data that are not in the mutation data
####side_bar_data = side_bar_data[side_bar_data['category'].isin(mutation_data['category'])]
# Create a new column 'total' that is the sum of 'founder', 'shared', and 'private' columns
side_bar_data['total'] = side_bar_data['Founder'] + side_bar_data['Shared'] + side_bar_data['Private']
print("side bar data",side_bar_data.head(2))
# order the columns as follows 
#side_bar_data = side_bar_data[['category','Founder','Shared','Subclonal','Private','total']]
side_bar_data = side_bar_data[['category','Founder','Shared','Private','total']]
# Sort the DataFrame by the 'total' column in descending order
side_bar_data = side_bar_data.sort_values(by='total', ascending=False)
#drop the 'total' column
side_bar_data = side_bar_data.drop(columns=['total'])
category_order = (side_bar_data['category'].unique())
# reverse the order of the genes
category_order = category_order[::-1]
###print("Category order is ......",category_order)

#Order the genes in the mutation data based on groups Tumor suppressor genes, Oncogenes, and other genes from the driver genes file
driver_genes = pd.read_csv(driver_genes_blca, sep = ',')
tsg=driver_genes[driver_genes['Category']=='Tumor Suppressor']
onco=driver_genes[driver_genes['Category']=='Oncogene']
exclude_categories = ['Tumor Suppressor', 'Oncogene','None']
other = driver_genes[~driver_genes['Category'].isin(exclude_categories)]
# write to a file
other.to_csv('../../results/other.csv', index=False)
##other = driver_genes[driver_genes['Update Made']=='Others'] # Include the genes that are not in the TSG and Oncogene category but literature curated and have some role in cancer

# Convert all the gene names to list
tsg = tsg['gene'].tolist();onco = onco['gene'].tolist();other = other['gene'].tolist()
# Sort the genes based on tumor suppressor genes, oncogenes and other genes
new_category_order = other +  onco + tsg 
# Remove the genes that are not in the mutation data
new_category_order = [gene for gene in new_category_order if gene in category_order]
###print("New category order is ......",new_category_order)

# Drop any genes that have Only No mutation and No sCNV
mutation_data = mutation_data.groupby('category').filter(
    lambda x: not all(value in ['No mutation', 'No sCNV'] for value in x['value']))
print("Filtered genes in TSG/Oncogene that have matched cna direction",mutation_data['category'].nunique())
values_to_count = ["Founder", "Shared", "Private", "Founder_Amp", "Shared_Amp", "Private_Amp", "Founder_Del", "Shared_Del", "Private_Del"]


# Filter mutation_data to include only rows where 'category' is in tsg and 'value' is in values_to_count
filtered_tsg_data = mutation_data[mutation_data['category'].isin(tsg) & mutation_data['value'].isin(values_to_count)]
# Group by 'category' and count the occurrences
gene_counts = filtered_tsg_data['category'].value_counts()
# Sort the genes based on the count in descending order
sorted_genes_tsg = gene_counts.index.tolist();sorted_genes_tsg = sorted_genes_tsg[::-1]
###print("TSG genes sorted based on the count",sorted_genes_tsg)
# Same for oncogenes
filtered_onco_data = mutation_data[mutation_data['category'].isin(onco) & mutation_data['value'].isin(values_to_count)]
gene_counts = filtered_onco_data['category'].value_counts()
sorted_genes_onco = gene_counts.index.tolist();sorted_genes_onco = sorted_genes_onco[::-1]
#reverse the order of the genes
###print("Oncogenes sorted based on the count",sorted_genes_onco)
# Same for other genes
filtered_other_data = mutation_data[mutation_data['category'].isin(other) & mutation_data['value'].isin(values_to_count)]
gene_counts = filtered_other_data['category'].value_counts()
sorted_genes_other = gene_counts.index.tolist();sorted_genes_other = sorted_genes_other[::-1]
###print("Other genes sorted based on the count",sorted_genes_other)
grouped_sorted_order = sorted_genes_other + sorted_genes_onco + sorted_genes_tsg #Included ONCOGENES, TSG and other genes
##grouped_sorted_order = sorted_genes_onco + sorted_genes_tsg #Included ONCOGENES and TSG only
###print("Grouped sorted order",grouped_sorted_order)
print("Final number of genes in the comut plot",len(grouped_sorted_order))
gene_colors = {"tsg": "blue", "onco": "red", "other": "black"}

# extract the first 4 characters of the sample name to get the patient ID
indicator_data['Analysis_ID' ] = indicator_data['sample'].str.slice(0, 4)
print(indicator_data.head(2))
patients_list = indicator_data['Analysis_ID'].unique()
print(patients_list)

# Drop the samples that did not pass the QC
samples_to_drop = ["P047Pri","P120Pri","P044Pri","P101Pri","P004Lym3","P004Lym2"]
mutation_data = mutation_data[~mutation_data['sample'].isin(samples_to_drop)]
indicator_data = indicator_data[~indicator_data['sample'].isin(samples_to_drop)]
variant_data = variant_data[~variant_data['sample'].isin(samples_to_drop)]
tissue_data = tissue_data[~tissue_data['sample'].isin(samples_to_drop)]

priority = ['Missense']
# Define the color mapping for the sources
#histology_colors = {'UC':'#815CA6','PUC':'#F9CADF','NE':'#84C2EB','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D'}
histology_colors = {'UC':'#815CA6','PUC':'#84C2EB','NE':'#E24D74','UC_SQFocal':'#ACD7A0','UC_SQ':'#43823D','SQ':'#43823D','UCSD':'#43823D','UCSD_focal':'#ACD7A0','UCSD_focal':'#ACD7A0','UC_Sarc':"#d2bea5"} # August 2024
####"Shared-Primary-Metastasis":"#95bbc6",
mapping_colors_mutations = {"Founder": "#FEC216", "Shared": "#2B3990", "Private": "#C91220", 'No mutation': {'facecolor': '#e9e9e9', 'alpha': 1.0},"No power": 'white',"Subclonal": "#2B3990","Below Threshold":"#a4d5d0","NoCoverage":"black",'No sCNV': {'facecolor': '#e9e9e9', 'alpha': 1.0},
                            "Indeterminate":"pink","Founder_Amp": "#FEC216", "Shared_Amp": "#2B3990", "Private_Amp": "#C91220","Founder_Del": "#FEC216", "Shared_Del": "#2B3990", "Private_Del": "#C91220",
                            "Founder_SV":{'facecolor': 'none', 'edgecolor':'#FEC216', 'linewidth': 3.0,'hatch':''}, "Shared_SV":{'facecolor':'none', 'edgecolor':'#2B3990', 'linewidth': 3.0,'hatch':''}, "Private_SV":{'facecolor':'none', 'edgecolor':'#C91220', 'linewidth': 3.0,'hatch':''},
                            "Founder_proportion": "#FEC216", "Shared_proportion": "#2B3990", "Private_proportion": "#C91220",}
tissue_colors_mapping = {'Bladder':'black','Primary':'black','Primary':'#750029','Liver':'#bc65b5','Lung':'#8d7f95','Lymph':'#f8aea7','Other':'#c1d6ef'} #'Bladder':'#750029'
side_bar_mapping = {'Num. of mutations': 'darkgrey'}
borders = ["Founder_SV", "Shared_SV", "Private_SV"]
TMB_mapping ={"value":"grey"}
# setting the heights of the plots
##heights = {'Same patient': 1, 'Histology': 1.0, 'Tissue': 1.0, 'Mutation burden': 1.75,'BLCA Del':1.75,'TMB':1.75,'TFx':1.75}
heights = {'Same patient': 1, 'Histology': 1.0, 'Tissue': 1.0, 'Mutation burden': 1.75,'BLCA Del':1.75}
mapping_clonal_evolution = {'Diminished':'#a1a59e','Continued':'black'}
# Add the data to the comut object
object_comut.add_sample_indicators(indicator_data, name = 'Same patient', plot_kwargs = indicator_kwargs)
#object_comut.add_continuous_data(ichor_data, name = 'Tfx',mapping='Reds',cat_mapping = cat_mapping, value_range = value_range)
mutation_data.to_csv('../../results/mutation_data_FINAL_comut.csv', index=False)
object_comut.add_categorical_data(mutation_data, name = 'Mutations',mapping=mapping_colors_mutations,category_order = grouped_sorted_order,value_order=value_order_mutations,borders=borders)
##object_comut.add_side_bar_data(side_bar_data, paired_name="Mutations", name = 'Frequency',mapping = side_bar_mapping, xlabel = 'Mutsig -log(Q)', position = 'right', bar_kwargs = bar_kwargs) # add the side bar data
object_comut.add_categorical_data(variant_data, name = 'Histology', mapping = histology_colors)
object_comut.add_categorical_data(tissue_data, name = 'Tissue', mapping = tissue_colors_mapping)
#object_comut.add_categorical_data(clonal_evolution_data, name = 'Clonal Divergence', mapping = mapping_clonal_evolution)
bar_kwargs_mutations = {'width': 0.8}
# Get a frequencey of the number of alterations founder, shared and private in each sample from mutation_data
proportions_df = calculate_proportions(mutation_data)
object_comut.add_bar_data(proportions_df, name = "Prop. Alterations", mapping = mapping_colors_mutations, stacked=True,bar_kwargs = bar_kwargs_mutations) # if using the color mapping
object_comut.add_bar_data(tmb_coding_data, name = 'TMB',ylabel = 'TMB in Mb(log)',mapping =  {'log_tmb': '#435b57'},bar_kwargs = bar_kwargs_mutations)
print("Tumor fraction data",tumor_Fraction_data.head(2))
print("TMB data",tmb_coding_data.head(2))
object_comut.add_bar_data(tumor_Fraction_data, name = 'TFx',ylabel = 'Tumor Fraction',mapping =  {'tfx': '#435b57'},bar_kwargs = bar_kwargs_mutations)


# Add the comut plot
comut = object_comut.plot_comut(figsize = (18, 30), x_padding = 0.00, y_padding = 0.00, tri_padding = 0.00, hspace = 0.01, heights=heights,wspace=0.04)
ax.set_xlabel('xlabel', fontsize=5)
ax.set_ylabel('ylabel', fontsize=3)
font = {'family': 'serif', 'color':  'darkred', 'weight': 'normal', 'size': 5}
ax.set_yticklabels(ax.get_yticklabels(), fontdict=font)
rcParams.update(custom_rcParams)
object_comut.axes['Same patient'].tick_params(axis = 'x', rotation=90,labelsize=15)
object_comut.axes['Mutations'].tick_params(axis='y', labelsize=22) # for gene label size
def color_gene_labels(ax, genes, color):
    ax.figure.canvas.draw()  # Ensure labels are rendered before modifying
    for label in ax.get_yticklabels():
        gene_name = label.get_text().strip()  # Strip spaces for matching
        ##print(f"Processing: {gene_name}")  # Debugging
        if gene_name in genes:
            label.set_color(color)
            label.set_fontstyle('normal')
    ax.figure.canvas.draw()  # Force redraw to apply changes
# color the gene labels based on the driver genes
color_gene_labels(object_comut.axes['Mutations'], onco, color='#8c1400')  # Tumor suppressor genes
color_gene_labels(object_comut.axes['Mutations'], tsg, color='#005ac8') # Oncogenes
color_gene_labels(object_comut.axes['Mutations'], other, color='black')   # Other genes

# Add a vertical line at each boundary on each axis
#Get the current x-axis labels
labels = [item.get_text() for item in object_comut.axes['Same patient'].get_xticklabels()]
# Identify the unique patient IDs and their positions
unique_patients, indices = np.unique([label[:4] for label in labels], return_index=True)
#increase the indices by another 1
indices = np.append(indices, len(labels))

# Add a rectangle for each patient on the figure
fig = object_comut.figure
center_list=[]
y = 0.15 # for the y position
#y = 0.5 # for the y position
# Add a rectangle for each patient on the figure which is four samples wide
count = 0
for i in range(0, len(indices)-3, 1):  # Adjust the range to len(indices)-3
     # Transform the x and y coordinates from the axes data coordinate system to the figure coordinate system
    trans_start_data = object_comut.axes['Same patient'].transData.transform([(indices[i], 0)])
    trans_end_data = object_comut.axes['Same patient'].transData.transform([(indices[i+3], 1)])
    trans_start = fig.transFigure.inverted().transform(trans_start_data)
    trans_end = fig.transFigure.inverted().transform(trans_end_data)
    # Draw the rectangle on the figure
    rect = patches.Rectangle((trans_start[0, 0], 0.120), width=(trans_end[0, 0]-trans_start[0, 0]), height=0.765, transform=fig.transFigure, color='#000000', linewidth=1.0, fill=False) # the value in trans_start adjusts the start of rectangle i.e. comut bottom and the height specifies the rectangle
    fig.patches.extend([rect]) #color='#7e7c7c'
    x = rect.get_x()
    print("rectangle",x,patients_list[count+1],rect)
    # Add text on top of the rectangle
    ##fig.text(trans_start[0, 0] + (trans_end[0, 0] - trans_start[0, 0]) / 2, 0.9, patients_list[count+1], ha='center', va='bottom', transform=fig.transFigure)
    count += 1

# Get current y-axis labels (gene names)
y_labels = [label.get_text() for label in ax.get_yticklabels()]
# Assign colors based on gene category
label_colors = []
for gene in y_labels:
    if gene in tsg:
        label_colors.append(gene_colors["tumor_suppressor"])
    elif gene in onco:
        label_colors.append(gene_colors["oncogene"])
    else:
        label_colors.append(gene_colors["other"])

# Apply colored y-axis labels
for label, color in zip(ax.get_yticklabels(), label_colors):
    label.set_color(color)
    
# Add space after each patient
# change the font of genes on the y-axis
font = {'family': 'serif', 'color':  'darkred', 'weight': 'normal', 'size': 5}
ax.set_yticklabels(ax.get_yticklabels(), fontdict=font)
ax.tick_params(axis='y', labelsize="12")
for ax in object_comut.axes.values():
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
# tight_layout() doesn't take into account the title, so we need to adjust the top

# Save the figure
object_comut.figure.savefig('../../results/Fig2_comut_withMutationsSVCNA' + extension, dpi=300,bbox_inches = 'tight')
#End of the script
