# ===============================================================
# Script       : SVClone_merge_git.R
# Purpose      : Script to generate a final SV file for each patient by iterating over per sample runs of the tumor and then merging the results.
# Dependencies : tools, ggforce, ggforce and dplyr libraries
# Author       : Developed by Alan Min and modified by Pushpa Itagi
# Note: This file has a dependency, please run the SVCloneFix_git.R script first to generate the required files.
#   - This script uses a input parameters that can be modified:
#       * input_file directories    – path to input data where we have per sample SV clone runs. Make sure the SVCloneFix_git.R hasa been run before as this script needs a file which is generated from the dependecny script. 
#       * output_file directories   – path to output data where files are stored
#       * Important variables that can be modified: threshold this determines the distance between two SVs to be considered the same. Useful for combining SVs that are the same between samples in a patient.

#Load required libraries
library(tools)
#library(ggforce)
#library(ggforce)
library(dplyr)


# Checks if start2 is within threshold of start1 and if end2 is within threshold of end1
within_threshold = function(start1, end1, start2, end2, threshold) {
  if (abs(start1 - start2) <= threshold & abs(end1 - end2) <= threshold) {
    return(T)
  } else {
    return(F)
  }
}

# This function takes a list of samples and a threshold, and returns a list of SVs per sample, a master list of SVs, and the updated sample list with SV names assigned. 
# ***threshold is important as it determines how close two SVs need to be to be considered the same.
call_svs_per_sample = function(sample_list, threshold) {
  # Master list of names for all identified SVs
  master_sv = data.frame(name=NULL, start=NULL, end=NULL)
  
  # Counter for giving unique names to each identified SV
  sv_counter = 1
  
  # Big list to keep track of what SVs were detected in which samples
  all_svs_per_sample = list()
  
  # Run through each tissue sample or cell-free sample
  for (list_name in names(sample_list)) {
    sample_sv = sample_list[[list_name]]
    # Running list of SVs detected in this sample
    sample_sv_names = c()
    # Add a name column to the sample_sv
    sample_list[[list_name]]$SV = "unassigned"
    # Run through all the SVs in sample
    for (i in 1:nrow(sample_sv)) {
      start1 = sample_sv$start[i]
      end1 = sample_sv$end[i]
      chrom_start1 = sample_sv$chrom_start[i]
      chrom_end1 = sample_sv$chrom_end[i]
      found = F
      if (nrow(master_sv) >= 1) {
        for (j in 1:nrow(master_sv)) {
          start2 = master_sv$start[j]
          end2 = master_sv$end[j]
          chrom_start2 = master_sv$chrom_start[j]
          chrom_end2 = master_sv$chrom_end[j]
          # If the sample matches something in the master list, take the name from the master list
          # Also double check that the chromosomes match
          if (within_threshold(start1, end1, start2, end2, threshold) &&
              chrom_start1 == chrom_start2 &&
              chrom_end1 == chrom_end2) {
            sv_name = master_sv$name[j]
            sample_list[[list_name]]$SV[i] = sv_name
            found = T
          }
        }
      }
      # If the sample is not found in the master list, make a new name
      if (!found) {
        sv_name = sprintf("SV%04d", sv_counter)
        sample_list[[list_name]]$SV[i] = sv_name
        sv_counter = sv_counter + 1
      }
      sample_sv_names = c(sample_sv_names, sv_name)
      master_sv = rbind(master_sv, data.frame(start=start1, end=end1, chrom_start = chrom_start1, chrom_end = chrom_end1, name=sv_name, sample=list_name))
    }
    all_svs_per_sample[[list_name]] = sample_sv_names
  }
  ret = list(all_svs_per_sample = all_svs_per_sample, master_sv = master_sv, sample_list = sample_list)
  return(ret)
}

# Finds all the samples in the folder with sample_name and reads them into a sample_list format
make_sample_list = function(folder, sample_name) {
  current_samples = list.files(folder)[grepl(sample_name, list.files(folder))]
  
  sample_list = list()
  print(sample_list)
  for (current_sample in current_samples) {
    data_list = read.table(paste0(folder, current_sample, '/ccube_out/', current_sample, '_cluster_certainty.txt'), header=T)
    more_info = read.table(paste0(folder, current_sample, '/', current_sample, '_filtered_svs.tsv'), header=T)
    more_info = more_info[, 8:ncol(more_info)]
    names(data_list)[1:5] = c('chrom_start', 'start', 'dir1', 'chrom_end', 'end')
    data_list = cbind(data_list, more_info)
    cluster_ccf = read.table(paste0(folder, '/', current_sample, '/ccube_out/', current_sample, '_subclonal_structure_with_CCF_summary.txt'), header=T) # Read this file which is genrated from the SVCloneFix_git.R script
    data_list$CCF = cluster_ccf$CCF[data_list$most_likely_assignment]
    sample_list[[current_sample]] = data_list
  }
  return(sample_list)
}

# This function makes a merged table of SVs across all samples in the sample_list.
make_merged_table = function(sample_list, threshold = 500) {
  ret = call_svs_per_sample(sample_list, threshold=threshold)
  max_sv = max(as.integer(substring(ret$master_sv$name, 3)))
  df = data.frame(SV = sprintf("SV%04d", 1:max_sv))
  for (sample_name in names(ret$sample_list)) {
    sample = ret$sample_list[[sample_name]]
    CCF = rep(NA, max_sv)
    clonality = rep(NA, max_sv)
    sv_ids = as.integer(substring(sample$SV, 3))
    for (i in 1:nrow(sample)) {
      sv_id = sv_ids[i]
      CCF[sv_id] = sample$CCF[i]
      if (sample$CCF[i] > 0.9) {
        clonality[sv_id] = "Clonal"
      } else {
        clonality[sv_id] = "Subclonal"
      }
    }
    df[paste0(sample_name, "CCF")] = CCF
    df[paste0(sample_name, "clonality")] = clonality
  }
  patient_clonality_func = function(i) {
    x = df[i, ]
    n_sample = (length(x) - 1) / 2
    if (sum(x == "Clonal", na.rm=T) == n_sample && sum(x == "Subclonal", na.rm=T) == 0) {
      return("Founder Clone")
    } else if (sum(x == "Clonal", na.rm=T) > 1 && sum(x == "Subclonal", na.rm=T) == 0) {
      return("Shared Clonal/Subclonal")
    } else if (sum(x == "Clonal", na.rm=T) == 1 && sum(x == "Subclonal", na.rm=T) == 0) {
      return("Private Clonal")
    } else if (sum(x == "Clonal", na.rm=T) > 0 && sum(x == "Subclonal", na.rm=T) > 0) {
      return("Shared Clonal/Subclonal")
    } else if (sum(x == "Subclonal", na.rm=T) > 1 && sum(x == "Clonal", na.rm=T) == 0) {
      return("Shared Subclone")
    } else if (sum(x == "Subclonal", na.rm=T) == 1 && sum(x == "Clonal", na.rm=T) == 0) {
      return("Private Subclone")
    } else {
      stop("THIS SHOULDNT HAPPEN")
    }
    return(sum(x == "AHHH THIS IS AN ERROR", na.rm=T) + sum(x == "Subclonal", na.rm=T))
  }
  df$patient_clonality = sapply(1:nrow(df), patient_clonality_func)
  
  all_samples = bind_rows(ret$sample_list, .id = "column_label")
  classification <- all_samples %>% 
    group_by(SV) %>%
    summarise(
      chrom_start = first(chrom_start),  # or use unique(chr_start) if it can vary
      chrom_end = first(chrom_end),
      #start2 = min(start2),
      #end2 = max(end2),# same as above
      min_start = min(min(start), min(end)),
      max_end = max(max(start), max(end)),
      classification_list = paste(classification, collapse = "|")
    )
  df = merge(df, classification)
  return(df)
}

folder = './TwoOverlap_SVClone_output' #very important this folder needs to have a per sample run form SVClone

# Get the unique sample names from the folder
sample_names = unique(substr(list.files(folder), 1, 6))
for (sample_name in sample_names) {
  print(sample_name)
  print("...")
  sample_list = make_sample_list(folder, sample_name)
  merged_table = make_merged_table(sample_list, threshold = 500)
  write.csv(merged_table, paste0("../../results/sv_clone/", sample_name, ".csv"), row.names = F)
  
}


# Make a plot comparing subclonal vs subclonal SVs per patient
output_folder = "../../results/sv_clone/"
output_files = list.files(output_folder, full.names = T)
big_plot_df = data.frame()
for (file in output_files) {
  # file = output_files[1]
  dat = read.csv(file,check.names = F)
  clonality = dat[, grepl("X.*clonality", names(dat))]
  clonal = rowSums(as.matrix(clonality) == "Clonal", na.rm = T)
  subclonal = rowSums(as.matrix(clonality) == "Subclonal", na.rm = T)
  n_sample = ncol(clonality)
  plot_df = data.frame(clonal, subclonal, n_sample, sample=file_path_sans_ext(basename(file)))
  radius = 0.4
  noise1 = runif(min = -radius, max = radius, n=nrow(plot_df))
  noise2 = runif(min = -sqrt(radius^2 - noise1^2), max = sqrt(radius^2 - noise1^2), n=nrow(plot_df))
  plot_df$clonal = plot_df$clonal + noise1
  plot_df$subclonal = plot_df$subclonal + noise2
  
  big_plot_df = rbind(big_plot_df, plot_df)
}

n_sample_df = big_plot_df %>% group_by(sample) %>% summarise(n_sample = unique(n_sample))
circle_df = data.frame()
for (i in 1:nrow(n_sample_df)) {
  n = n_sample_df$n_sample[i]
  sample = n_sample_df$sample[i]
  circles = expand.grid(0:n, 0:n)
  circles = circles %>% filter((Var1 + Var2) <= n)
  circles$sample = sample
  
  circle_df = rbind(circle_df, circles)
}
ggplot(big_plot_df) + geom_point(aes(x=subclonal, y=clonal), size = 0.5) + geom_abline(aes(intercept=n_sample, slope=-1), data=n_sample_df) + geom_circle(aes(x0=Var1, y0=Var2, r=0.42), data=circle_df, color="#69bcf0") + facet_wrap(vars(sample))
ggsave("../../results/sv_clone/subclone_vs_clone.pdf", width=9, height=8)
write.table(big_plot_df,file="../../results/sv_clone/subclone_vs_clone.txt",row.names = F,col.names = T,quote = F,sep="\t")

# Make a plot comparing types of SVs to clonality
output_folder = "../../results/sv_clone/"
output_files = list.files(output_folder, full.names = T)
big_plot_df = data.frame()
for (file in output_files) {
  dat = read.csv(file)
  plot_dat = data.frame(patient_clonality = dat$patient_clonality, classification=dat$classification_list, patient=file_path_sans_ext(basename(file)), sv_size = dat$max_end - dat$min_start)
  big_plot_df = rbind(big_plot_df, plot_dat)
}

silly_unique = function(x) {
  svtypes = unlist(strsplit(x, "\\|"))
  if (!(sum(svtypes == svtypes[1]) == length(svtypes))) {
    # print(svtypes)
    warning("There is an example of where not all the entries of classification_list are the same")
    return("Mixed")
  }
  return(svtypes[1])
}
big_plot_df$classification = sapply(big_plot_df$classification, silly_unique)

width = 5
height = 4

# SV Type by Clonality Type
tbl = table(big_plot_df$patient_clonality, big_plot_df$classification)
write.table(tbl,file="../../results/sv_clone/type_of_sv_proportion_value_value2.txt",row.names = F,col.names = T,quote = F,sep="\t")

ggplot(as.data.frame(tbl)) + geom_col(aes(x=Var1, y=Freq, fill=Var2), position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + ylab("Proportion") + guides(fill=guide_legend(title="Type of SV"))
ggsave("../../results/sv_clone/type_of_sv_proportion.pdf", width=width, height=height)

ggplot(as.data.frame(tbl)) + geom_col(aes(x=Var1, y=Freq, fill=Var2), position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + guides(fill=guide_legend(title="Type of SV"))
ggsave("../../results/sv_clone/type_of_sv_value2.pdf", width=width, height=height)

ggplot(as.data.frame(tbl)) + geom_col(aes(x=Var1, y=Freq, fill=Var2)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + ylab("Proportion") + guides(fill=guide_legend(title="Type of SV"))
ggsave("../../results/sv_clone/type_of_sv_value.pdf", width=width, height=height)

clonal_generalize_convert = function(x) {
  if (x == 'Private Clonal' || x == 'Founder Clone') {
    return('Clonal')
  } else if (x == 'Private Subclone' || x == 'Shared Subclone') {
    return('Subclonal')
  } else if (x == "Shared Clonal/Subclonal"){
    return('Shared Clonal/Subclonal')
  } else {
    stop("Clonality is one of these options")
  }
}
big_plot_df$general_clonality = sapply(big_plot_df$patient_clonality, clonal_generalize_convert)
write.table(big_plot_df,file="../../sv_clone/results/type_of_sv.txt",row.names = F,col.names = T,quote = F,sep="\t")

no_mixed_big_plot_df = big_plot_df %>% filter(general_clonality != 'Mixed')
tmp = as.data.frame(prop.table(table(no_mixed_big_plot_df$general_clonality, no_mixed_big_plot_df$classification), margin = 1))
write.table(tmp,file="../../results/sv_clone/generalized_clonality_proportion.txt",row.names = F,col.names = T,quote = F,sep="\t")

ggplot(tmp) + geom_col(aes(x=Var1, y=Freq, fill=Var2), position='dodge') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + ylab('Proportion') + guides(fill=guide_legend(title="Type of SV"))
ggsave("../../results/sv_clone/generalized_clonality_proportion.pdf", width=width, height=height)

# Chi Square Test to see if there is a difference in DEL and INV in clonal vs subclonal
t = table(no_mixed_big_plot_df$general_clonality, no_mixed_big_plot_df$classification)
mat = matrix(0, nrow=2, ncol=2)
mat[1,1] = t['Clonal', 'DEL']
mat[1,2] = sum(t['Clonal', ]) - t['Clonal', 'DEL']
mat[2,1] = t['Subclonal', 'DEL']
mat[2,2] = sum(t['Subclonal', ]) - t['Subclonal', 'DEL']
chisq.test(mat)

mat = matrix(0, nrow=2, ncol=2)
mat[1,1] = t['Clonal', 'INV']
mat[1,2] = sum(t['Clonal', ]) - t['Clonal', 'INV']
mat[2,1] = t['Subclonal', 'INV']
mat[2,2] = sum(t['Subclonal', ]) - t['Subclonal', 'INV']
chisq.test(mat)

ggplot(big_plot_df) + geom_bar(aes(x=patient_clonality, fill=classification), position="fill") + facet_wrap(vars(patient)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + ylab("Proportion") + guides(fill=guide_legend(title="Type of SV"))
ggsave("../../results/sv_clone/type_of_sv_proportion_by_patient.pdf", width=width, height=height)

ggplot(big_plot_df) + geom_bar(aes(x=patient_clonality, fill=classification)) + facet_wrap(vars(patient)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + guides(fill=guide_legend(title="Type of SV"))
ggsave("../../results/sv_clone/type_of_sv_value_by_patient.pdf", width=width, height=1.5*height)

ggplot(big_plot_df) + geom_bar(aes(x=patient_clonality, fill=classification), position='dodge') + facet_wrap(vars(patient)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + guides(fill=guide_legend(title="Type of SV"))
ggsave("../../results/sv_clone/type_of_sv_value_by_patient.pdf", width=width, height=1.5*height)

ggplot(big_plot_df) + geom_bar(aes(x=patient_clonality)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality Type")
ggsave("../../results/sv_clone/clonality.pdf", width=width, height=height)

ggplot(big_plot_df) + geom_bar(aes(x=classification)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("SV Type")
ggsave("../../results/sv_clone/sv_type.pdf", width=width, height=height)

ggplot(big_plot_df) + geom_bar(aes(x=patient, fill=classification), position="fill")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../../results/sv_clone/classification_by_patient.pdf", width=width, height=height)
table(with(big_plot_df, prop.table(table(patient, classification), margin=1)))

ggplot(big_plot_df) + geom_bar(aes(x=patient, fill=patient_clonality), position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../../results/sv_clone/clonality_by_patient.pdf", width=width, height=height)

ggplot(big_plot_df) + geom_histogram(aes(sv_size), color="white") + xlab("Size of SV")
ggsave("../../results/sv_clone/sv_size_histogram.pdf", width=width, height=height)

ggplot(big_plot_df,aes(y=sv_size, x=patient_clonality)) + geom_boxplot() + xlab("Size of SV") + scale_y_continuous(trans='log10') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Size of SV (log scale)") +
  stat_compare_means(method = "t.test",size=3, label = "p.format", comparisons = combn(unique(big_plot_df$patient_clonality), 2, simplify = FALSE)) + geom_jitter(color="black", size=0.4, alpha=0.9)
ggsave("../../results/sv_clone/sv_size_histogram_by_patient_clonality_t.test.pdf", width=width, height=height)

## Wilcoxon Rank-Sum Test
pairwise_results <- pairwise.wilcox.test(big_plot_df$sv_size, big_plot_df$patient_clonality, p.adjust.method = "bonferroni")
comparisons1 <- combn(unique(big_plot_df$patient_clonality), 2, simplify = FALSE)
ggplot(big_plot_df,aes(y=sv_size, x=patient_clonality)) + geom_boxplot(alpha = 0.6) + geom_jitter(color="black", size=0.4, alpha=0.9) + stat_compare_means(method = "wilcox.test", comparisons = comparisons1,size=3) +
  theme_minimal() + labs(title = "Pairwise Wilcoxon Rank-Sum Test") + theme(legend.position = "none") + ylab("Size of SV (log scale)") 
ggsave("../../results/sv_clone/sv_size_histogram_by_patient_clonality_wilcox.test.pdf", width=width, height=height)

#ggplot(big_plot_df, aes(y=sv_size, x=classification)) + geom_boxplot() + xlab("Size of SV") + scale_y_continuous(trans='log10') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Size of SV (log scale)") + geom_signif(comparisons = list(c("DEL", "DUP")), map_signif_level = T)
ggplot(big_plot_df, aes(y=sv_size, x=classification)) + geom_boxplot() + xlab("Size of SV") + scale_y_continuous(trans='log10') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Size of SV (log scale)") + 
  stat_compare_means(method = "t.test",size=3, label = "p.format", comparisons = combn(unique(big_plot_df$classification), 2, simplify = FALSE)) + geom_jitter(color="black", size=0.4, alpha=0.9)
ggsave("../../results/sv_clone/sv_size_histogram_by_sv_type.pdf", width=width, height=height)

## Wilcoxon Rank-Sum Test
pairwise_results <- pairwise.wilcox.test(big_plot_df$sv_size, big_plot_df$classification, p.adjust.method = "bonferroni")
comparisons2 <- combn(unique(big_plot_df$classification), 2, simplify = FALSE)
ggplot(big_plot_df,aes(y=sv_size, x=classification)) + geom_boxplot(alpha = 0.6) + geom_jitter(color="black", size=0.4, alpha=0.9) + stat_compare_means(method = "wilcox.test", comparisons = comparisons2,size=3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(title = "Pairwise Wilcoxon Rank-Sum Test") + theme(legend.position = "none") + scale_y_continuous(trans='log10') + ylab("Size of SV (log scale)") + xlab("Size of SV")
ggsave("../../results/sv_size_histogram_by_sv_wilcox.test.pdf", width=width, height=height)

# Function to plot the proportion of generalized clonality by size cutoff
size_cutoff = function(cutoff) {
  cutoff_no_mixed_big_plot_df = no_mixed_big_plot_df %>% filter(sv_size > cutoff)
  tmp = as.data.frame(prop.table(table(cutoff_no_mixed_big_plot_df$general_clonality, cutoff_no_mixed_big_plot_df$classification), margin = 1))
  ggplot(tmp) + geom_col(aes(x=Var1, y=Freq, fill=Var2), position='dodge') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Clonality") + ylab('Proportion') + guides(fill=guide_legend(title="Type of SV"))
}
size_cutoff(10000)

patient_details = read.csv("patient_BLCA_TAN_details.csv")
names(patient_details)[1] = "patient"

big_plot_df = merge(big_plot_df, patient_details, by = "patient")
ggplot(big_plot_df) + geom_bar(aes(x=histology, fill=classification), position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Histology")

ggplot(big_plot_df) + geom_bar(aes(x=histology, fill=patient_clonality), position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Histology")

ggplot(big_plot_df) + geom_bar(aes(x=histology, fill=general_clonality), position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Histology")

ggplot(big_plot_df) + geom_violin(aes(x=histology, y=sv_size), position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Histology")

ggplot(big_plot_df) + geom_bar(aes(x=cisplatin, fill=classification), position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Cisplatin")

ggplot(big_plot_df) + geom_bar(aes(x=cisplatin, fill=patient_clonality), position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Cisplatin")

ggplot(big_plot_df) + geom_bar(aes(x=cisplatin, fill=general_clonality), position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Cisplatin")

ggplot(big_plot_df) + geom_violin(aes(x=cisplatin, y=sv_size), position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Cisplatin")

