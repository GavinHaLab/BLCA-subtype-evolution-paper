
# ===============================================================
# Script       : complexmap_final_git.R
# Purpose      : Script to generate a heatmap for the TFBS nucleosome coverage data
# Usage        : Run in command line RStudio or  with `Rscript complexmap_final_git.R` 
# Dependencies : ComplexHeatmap, circlize, RColorBrewer,viridis
# Author       : Pushpa Itagi
# ===============================================================
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
# Sample data: 11 genes x 8 samples
set.seed(1)
mat <- as.matrix(read.csv("../data/heatmap_final_list_paper_zscore_git.csv",row.names=1))
print(dim(mat)) 
rownames(mat) <- paste0("Gene", 1:9) #~~~~~~~ Change this to match your gene names
colnames(mat) <- paste0("Sample", 1:8) #~~~~~~~ Change this to match your sample names

rownames(mat) <- c(
  "ASCL1", "SOX2", "HES1", "RARG", "GRHL2",
  "PPARG", "ZNF217", "ZNF146","HMGA1")
colnames(mat) <- c(
  "15-109", "19-037", "19-001", "16-070","19-022", "16-079", 
  "16-097", "18-101")
# Annotation data
tumor_location <- c( "LuminalU","NE", "LuminalP", "LuminalP","Mixed", "Basal_Squamous", "Basal_Squamous","Basal_Squamous")
cancer_type <- c("PUC",  "NE", "UC", "UCSD", "UCSD","UC","UCSD","UCSD")
tumor_fraction <- c( 0.5, 0.85,0.65, 0.20, 0.08,0.44,0.26,0.27)

# Top annotation
ha <- HeatmapAnnotation(
  Histology = cancer_type,
  Molecular_Subtype = tumor_location,
  Fraction = anno_barplot(tumor_fraction, gp = gpar(fill = "grey")),
  col = list(
    Molecular_Subtype = c("NE" = "#911EB4", "LuminalP" = "#FFE119", "Basal_Squamous" = "#E6194B","LuminalU" = "#4363D8","Mixed" = "#92ee02","Healthy" = "#A9A9A9"),
    Histology = c("UC" = "#815CA6", "PUC" = "#84C2EB", "UCSD" = "#43823D","NE"="#E24D74","Healthy"="#A9A9A9")
  ),
  annotation_name_side = "left"
)

# Heatmap color function

col_fun <- colorRamp2(c(-1.8, 0.1, 1.8), c("red", "black", "black"))

# Define gaps after gene 4 and gene 8 (between 4+4+3 groups)
row_gaps <- unit.c(unit(2, "mm"), unit(2, "mm"), unit(2, "mm"), unit(2, "mm"))
# print the number of rows in the matrix
print(paste("Number of rows in the matrix:", nrow(mat)))

# open a pdf to save
pdf("heatmap_final_list_paper_zscore_git.pdf", width = 9, height = 7)
# drop the Sample 7 column from the mat

# Draw heatmap
Heatmap(
  mat,
  name = "Nucleosome coverage",
  top_annotation = ha,
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_split = factor(c(rep("Group1", 2), rep("Group2", 1), rep("Group3", 3), rep("Group4", 2), rep("Group5", 1)), levels = c("Group1", "Group2", "Group3","Group4","Group5")),
  gap = row_gaps,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Samples",
  border = TRUE,
  row_title = "Transciption factors",
)
message("finished drawing heatmap")
dev.off()
