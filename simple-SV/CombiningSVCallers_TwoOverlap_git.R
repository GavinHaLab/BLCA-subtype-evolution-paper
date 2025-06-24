
# ===============================================================
# Script       : CombiningSVCallers_TwoOverlap.R
# Purpose      : Script to generate a consenus SV call from three SV callers (MANTA, GRIDDS, and SVABA) by finding overlaps between them
# Usage        : Run in command line RStudio or  with `Rscript CombiningSVCallers_TwoOverlap.R` 
# Dependencies : GenomicRanges, stringr, dplyr,plyr,GenomeInfoDb,data.table
# Author       : Manasvita Vashisth modified by Pushpa Itagi
# Note:
#   - This script uses a few hard-coded parameters that can be modified:
#       * input_file directories    – path to input data where MANTA, GRIDDS, and SVABA VCF files are stored
#       * output_file directories   – path to output data where results will be saved 
#       * sample_names file         – path to a text file containing sample names to be processed
#       * maxgap1 = 1000 # maximum gap size for finding overlaps this is the size between two breakpoints to be considered as an overlap

#   - These are defined near the top of the script under “User Settings
# ===============================================================

# Clear all variables from the workspace
rm(list=ls())

# Load necessary libraries
library(data.table)
library(GenomicRanges)
library(VariantAnnotation)
library(stringr)
library(plyr)
library(vcfR)
library(dplyr)


# Set working directory
  # [USER-DEFINED] Input file paths
setwd("/path/to/thisScriptDirectory/")
sample_names = readLines("samples.txt") # a file with sample names, one per line and each sample should have gridds, manta, and svaba VCF files in the specified directories

# Loop over all sample names
for(i in 1:length(sample_names))
{
  # Read VCF files for each sample
  # [USER-DEFINED] Input file path
  gridds=VariantAnnotation::readVcf(paste('/path/to/gridds/',sample_names[i],'.gripss.filtered.vcf',sep=''),genome = 'hg38')
  manta=VariantAnnotation::readVcf(paste('/path/to/manta/',sample_names[i],'/results/variants/somaticSV.vcf',sep=''),genome = 'hg38')
  svaba=VariantAnnotation::readVcf(paste('/path/to/svaba/',sample_names[i],'/',sample_names[i],'.svaba.somatic.sv.vcf',sep=''),genome = 'hg38')
  # remove any lines that contain chrUn
  ##gridds=gridds[!grepl("chrUn", gridds$seqnames),]
  ##manta=manta[!grepl("chrUn", manta$seqnames),]
  ##svaba=svaba[!grepl("chrUn", svaba$seqnames),]
  #manta=manta[abs(info(manta)$SVLEN)>1000]
  # Convert VCF to GenomicRanges
  gridds_gr <- breakpointRanges(gridds)
  gridds_gr=gridds_gr[seqnames(gridds_gr) %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'),]
  gridds_gr=gridds_gr[hasPartner(gridds_gr)]
  gridds_gr <- partner(gridds_gr)
  manta_gr <- breakpointRanges(manta)
  manta_gr=manta_gr[seqnames(manta_gr) %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'),]
  manta_gr=subset(manta_gr,abs(manta_gr$svLen)>1000) # TAKE ONLY SVs WITH LENGTH > 1000bp
  manta_gr=manta_gr[hasPartner(manta_gr)]
  manta_gr <- partner(manta_gr,selfPartnerSingleBreakends = FALSE)
  svaba_gr <- breakpointRanges(svaba)
  svaba_gr=svaba_gr[seqnames(svaba_gr) %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'),]
  svaba_gr=svaba_gr[hasPartner(svaba_gr)]
  svaba_gr <- partner(svaba_gr)
  
  maxgap1=1000 # [EDITABLE] CNV segment size threshold (bp)
  # Find overlaps between different SV callers
  overlap1=findBreakpointOverlaps(svaba_gr,gridds_gr,maxgap = maxgap1)
  overlap_ms=findBreakpointOverlaps(manta_gr,svaba_gr,maxgap = maxgap1)
  overlap_gm=findBreakpointOverlaps(gridds_gr,manta_gr,maxgap = maxgap1)
  
  gridds_overlap=svaba_gr[queryHits(overlap1)]
  a=as.data.frame(unique(gridds_overlap))
  ms=as.data.frame(unique(svaba_gr[subjectHits(overlap_ms)]))
  gm=as.data.frame(unique(gridds_gr[queryHits(overlap_gm)]))
  #overlap2=findBreakpointOverlaps(gridds_overlap,manta_gr,maxgap = 2000)
  #intersection_gr=gridds_overlap[queryHits(overlap2)]
  #intersection1=as.data.frame(intersection_gr)
  #merged_gr <- c(svaba_gr,gridds_gr)
  #merged_gr1=reduce(c(merged_gr,manta_gr),maxgap=2000) #union2 GR
  #sm_sg=findBreakpointOverlaps(svaba_gr[subjectHits(overlap_ms)],gridds_overlap,maxgap=maxgap1)
  #s_overlap=svaba_gr[queryHits(sm_sg)]
  #s_
  twooverlap=c(svaba_gr[subjectHits(overlap_ms)],gridds_overlap)
  twooverlap=unique(c(twooverlap,gridds_gr[queryHits(overlap_gm)]))
  twooverlap=twooverlap[hasPartner(twooverlap)]
  twooverlap1=as(twooverlap,'data.frame')
  bed=breakpointgr2bedpe(twooverlap)
  bed$start1=bed$start1+1 #adjusting for the bedpe format where chromosome start position is 0
  bed$start2=bed$start2+1
  bed$end1=bed$end1+1
  bed$end2=bed$end2+1
  bed$SPAN=ifelse(bed$chrom1==bed$chrom2,abs(bed$start1-bed$start2),-1)
  #bed=bed[bed$chrom1 %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'),]
  #bed=bed[bed$chrom2 %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'),]
  
  #bed$chromosome_1=str_extract(bed$chrom1, "(?<=chr)\\w+")
  #bed$chromosome_2=str_extract(bed$chrom2, "(?<=chr)\\w+")
  bed$chromosome_1=bed$chrom1
  bed$chromosome_2=bed$chrom2
  bed$orient_1=ifelse(bed$strand1=='-','fwd','rev')
  bed$orient_2=ifelse(bed$strand2=='-','fwd','rev')
  bed$DR=ifelse(bed$name %in% rownames(geno(svaba)$DR),geno(svaba)$DR[match(bed$name,rownames(geno(svaba)$DR)),ncol(geno(svaba)$DR)],NA) #discordant reads
  bed$SR=ifelse(bed$name %in% rownames(geno(svaba)$SR),geno(svaba)$SR[match(bed$name,rownames(geno(svaba)$SR)),ncol(geno(svaba)$SR)],NA) #split reads
  #bed$AD=ifelse(bed$name %in% rownames(geno(svaba)$DP),geno(svaba)$AD[match(bed$name,rownames(geno(svaba)$DP)),ncol(geno(svaba)$DP)],NA) #read depth
  bed$DP=ifelse(bed$name %in% rownames(geno(svaba)$DP),geno(svaba)$DP[match(bed$name,rownames(geno(svaba)$DP)),ncol(geno(svaba)$DP)],NA) #read depth
  bed$GT=ifelse(bed$name %in% rownames(geno(svaba)$GT),geno(svaba)$GT[match(bed$name,rownames(geno(svaba)$GT)),ncol(geno(svaba)$GT)],NA) #genotype
  #bed$SCTG=ifelse(bed$name %in% rownames(svaba),info(svaba)$SCTG[match(bed$name,rownames((svaba)))],NA)
  bed$DR=ifelse(bed$name %in% rownames(geno(gridds)$RP),geno(gridds)$RP[match(bed$name,rownames(geno(gridds)$RP)),ncol(geno(gridds)$RP)],bed$DR)
  bed$SR=ifelse(bed$name %in% rownames(geno(gridds)$SR),geno(gridds)$SR[match(bed$name,rownames(geno(gridds)$SR)),ncol(geno(gridds)$SR)],bed$SR)
  bed$DP=ifelse(bed$name %in% rownames(info(gridds)),info(gridds)$VF[match(bed$name,rownames(info(gridds)))]/info(gridds)$TAF[match(bed$name,rownames(info(gridds)))],bed$DP)
  bed$GT=ifelse(bed$name %in% rownames(geno(gridds)$GT),geno(gridds)$GT[match(bed$name,rownames(geno(gridds)$GT)),ncol(geno(gridds)$GT)],bed$GT)
  bed$vaf=ifelse(bed$SR!=0,bed$SR/bed$DP,bed$DR/bed$DP)
  bed$FILTER=rep('PASS',nrow(bed))
  bed$alt_1=ifelse(bed$name %in% twooverlap1$sourceId,twooverlap1$ALT[match(bed$name,twooverlap1$sourceId)],NA)
  bed$alt_2=ifelse(bed$name %in% twooverlap1$partner ,twooverlap1$ALT[match(bed$name,twooverlap1$partner)],NA)
  bed$mateID=ifelse(bed$name %in% twooverlap1$sourceId,twooverlap1$partner[match(bed$name,twooverlap1$sourceId)],NA)
  bed$REF=ifelse(bed$name %in% twooverlap1$sourceId,twooverlap1$REF[match(bed$name,twooverlap1$sourceId)],NA)
  bed$svaba=ifelse(bed$name %in% rownames(info(svaba)),1,0)
  bed$gridss=ifelse(bed$name %in% rownames(info(gridds)) | bed$name %in% rownames(a),1,0 )
  bed$manta=ifelse(bed$name %in% rownames(ms) | bed$name %in% rownames(gm),1,0)
  bed=bed[order(bed$chrom1,bed$start1),]
  bed=as.data.frame(bed)
  j=1
  while(j<(nrow(bed)-1)) {
    if((bed$chrom1[j] == bed$chrom1[j+1]) & (bed$start1[j+1] < (bed$start1[j] + maxgap1))) {
      # Calculate the sum for the current and next rows
      sum_current <- bed$gridss[j] + bed$manta[j] + bed$svaba[j]
      sum_next <- bed$gridss[j+1] + bed$manta[j+1] + bed$svaba[j+1]
      
      # Use filter to keep the row with the higher sum
      if (sum_current > sum_next) {
        bed <- bed[-(j+1),]
      } else {
        bed <- bed[-j,]
      }
      j=j-1
    }
    j=j+1
  }
 # Filter out rows where the name contains 'chr_Un'
  bed = bed[!grepl("chr_Un", bed$name), ]
  # Write the filtered data frame to a file
  write.table(as.data.frame(bed), file = paste('results/', sample_names[i], '_twooverlap_chr_1000bmaxgap.txt', sep=''), sep='\t', col.names=T, row.names = F,quote=FALSE) # provide your output file path here

}
