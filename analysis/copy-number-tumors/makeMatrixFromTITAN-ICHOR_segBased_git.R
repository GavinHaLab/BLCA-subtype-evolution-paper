library(stringr)
library(data.table)
setDTthreads(threads = 2)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#install.packages(GenomicRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(dplyr)
source("/fh/fast/ha_g/user/gha/software/git/TitanCNA/R/utils.R")

args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE, width=180)

inDir <- args[1] #list of files for input
geneFile <- args[2] #references/GRCh38/gencode/GRCh38.p12.ensembl.gene.annotations.sorted.txt
sampleFile <- args[3]

print(sampleFile)
#method <- args[4] #{'common','severity','complete'}; default is 'severity'
method <- "common"
#headerType <- args[5] #{'Gene','chrPosn'}
headerType <- "Gene"
#filterLen <- as.numeric(args[6])
filterLen <- as.numeric(1000)
outPrefix <- args[7]
outImage <- paste(outPrefix,"_geneMats.RData",sep="")

#cnExt <- ".titan.ichor.seg.txt"
cnExt <- ".titan.ichor.seg.noSNPs.txt"
cnExt <- ".titan.ichor.seg.noSNPs_Modified.txt"  
overlapType <- "any"
cnCol <- "Corrected_Copy_Number"
callCol <- "Corrected_Call"
filterSnp <- 0
genomeStyle <- "UCSC"
genomeBuild <- "hg38"
print("here after genomeBuild")
seqinfo <- Seqinfo(genome=genomeBuild)
print(" after genomeBuild")
chrs <- c(1:22,"X")
seqlevelsStyle(chrs) <- genomeStyle
#print("here seqlevel style")
#print((chrs))
# if the sample is not found then skip to the next sample

samples <- fread(sampleFile)
print(samples)
genes <- fread(geneFile)
print("error gene")
if (headerType == "Gene"){
	setnames(genes, c("cdsStart", "cdsEnd"), c("Start", "End"))
	genes <- genes[, .(Chr = Chr[1], Start = min(Start), End = max(End), Karyotype_band=Karyotype_band[1], strand=strand[1]), by=Gene]
	genes[, Length := End - Start]
	genes <- genes[Length > 1]
	genes[, chrPosn := paste0(Chr,":", Start, "-", End)]
}else{ #if (!headerType %in% colnames(genes)){
   colnames(genes)[1:3] <- c("Chr", "Start", "End")
}
seqlevelsStyle(genes$Chr) <- genomeStyle
genes <- genes[Chr %in% chrs]
genes$Chr <- factor(genes$Chr, levels = chrs)
genes <- genes[order(Chr, Start)]
numGenes <- nrow(genes)
#setkey(genes, chrPosn)
#print(numGenes)
stateKey <- array(0:9); names(stateKey) <- c('HOMD','DLOH','HET','NLOH','GAIN','ALOH','BCNA','UBCNA','ASCNA','HLAMP')
severity <- array(c(6,5,4,5,1,2,3,3,4,5)); names(severity) <- c('HOMD','DLOH','NLOH','ALOH','HET','GAIN','BCNA','UBCNA','ASCNA','HLAMP')
print("error prints")
save.image(outImage)
#########################################################################################
####################### FUNCTION: FIND OVERLAP OF SEGMENT AND GENE ######################
#########################################################################################"

## input: lohHits = loh rows that overlap region of interests
# start = start coordinate of region of interest
# end = end coordinate of region of interest
## NOT USED
getOverlapLength <- function(lohHits, start, end){
	coords <- cbind(lohHits[, c("Start","Stop")], as.numeric(start), as.numeric(end))
	coordsSort <- t(apply(coords, 1, sort))
	dist <- coordsSort[, 3] - coordsSort[, 2] + 1
	return(dist)
}

#Input: states is an array of state names (e.g. (DLOH,NLOH,...,))
#Uses global variable "severity"
## NOT USED
getMostSevereState <- function(states){	
	severityValue <- 0
	severeState <- states[1]
	for (i in states){
		if (severity[i] > severityValue){
			severeState <- i
			severityValue <- severity[i]
		}
	}
	return(severeState)
}

# output the matrix to file
#Input: Matrix to output; output file name
writeMatrixToFile <- function(mat,outfile){
	outMat <- cbind(rownames(mat),mat)
	if (!is.null(colnames(outMat))){
		colnames(outMat)[1] <- "Sample"
	}
	write.table(outMat,file=outfile,row.names=F,col.names=T,quote=F,na="NaN",sep="\t")
}
print("error")
## finds the bin-level overlap
# bins from CN data that overlap a gene 
getBinCNOverlap <- function(x, y, func = "mean", type = "any", colToReturn = "Copy_Number"){
	cn <- rep(NA, nrow(x))
	x.gr <- as(x, "GRanges")
	y.gr <- as(y, "GRanges")
	hits <- findOverlaps(query = y.gr, subject = x.gr, type = type)
	hits.dt <- as.data.table(hits)
	avgBinCN <- hits.dt[, match.fun(func)(y[queryHits, round(get(colToReturn))], na.rm=T), by = subjectHits]
	cn[avgBinCN$subjectHits] <- avgBinCN$V1
	return(cn)
}


## finds the bin-level overlap
# bins from CN data that overlap a gene 
getBinCNOverlap2 <- function(x, y, func = "mean", type = "any", colToReturn = "Corrected_MajorCN"){
        cn <- rep(NA, nrow(x))
        x.gr <- as(x, "GRanges")
        y.gr <- as(y, "GRanges")
        hits <- findOverlaps(query = y.gr, subject = x.gr, type = type)
        hits.dt <- as.data.table(hits)
        avgBinCN <- hits.dt[, match.fun(func)(y[queryHits, round(get(colToReturn))], na.rm=T), by = subjectHits]
        cn[avgBinCN$subjectHits] <- avgBinCN$V1
        return(cn)
}

print(inDir)
files <- list.files(inDir, pattern = cnExt, full.names = TRUE)
message("files",files)
ids <- str_extract(basename(files), "[0-9A-Za-z-_]+")
message("ids ",ids)
names(files) <- ids
files <- files[samples[[1]]]
message("files ",files[1:2])
numSamples <- length(files)
message("numSamples ",numSamples)
numGenes <- nrow(genes)


geneCallmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneCNmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneLOHmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneARmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneHRmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneLogRmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneCFmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneCCmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneLenmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneBinCNmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneMajorCNmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
#newly added
geneCorrectedCNMat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneCorrectedCallMat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))

#print(geneBinCNmat)
for (i in 1:numSamples){
	caseId <- names(files[i])
	cat("Analyzing sample: ",caseId," for file: ",files[i],"...\n")
	cat("files",files[i])
	# Check if the file exists added on 10/22/2024
    if (!file.exists(files[i])) {
        cat("No data found for sample: ", caseId, " for file: ", files[i], "\n")
        next
    }
	print(c("Case",caseId))
	#LOH
	loh <- fread(files[i])
	#print(loh)
	colnames(loh)[c(2,3,4)] <- c("Chr","Start","End")
	loh[, Length.bp := End - Start + 1]
	## filter by length threshold ##
	#print(loh)
	loh <- loh[Length.bp>=filterLen & Length.snp.>=filterSnp, ]
	loh[, LOH := as.numeric(MinorCN == 0)]
        #print(files[i])	
	# CNA bin file to use for chrX CN
	bin <- fread(gsub("seg.txt", "cna.txt", files[i]))
	bin2 <- fread(gsub("seg.txt", "seg.txt", files[i]))
	print(bin2)
	#binChrXind <- bin[Chr == "chrX", which = TRUE]
	#geneChrXind <- genes[Chr == "chrX", which = TRUE]
	geneBinCN <- getBinCNOverlap(x=genes, y=bin, func = "mean", type = "any", colToReturn = "logR_Copy_Number")
	geneMajorCN <- getBinCNOverlap2(x=genes, y=bin2, func = "mean", type = "any", colToReturn = "Corrected_MajorCN")
	#geneMajorCN <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Corrected_MajorCN")
	#convert to gene scaffold
	geneCN <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn=cnCol)
	geneLOH <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="LOH")
	geneLogR <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Median_logR")
	geneAR <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Median_Ratio")
	if ("Median_HaplotypeRatio" %in% colnames(loh)){
    	geneHR <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Median_HaplotypeRatio")
  	}
	geneCF <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Cellular_Prevalence")
	geneCC <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Clonal_Cluster")
	geneState <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="TITAN_state")
	geneLen <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Length.bp")
	#geneMajorCN <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Corrected_MajorCN")
	geneCorrectedCN <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Corrected_Copy_Number")
	geneCall <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Corrected_Call")
	#build matrices
	geneCallmat[caseId,] <- geneState  #call matrix
	geneCNmat[caseId,] <- geneCN #call matrix
	geneLOHmat[caseId,] <- geneLOH #LOH matrix
	geneARmat[caseId,] <- geneAR   #AR matrix
	geneLogRmat[caseId,] <- geneLogR   #logR matrix
	geneCFmat[caseId,] <- geneCF  #cellular prevalence matrix
	geneCCmat[caseId,] <- geneCC  #clonal cluster matrix
	geneBinCNmat[caseId, ] <- geneBinCN # average bin-level CN 
	geneMajorCNmat[caseId, ] <- geneMajorCN # Major CN 
	geneCorrectedCNMat[caseId, ] <- geneCorrectedCN #corrected total copy number
	geneCorrectedCallMat[caseId, ] <- geneCall #corrected type of copy number - AMP, DEL etc

	if ("Median_HaplotypeRatio" %in% colnames(loh)){
		geneHRmat[caseId,] <- geneHR   #haplotypeRatio matrix
	}
	geneLenmat[caseId, ] <- geneLen
	save.image(outImage)
}
save.image(outImage)


outCall <- paste(outPrefix,"_geneCalls.txt",sep="")
writeMatrixToFile(geneCallmat,outCall)
outCN <- paste(outPrefix,"_geneCN.txt",sep="")
writeMatrixToFile(geneCNmat,outCN)
outLOH <- paste(outPrefix,"_geneLOH.txt",sep="")
writeMatrixToFile(geneLOHmat,outLOH)
outAR <- paste(outPrefix,"_geneAR.txt",sep="")
writeMatrixToFile(geneARmat,outAR)
outLogR <- paste(outPrefix,"_geneLogR.txt",sep="")
writeMatrixToFile(geneLogRmat,outLogR)
outCF <- paste(outPrefix,"_geneCF.txt",sep="")
writeMatrixToFile(geneCFmat,outCF)
outCC <- paste(outPrefix,"_geneCC.txt",sep="")
writeMatrixToFile(geneCCmat,outCC)
outHR <- paste(outPrefix,"_geneHR.txt",sep="")
writeMatrixToFile(geneHRmat,outHR)
outLen <- paste(outPrefix,"_geneLength.txt",sep="")
writeMatrixToFile(geneLenmat,outLen)
outBinCN <- paste(outPrefix,"_geneBinCN.txt",sep="")
writeMatrixToFile(geneBinCNmat,outBinCN)
outMajorCN <- paste(outPrefix,"_majorCN.txt",sep="")
writeMatrixToFile(geneMajorCNmat,outMajorCN)
outTotalCN <- paste(outPrefix,"_TotalCN.txt",sep="")
writeMatrixToFile(geneCorrectedCNMat,outTotalCN)
outCallType <- paste(outPrefix,"_CallType.txt",sep="")
writeMatrixToFile(geneCorrectedCallMat,outCallType)


save(geneCallmat, geneCNmat, geneLOHmat, geneARmat, geneLogRmat, geneCFmat, geneCCmat, geneHRmat, geneLenmat,geneCorrectedCNMat,geneCorrectedCallMat,file=outImage)
#save.image(outImage)
