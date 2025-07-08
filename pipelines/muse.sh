#!/bin/bash
# Bash script to launch a job with one tumor and normal pair bam file for running MuSE caller. This script was used in a wrapper for running multiple samples on the HPC cluster. 
#SBATCH -N1 -n4 -t 3-0

MAINDIR=/path/to/Muse/Directory

SAMPLENAME=Sample-Name1

NORMAL=/path/to/samples/Sample-Name1_normal.bam
TUMOR=/path/to/samples/Sample-Name1_tumor.bam


HG38FA=/path/to/HG38.fasta
SNP144=/path/to/dbsnp_144.hg38.vcf.gz
MuSE=/path/to/MuSEInstallation/

resdir=$MAINDIR/muse/$SAMPLENAME
mkdir -p $resdir
cd $resdir
$MuSE call -O ${SAMPLENAME}  -f $HG38FA $TUMOR $NORMAL
$MuSE sump -I ${SAMPLENAME}.MuSE.txt -G -O ${SAMPLENAME}.vcf -D ${SNP144}
