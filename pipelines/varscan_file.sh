#!/bin/bash
#SBATCH -N1 -n4-t 3-0

MAINDIR=/path/to/Muse/Directory

SAMPLENAME=Sample-Name1

NORMAL=/path/to/samples/Sample-Name1_normal.bam
TUMOR=/path/to/samples/Sample-Name1_tumor.bam

HG38FA=/path/to/HG38.fasta
SNP144=./dbsnp_144.hg38.vcf.gz
varscan=/path/to/MuSEInstallation/VarScan.v2.4.4.jar

resdir=$MAINDIR/varscan2/$SAMPLENAME
mkdir -p $resdir
cd $resdir
samtools mpileup -f $HG38FA -q 1 -B $NORMAL > ${SAMPLENAME}.normal.pileup
samtools mpileup -f $HG38FA -q 1 -B $TUMOR > ${SAMPLENAME}.tumor.pileup
java -jar $varscan somatic ${SAMPLENAME}.normal.pileup ${SAMPLENAME}.tumor.pileup ${SAMPLENAME} --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6 --min-var-freq 0.10 --min-freq-for-hom 0.75 --normal-purity 1.0 --tumor-purity 1.00 --p-value 0.99 --somatic-p-value 0.05 --strand-filter 0 --output-vcf
