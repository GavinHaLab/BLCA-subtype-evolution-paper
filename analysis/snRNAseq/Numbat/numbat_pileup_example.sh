#!/bin/bash 

### RUN THIS FIRST: micromamba activate cellsnp-lite-env

# Load modules
module load SAMtools/1.19.2-GCC-13.2.0
module load R/4.4.0-gfbf-2023b

# Set up PATH for Micromamba and Eagle2
export PATH=~/micromamba/envs/numbat_env/bin:$PATH
export PATH=Numbat/Eagle_v2.4.1/eagle:$PATH

# Change to the directory with the R script or environment
cd /snRNAseq/Numbat

# Define variables specific to patient
BAM_FILES=(
"ML_16_079H1_possorted_genome_bam.bam"
"ML_16_079N3_possorted_genome_bam.bam"
"ML_16_079P4_possorted_genome_bam.bam"
)
BARCODE_FILES=(
"ML_16_079H1/barcodes.tsv.gz"
"ML_16_079N3/barcodes.tsv.gz"
"ML_16_079P4/barcodes.tsv.gz"
)

# Set other required variables (shared across patients)
GMAP="Numbat/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
EAGLE="Numbat/Eagle_v2.4.1/eagle"
SNPVCF="Numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
PANELDIR="Numbat/1000G_hg38"
OUTDIR_BASE="Numbat/results/merged_pileup_and_phase"
NCORES=16
UMITAG="None"
CELLTAG="None"

# Define the patient ID and label
PATIENT_ID="16_079"
LABEL="16_079_test1"
OUTDIR="${OUTDIR_BASE}/${PATIENT_ID}"

# Ensure the output directory exists
mkdir -p "$OUTDIR"

# Convert arrays to comma-separated strings
BAMS=$(IFS=,; echo "${BAM_FILES[*]}")
BARCODES=$(IFS=,; echo "${BARCODE_FILES[*]}")
SAMPLES="16_079H1,16_079N3,16_079P4"

# Log the key variables for Pt32
echo "Running Numbat for ${PATIENT_ID} with the following settings:"
echo "Label: $LABEL"
echo "BAM files: $BAMS"
echo "Barcodes: $BARCODES"
echo "Samples: $SAMPLES"
echo "Genetic Map: $GMAP"
echo "Eagle2 Path: $EAGLE"
echo "SNP VCF: $SNPVCF"
echo "Phasing Reference Panel: $PANELDIR"
echo "Output Directory: $OUTDIR"
echo "Number of Cores: $NCORES"
echo "UMI Tag: $UMITAG"
echo "Cell Tag: $CELLTAG"

# Run the pileup_and_phase.R script for Pt32 with the modified samples argument
Rscript Numbat/renv/library/linux-ubuntu-bionic/R-4.4/x86_64-pc-linux-gnu/numbat/bin/pileup_and_phase.R \
    --label $LABEL \
    --samples "$SAMPLES" \
    --bams "$BAMS" \
    --barcodes "$BARCODES" \
    --gmap "$GMAP" \
    --eagle "$EAGLE" \
    --snpvcf "$SNPVCF" \
    --paneldir "$PANELDIR" \
    --outdir "$OUTDIR" \
    --ncores "$NCORES"

# Notify completion
echo "Script completed for 16-079" | mail -s "Notification: Numbat Script Finished for 16-079" xxx@xx.org
