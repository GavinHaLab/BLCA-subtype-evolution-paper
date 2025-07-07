# cfDNA Mutation Overlap Analysis

## Analysis Steps

1. **Generate BED files of tumor SNVs**  
   Use `CollectAllelicCounts.Rmd` to make BED files of tumor SNV loci and prepare for running `CollectAllelicCounts`.

2. **Run GATK CollectAllelicCounts**  
   Execute the pipeline to count cfDNA reads at tumor SNV loci using the `gatk.sh` script.

3. **Create DataFrames of cfDNA Counts**  
   Use `CollectAllelicCounts.Rmd` to read GATK CollectAllelicCounts output and make dataframes of cfDNA reads at tumor SNV loci.

4. **Analyze and Visualize Mutation Overlap**  
   Use `Overlap_Analysis.Rmd` to:
   - Graph mutation overlap
   - Correlate overlap with other metrics (e.g., TFx, purity, tumor volume, and CP)

---

## Additional Analysis

- **TFx_determinants.Rmd**  
  This script:
  - Counts liver metastases per patient  
  - Compares TFx values based on liver metastasis presence
