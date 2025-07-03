Steps to run Numbat:
user guide: https://kharchenkolab.github.io/numbat/
1- numbat_pileup_example.sh: example of running the first pileup and phase step of Numbat
2- Numbat_Step2.R: run Numbat CNA calling from pileup output

Numbat analysis scripts:
 -Numbat_analysis.Rmd: graphing tumor vs normal clones
 -Numbat_GeneSpecific_makeDFs.Rmd: make dataframes of clonality status for individual genes
 -Numbat_GeneSpecific.Rmd: use above dataframes and graph/analyze clonality

Preparing TITAN bulk CNA calls for Numbat input:
1- Numbat_TITAN_CNV_to_Numbat_input.Rmd: put TITAN segs.txt output file in the Numbat input format
2- Numbat_combine_segs.Rmd: for making a "union" CNA profile from multiple tumors, after converting to the Numbat input format above
