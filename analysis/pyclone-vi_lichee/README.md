# Pyclone-VI and LICHeE Analysis

This directory contains custom scripts for generating, running pyclone-VI and LICHeE tools. 
Original tool github : https://github.com/Roth-Lab/pyclone-vi and https://github.com/viq854/lichee/tree/master

## Overview
- Steps to run Pyclone-Vi and LICHeE
- Scripts used in the analysis which are avaiable for reference/use

## Key Scripts
- `masterScript_Pyclone_git.py`: Used for generating the inputs for pylcone-VI, requires per sample
- `Metastatic_seeding_analysis.Rmd`: Parses and summarizes seeding results for integration with clinical and genomic data.

##Resources
## Example: Running MACHINA on Sample `19-001`
This repository includes an example folder, `19-001/`, containing output files generated using the script `clone_tree_generation_for_machina_BLCA.R`.
Input files and expected results are provided in the `19-001_input_example/` directory.

### MACHINA Commands

**pmh**
```bash
pmh -p P001PriFFPE \
    -c 19_001/19_001_colorfile.tsv \
       19_001/19_001_tree.tsv \
       19_001/19_001_labeling.tsv \
    -o 19_001/pmh/ > 19_001/pmh_result.txt
**pmh_tr**
pmh_tr -p P001PriFFPE \
       -c 19_001/19_001_colorfile.tsv \
          19_001/19_001_tree.tsv \
          19_001/19_001_labeling.tsv \
       -o 19_001/pmh_tr/ > 19_001/pmh_tr_result.txt

**DOT to PDF Conversion**
MACHINA generates migration tree visualizations as .dot files in the results/ directory. To convert these files into PDF format for easier viewing and sharing, use the following command:
dot -Tpdf sampleName-S.dot -o sampleName-S.pdf
