# Use the custom script to generate LICHeE inputs
python3 analysis/pyclone-vi/pycloneToLichee.py data/19-001/19-001_output.tsv./data/19-001/ 19-001
# Run the LICHeE tool once we have the inputs required
/fh/fast/ha_g/user/pitagi/tools/lichee/LICHeE/release/lichee -build -i ./data/19-001//19-001_lichee_sSNV.txt -cp -sampleProfile -minRobustNodeSupport 2 -minClusterSize 2 -maxClusterDist 0.2 -minPrivateClusterSize 1 -e 0.1 -o ./data/19-001//19-001_lichee_trees.txt -dotFile  ./data/19-001//19-001_lichee_tree.dot -color -dot
dot -Tpdf  ./data/19-001//19-001_lichee_tree.dot -O
