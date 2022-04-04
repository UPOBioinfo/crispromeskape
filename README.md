# crispromeskape

### Scripts for the paper "An outstanding bacterial membranome requires CRISPR-Cas systems to avoid the intervention of phages"

-- add doi of arXiv

Add some text to describe the project. 


## Python scripts
### Random forest inference
The script utilizes scikit-learns mean decrease in impurity feature importance measre to infer genes that set crispr containing strains apart from those without crispr systems. Cas-genes have been removed in the data, to allow to detect non Crispr-Cas related genes. 
### Jaccard distances
This script computes the Jaccard dissimilarity between the MLST groups which is used as distance measurement for a hierarchical clustering and dendrograms available in the papers supplementary data.
### Correlation plot
This script uses the Jaccard indices of the jaccard-panel.py script, so the latter needs to be run beforehand. It plots the correlation plot featured in the paper. 
