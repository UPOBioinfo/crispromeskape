# crispromeskape

### Scripts for the paper "An outstanding bacterial membranome requires CRISPR-Cas systems to avoid the intervention of phages"

-- add doi of arXiv

Add some text to describe the project. 


## Python scripts
Package versions:
matplotlib 3.5.0
numpy 1.19.2
pandas 1.3.5
scikit-learn 1.0.2
scipy 1.6.2
seaborn 0.11.2
### Random forest inference
The script utilizes scikit-learns mean decrease in impurity feature importance measre to infer genes that set crispr containing strains apart from those without crispr systems. Cas-genes have been removed in the data, to allow to detect non Crispr-Cas related genes. 
### Jaccard distances
This script computes the Jaccard dissimilarity between the MLST groups which is used as distance measurement for a hierarchical clustering and dendrograms available in the papers supplementary data.
### Correlation plot
This script uses the Jaccard indices of the jaccard-panel.py script, so the latter needs to be run beforehand. It plots the correlation plot featured in the paper. 

## Perl scripts
**addColumnFromABfile.pl**  
Create an output to add to the metadata table, from a file with 2 columns: metadata ID.

**addColumnsFromABfile.pl**  
Creates an output with several columns to add to the metadata table, from a list of files, each with a list of IDs.

**analyseSpacers.pl**  
Calculate proportions of genomes with unique spacer/phages in two different clusters.

**calculatePercentGenesInStrainsID.pl**  
Calculate both absolute and relative frequency of pangenome genes in two clusters of genomes.

**countGenesAndShareGenesFromPangenome.pl**  
Count genes and average number of shared genes for each genome in the pangenome.

**createComparisonMatrix.pl**  
Create matrix for heatmaps comparing genes vs genomes.

**createMatrixPresenceAbsence.pl**  
Create gene presence/absence matrix.

**extractSpacers.pl**  
Extract spacers with evidence code 4 from CRISPRCasFinder results.

**filterBlast_vs_db.pl**  
Filter Blast hits with a given identity and subject coverage.

**getGroupsFromPangenome.pl**  
Count appearances of each gene in the pangenome.

**joinGFFwithCRISPR_v2.pl**  
Prepare contigs from each genome as a string of gene names.

**putSma3sGenenames2roaryClusters.pl**  
Modify gene names assigned by Sma3s to avoid redundant names.

**removeAssemblies.pl**  
Check assembly files and remove those from GenBank when RefSeq is available. 
