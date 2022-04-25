# AJPerez, 2022
library(pheatmap)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

sp <- "ec" # faltan OMPs para ec, ef y sa
type <- "if"
#MIN_STRAIN <- 1 # 33% (ab=400; pa=829 si ifs_only) ifa=111; ifb=303 i=1184 if=940
CLUSTERS <- 2 # 5 con pa, 2 con ab

setwd(paste0("./", sp))
FILE <- paste0("max/gnmatrix_", sp, ".tsv") # for memory RAM free: _membrane_omp.tsv <--------------------##############3
METADATA <- paste0("max/metadata_", sp, "_def.tsv") # metadata_crispr.tsv _def_extended.tsv #<-------
GENES <- paste0("membrane_omp_", type, ".id") # <---------- ###################3 vir_phages.id ifab_mlst_-07_07.id mbocps_ifa_st79_-03_03.id pan_omp_9090.id membrane.id pan_biofilm_9090.id signalp_ab.id conjugation.id plasmid.id 
VIRUS <- "virus_sma3s_db_btw1_0.id" # OR vir_phages.id OR vir_phagues3.id <------------------------------###################3
STRAINS <- paste0("types/", type, ".ab") # ifs_only.ab st1.ab(b) st79.ab(a) max/strains.ab crisprif.mlst.ab crisprif.mlst_noST2.ab ifs.ab ifs_complete.ab
COLUMN_ORDER <- 2 # 1=MLST, 2=I-Fs, 6=Nspacers
CLUSTERING <- "FALSE" # Row clustering

# Get gene mapping
mapping <- read.csv(paste0("roary/", sp, "2/pangenome_references_", sp, "_uniprot_bacteria_go.tsv"), sep = "\t", header = T)[,1:2]

# Matrix and Metadata
gnmatrix <- as.matrix(read.csv(FILE, header = T, sep = '\t', row.names = 1))
metadata <- read.csv(METADATA, sep = "\t", header = T, row.names = 1)
#metadata <- metadata[,c(2:5,8:8)] #<---------
metadata <- metadata[,c(1,2)]

#metadatag <- read.csv("pangenome_annotations.tsv", sep = "\t", header = T, row.names = 1) #<-------
#metadatag <- metadatag[,c(1:3,6)] #<---------

# Filter by genes
virus <- readLines(VIRUS)
gnmatrix2 <- subset(gnmatrix, !rownames(gnmatrix) %in% virus)
genes <- readLines(GENES)
gnmatrix2 <- subset(gnmatrix2, rownames(gnmatrix2) %in% genes)
gnmatrix2 <- gnmatrix2[order(row.names(gnmatrix2)), ]
mapping <- subset(mapping, mapping[,1] %in% rownames(gnmatrix2))
mapping <- mapping[order(mapping$X.ID), ]
rownames(gnmatrix2) <- mapping[,2]
gnmatrix2 <- t(gnmatrix2)

# Filter strains
strains <- readLines(STRAINS)
MIN_STRAIN <- length(strains)*0.1
gnmatrix2 <- subset(gnmatrix2, rownames(gnmatrix2) %in% strains)
metadata <- subset(metadata, rownames(metadata) %in% rownames(gnmatrix2))
#metadatag <- subset(metadatag, rownames(metadatag) %in% colnames(gnmatrix2)) #<-------

# Heatmap
gnmatrix2 <- gnmatrix2[,colSums(gnmatrix2) >= MIN_STRAIN]
gnmatrix2 <- gnmatrix2[order(metadata[,COLUMN_ORDER], metadata[,2]), ]
p <- pheatmap(gnmatrix2, show_rownames = F, show_colnames = T, treeheight_col = 30, legend = F, angle_col = 90,
         cluster_rows = T, cex = 1, treeheight_row = 20, annotation_row = metadata, #annotation_col = metadatag, # <-------
         fontsize = 8, color = c("#F9F2EA", "#804AD5"), #clustering_distance_rows = "manhattan",
         cutree_rows = CLUSTERS)
#print(p)
row_clusters <- data.frame(cutree(p$tree_row, k = CLUSTERS))
write.table(row_clusters, paste0("clusters_row_", type, ".tsv"), sep = "\t", quote = F, col.names = F)

# Groups/Clusters
clusters <- data.frame(cutree(p$tree_row, k=2))
colnames(clusters) <- c("Cluster")
clusters$Cluster <- as.character(clusters$Cluster)
metadata2 <- merge(clusters, metadata, by='row.names',all=TRUE)
rownames(metadata2) <- metadata2[[1]]
metadata2 <- metadata2[-1]
metadata2 <- metadata2[,c(3,2,1)]
p2 <- pheatmap(gnmatrix2, show_rownames = F, show_colnames = T, treeheight_col = 0, legend = F, angle_col = 90,
               cluster_rows = T, cex = 1, treeheight_row = 0, annotation_row = metadata2, #annotation_col = metadatag, # <-------
               fontsize = 8, color = c("#F9F2EA", "#804AD5"), cutree_rows = CLUSTERS)

tiff(paste0("/home/ajperez/Documentos/Articulos/CRISRPRomESKAPE/Figures/heatmap_", sp, "_", type, ".tif"), width = 6, height = 6, units = "in", compression = "lzw", res = 300)
print(p2)
dev.off()
