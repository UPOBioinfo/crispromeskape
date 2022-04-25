library(ggplot2)
library(ggtree)
library(cowplot)
library(ape)
library(phytools)
library(treeio)

setwd("C:/Users/Usuario/Desktop/phylo")

tree<-read.raxml("RAxML_bipartitionsBranchLabels.SEQMLST_red") 

tree@phylo$tip.label

tree@data$bootstrap

roottree <- root(tree, "ab00001", node=NULL, resolve.root=FALSE)

heatmap_tabla <- read.table("clade_fig6_metadata.tsv", sep = "\t", header = T, quote = "", 
                            row.names = 1, stringsAsFactors = FALSE)

p <- ggtree(roottree, layout = 'circular', size = 0.5) + geom_point2(aes(label=bootstrap, subset=bootstrap>70), size=1.5, color = "red")

phylo2 <- gheatmap(p, heatmap_tabla, offset = -0.001,  width = 0.09, font.size = 0, colnames_angle = 90, hjust = 0, colnames = F) +
  scale_fill_brewer(palette = "Set1", name = "") + theme(legend.text  = element_text (size = 20))

pdf("fig6.pdf", width =20, height=10, paper='special')
print(phylo2)
dev.off()
