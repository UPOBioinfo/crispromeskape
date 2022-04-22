library(topGO)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)

setwd("./")
mybreaks1 = c(1, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 0) #C
mybreaks2 = c(1, 1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 0) # P
mybreaks3 = c(1, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 0) #F

GO_PARAM1 <- "MF"     # CC | BP | MF
GO_PARAM2 <- "F"      # C | P | F
MYBREAKS <- mybreaks3 # mybreaks1 | mybreaks2 | mybreaks3
MAX <- 1e-7           # 1e-7(CC) | 1e-10(BP) | 1e-7(MF) Change mybreak123 above

fig <- list()
mylabels <- list()
for (pre in c("ab_I-Fa", "ab_I-Fb", "ec_I-F", "ef_II-A", "kp_I-E", "pa_I-F", "pa_I-E", "pa_I-C", "sa_III-A", "kp_IV-A3", "pa_IV-A1")) {
  p <- list()
  ab <- substr(pre, 1, 2)
  file_background <- paste0("../", ab, "/roary/", ab, "2/pangenome_references_", ab , "_uniprot_bacteria_go.tsv")
  Nodes <- 6 # number of processes to show
  INTERVAL <- 10
  
  for (F in c("v1")) {
    for (type in c("crispr")) { #, "nocrispr")) {
      # Files
      allRes <- data.frame()
      file_genes <- paste0(F, "/", pre, "_", type, ".id") # file_genes <- paste0(F, "/", pre, "_freq3_", type, ".id")
      i <- paste0(F, "_", type)

  ONT <- GO_PARAM1
  letter <- substr(ONT, 2, 2)
  Ontology <- paste0("GO.", GO_PARAM2, ".ID") #PFC (BP MF CC)
    
  #Create temp file
  data <- read.csv(file_background, sep = "\t", header = TRUE, row.names = NULL)[,(c('X.ID', Ontology))]
  data[[Ontology]] <- as.character(gsub(';', ', ', data[[Ontology]]))
  file_temp <- paste0(file_background,"2")
  write.table(data, file = file_temp, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  genes <- readLines(file_genes)

  # Get background annotation
  GOesByID <- readMappings(file = file_temp)
  bg_genes <- names(GOesByID)

  compared_genes <- factor(as.integer(bg_genes %in% genes))
  names(compared_genes) <- bg_genes
  
  # Create topGO object
  GOdata <- new("topGOdata", ontology = ONT, allGenes = compared_genes,
                annot = annFUN.gene2GO, gene2GO = GOesByID)
  asd <- unlist(Term(GOTERM))
  
  # Run Fisher test
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Create and print table with enrichment result
  allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = Nodes)
  
  # Different palettes
  palette <- c("#F52A2A", "#D561EA", "#61B0EA", "#47CD30", "#E89B57", "#E4EA61", "white") # alternative palette
  myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
  
  # Figure
  ########
  layout(t(1:2), widths=c(8,3))
  par(mar=c(4, .5, .7, .7), oma=c(3, 15, 3, 4), las=1)
  
  pvalue <- as.numeric(gsub("<", "", allRes$classicFisher)) # remove '<' symbols
  allRes$classicFisher <- pvalue
  max_value <- as.integer(max(-log(pvalue)))+1
  pv_range <- exp(-seq(max_value, 0, -1))
  allRes <- mutate(allRes, plot_id = paste(GO.ID, Term, sep = " - "))

###  
  pvalue <- as.numeric(gsub("<", "", allRes$classicFisher)) # remove '<' symbols
  allRes$classicFisher <- pvalue
  max_value <- as.integer(max(-log(pvalue)))+1
  pv_range <- exp(-seq(max_value, 0, -1))
  allRes = mutate(allRes, plot_id = paste(GO.ID, Term, sep = " - "))
  mylabels <- paste(allRes$GO.ID, "-",  asd[allRes$GO.ID])
###
  
  if ( max(allRes$Significant) < 10) { INTERVAL <- 2 }
  if ( max(allRes$Significant) == 1) { INTERVAL <- 1 }
  p[[i]] <- 
    ggplot(data=allRes, aes(x=reorder(plot_id, Significant), y = Significant)) +
    geom_bar(stat="identity", color="black", aes(fill=as.numeric(log(classicFisher))), size = 0.3)+
    geom_text(aes(label=mylabels), position=position_fill(vjust=0), hjust=0, fontface="bold", size = 6) +
    coord_flip() + 
    theme(panel.background = element_blank(), panel.grid.major.x = element_line(colour = "darkgrey", size = 0.75),
          panel.grid.minor.x = element_line(colour = "grey", size = 0.75), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none", plot.margin = unit(c(1, 0, 0, 0), "cm"),
          axis.ticks.x =element_blank(), axis.line.y=element_blank(), axis.text=element_text(size = 12)) +
    ylab("Number of genes") +
    guides(fill = guide_colourbar(barheight = 12, reverse=T)) +
    scale_fill_gradientn(name = "p-value", colours = palette, limits = log(c(MAX, 1)), breaks = log(MYBREAKS), # 1e-5 | 1e-12
                         guide = guide_colourbar(reverse = TRUE), labels = MYBREAKS) +
    scale_y_continuous(breaks = seq(0, max(allRes$Significant), by = INTERVAL)) # 10 | 5
  fig[[pre]] <- ggarrange(p[[i]])
}
}
}

leg <-
  ggplot(data=allRes, aes(x=reorder(plot_id, Significant), y=Significant) ) +
  geom_bar(stat="identity", color="black", aes(fill=as.numeric(log(classicFisher))), size = 0.3)+
  ylab("number of genes") + 
  guides(fill = guide_colourbar(barheight = 10, barwidth = 1.6, reverse=T)) +
  scale_fill_gradientn(name = "p-value", colours = palette, limits = log(c(MAX, 1)), breaks = log(MYBREAKS), # 1e-5 | 1e-12
                       guide = guide_colourbar(reverse = TRUE), labels=MYBREAKS)
my_legend <- get_legend(leg)
legend <- as_ggplot(my_legend)
empty <- ggplot() + theme_minimal()

out <- ggarrange(fig[["ef_II-A"]], fig[["sa_III-A"]], legend,
                 fig[["kp_I-E"]], fig[["ab_I-Fa"]], fig[["ab_I-Fb"]],
                 fig[["pa_I-C"]], fig[["pa_I-E"]], fig[["pa_I-F"]], 
                 fig[["ec_I-F"]], fig[["pa_IV-A1"]], fig[["kp_IV-A3"]],
                 labels = c(
                   "E. faecium (II-A)", "S. aureus (III-A)", "",
                   "K. pneumoniae (I-E)", "A. baumannii (I-Fa)","A. baumannii (I-Fb)",
                   "P. aeruginosa (I-C)","P. aeruginosa (I-E)", "P. aeruginosa (I-F)",
                   "E. cloacae (I-F)", "P. aeruginosa (IV-A1)", "K. pneumoniae (kp_IV-A3)"
                 ), vjust = 1.8, hjust = -0.1,
                 font.label = list(size = 16, color = "#000276", face = c("bold.italic")), nrow = 4, ncol = 3)

tiff("/home/ajperez/Documentos/Articulos/CRISRPRomESKAPE/Figures/fig3.tiff", width = 19, height = 12, units = "in", compression = "lzw", res = 300)
print(out)
dev.off()

