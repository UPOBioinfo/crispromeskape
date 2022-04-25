library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

setwd("./")
fig <- list()

for (PRE in c("ab_ifa", "ab_ifb", "ec_if", "ef_iia", "kp_ie", "pa_if", "pa_ie", "pa_ic", "sa_iiia", "kp_iva3", "pa_iva1")) {
  df <- read.table(paste0("./", PRE, "_0_suppl.tsv"), header = T, nrows = 5)
  df$Strains <- factor(df$Strains, level = df$Strains)

  label1 <- paste0("Cluster 1\n(", sum(df$cluster_1), ")")
  label2 <- paste0("Cluster 2\n(", sum(df$cluster_2), ")")
  
  df1 <- df %>% select(Strains, cluster_1, cluster_2) %>% gather(Group, Frequency, cluster_1, cluster_2)
  getPalette = c("#9F9F9F","#753FC3","#009CBF","#36D85F","#E74545")
  out1 <- ggplot(df1, aes(x = Group, y = Frequency, fill = Strains)) +
    geom_bar(position="fill", stat="identity", width = 0.8) + 
    ylab("Strains (%)") + 
    annotate("text", x = c(1, 2), y = 1.05, label = c(label1, label2), lineheight = 0.7) +
    scale_fill_manual("CRISPR-Cas type", values = getPalette) + 
    scale_y_continuous(labels = scales::percent) + theme_minimal() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
          axis.ticks.y=element_blank(), legend.position = "none")
 
  out <- ggarrange(out1, font.label = list(size = 16, color = "black", face = "bold"), nrow = 1, ncol = 1)  
  fig[[PRE]] <- out
}

# Final figure
leg <-
  ggplot(df1, aes(x = Group, y = Frequency, fill = Strains)) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(values = getPalette, name = "Legend") +
  theme(legend.key.size = unit(1, 'cm'), legend.text = element_text(size = 16), legend.title = element_text(size = 18))
my_legend <- get_legend(leg)
legend <- as_ggplot(my_legend)
empty <- ggplot() + theme_minimal()

out <- ggarrange(fig[["kp_ie"]], fig[["ab_ifa"]], fig[["ab_ifb"]],
                fig[["pa_ic"]], fig[["pa_ie"]], fig[["pa_if"]],
                fig[["ec_if"]], fig[["ef_iia"]], fig[["sa_iiia"]],
                fig[["kp_iva3"]], fig[["pa_iva1"]], legend,
                labels = c("K. pneumoniae (I-E)", "A. baumannii (I-Fa)", "A. baumannii (I-Fb)", 
                           "P. aeruginosa (I-C)", "P. aeruginosa (I-E)", "P. aeruginosa (I-F)",
                           "E. cloacae (I-F)", "E. faecium (II-A)", "S. aureus (III-A)",
                           "K. pneumoniae (IV-A3)", "P. aeruginosa (IV-A1)", ""),
                vjust = 1.02, hjust = -0.1, font.label = list(size = 12, color = "#000276", face = c("bold.italic")), nrow = 2, ncol = 6)
print(out)
tiff("Figures/fig5.tiff", width = 14, height = 9, units = "in", compression = "lzw", res = 300)
print(out)
dev.off()

exit()

ggarrange(fig[["kp_ie"]], fig[["ab_ifa"]], fig[["ab_ifb"]],
          fig[["pa_ic"]], fig[["pa_ie"]], fig[["pa_if"]],
          labels = c(
            "K. pneumoniae (I-E)", "A. baumannii (I-Fa)", "A. baumannii (I-Fb)", 
            "P. aeruginosa (I-C)", "P. aeruginosa (I-E)", "P. aeruginosa (I-F)"
          ), vjust = 0.7, hjust = -0.1,
          font.label = list(size = 14, color = "#000276", face = c("bold.italic")), nrow = 2, ncol = 3)

ggarrange(fig[["ec_if"]], fig[["ef_iia"]], fig[["sa_iiia"]],
          fig[["kp_iva3"]], fig[["pa_iva1"]],
          labels = c(
            "E. cloacae (I-F)", "E. faecium (II-A)", "S. aureus (III-A)", 
            "K. pneumoniae (IV-A3)", "P. aeruginosa (IV-A1)"
          ), vjust = 0.7, hjust = -0.1,
          font.label = list(size = 14, color = "#000276", face = c("bold.italic")), nrow = 2, ncol = 3)
