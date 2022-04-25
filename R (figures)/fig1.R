library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

setwd("./")

colores <- c('#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F', '#94D4C0','#FDAF91','#AFBDDB','#EEADD5','#C1E487','#FFE56D')
name_colors <- c("#01AA18","#01AA18","#BF0093","#BF0093","#BF0093","#BF0093")

df <- as.data.frame(readxl::read_excel("SupplFiles/SupplTableS1.xlsx", sheet = "Summary"))
df$Bacteria <- factor(df$Bacteria, levels = df$Bacteria)

df1 <- df %>% select(Bacteria, Strains, Strain_crispr) %>% gather(Total, Value, -Bacteria)
fa <- ggplot(df1, aes(Bacteria, Value, fill = rev(interaction(Bacteria, Total)))) +
  geom_col(position = "dodge") +geom_text(aes(label = Value), vjust = -0.3, size = 5, position = position_dodge(width = .9)) +
  scale_fill_manual(values = colores) + ggtitle("Total number of genomes / Number of genomes with CRISPR-Cas    ") + #scale_fill_brewer(palette = "Set2") + 
  theme_minimal() + theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, color = name_colors), legend.position = "none", axis.title.x=element_blank(),
                          plot.title = element_text(hjust = 0.5, size = 18), text = element_text(size = 20)) +
  ylab("Number of genomes") # + theme(legend.position = "none")
print(fa)

df2 <- df %>% select(Bacteria, Pangenome, Core_99) %>% gather(Total, Value, -Bacteria)
fb <- ggplot(df2, aes(Bacteria, Value, fill = rev(interaction(Bacteria, Total)))) +
  geom_col(position = "dodge") + geom_text(aes(label = Value), vjust = -0.3, size = 5, position = position_dodge(width = .9)) +
  scale_fill_manual(values = colores) + ggtitle("Number of pangenome genes / Number of core genes") + #  scale_fill_brewer(palette = "Set2") +
  theme_minimal() + theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, color = name_colors), legend.position = "none", axis.title.x=element_blank(),
                          plot.title = element_text(hjust = 0.5, size = 18), text = element_text(size = 20)) +
  ylab("Number of genes") # + theme(legend.position = "none")
print(fb)

df3 <- df %>% select(Bacteria, Ngenes_desv, Nshared_desv, Ngenes_avg, Nshared_avg) %>% gather(Group, Value, Ngenes_avg, Nshared_avg, -Bacteria)
fc <- ggplot(df3, aes(Bacteria, Value, fill = interaction(Bacteria, Group))) +
  geom_col(position = "dodge") + geom_text(aes(label = Value), vjust = -1.7, size = 5, position = position_dodge(width = .9)) +
  scale_fill_manual(values = rev(colores)) + ggtitle("Average number of genes / Average number of shared genes") + #scale_fill_brewer(palette = "Set2") + 
  theme_minimal() + theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, color = name_colors), legend.position = "none", axis.title.x=element_blank(),
                          plot.title = element_text(hjust = 0.5, size = 18), text = element_text(size = 20)) + # + theme(legend.position = "none")
  geom_errorbar(aes(ymin = Value - Ngenes_desv, ymax = Value + Ngenes_desv), width=.2, position = position_dodge(.9)) +
  ylab("Number of genes")
print(fc)

df4 <- as.data.frame(readxl::read_excel("SupplFiles/SupplTableS1.xlsx", sheet = "CRISPR"))
colourCount = nrow(unique(df4 %>% select(Group)))
df4$Bacteria <- factor(df4$Bacteria, level = df$Bacteria)
getPalette = c("darkred","#DA5D5D","#FCA9B0","darkgreen","#87DE5E","orange","darkblue","#54C8FF","#A5017F","#DC4CBB","#F4A6E2","#9F9F9F")
fd <- ggplot(df4, aes(x = Bacteria, y = Frequency, fill = Group)) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual("CRISPR-Cas", values = getPalette) + theme_minimal() +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, color = name_colors), text = element_text(size = 20),
        legend.text = element_text(size = 18), legend.key.size = unit(1, 'cm'))
print(fd)

df5 <- as.data.frame(readxl::read_excel("SupplFiles/SupplTableS1.xlsx", sheet = "spacers"))
df5$Bacteria <- factor(df5$Bacteria, level = df5$Bacteria)
df5b <- df5 %>% select(Bacteria, Phage, Plasmid, Phage_Plasmid, Other, No_hit, Nspacers) %>% 
  gather(Group, Proportion, Phage, Plasmid, Phage_Plasmid, Other, No_hit)
getPalette = c("#9F9F9F","darkgreen","darkred","#944DFF","#10ACFF","#E90FF0","#AA0083","#780043","#9F9F9F")

label1 <- paste0(df5b$Bacteria[seq(1,6)], "\n(", df5b$Nspacers[seq(1,6)], ")")
fe <- ggplot(df5b) +
  geom_bar(aes(x = forcats::fct_rev(Bacteria), y = Proportion, fill = Group), position="fill", stat="identity", width = 0.8) +
  coord_flip() + scale_fill_manual("Protospacers", values = getPalette) + theme_minimal() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L), breaks = seq(0, 1, .1)) +
  scale_x_discrete(labels = rev(forcats::fct_rev(label1))) +
  ylab("Proportion of spacers") + xlab("Species") +
  theme(axis.text.y = element_text(face = "italic", hjust = 1, color = rev(name_colors)), text = element_text(size = 20), 
        legend.title = element_blank(), legend.position="top", legend.text = element_text(margin = margin(r = 0, unit = 'cm')),
        axis.text.x=element_text(angle = 45, hjust=1), legend.justification = "right") + guides(fill = guide_legend(reverse = TRUE))

df6 <- as.data.frame(readxl::read_excel("SupplFiles/SupplTableS1.xlsx", sheet = "virusplasmid"))
df6$Bacteria <- factor(df6$Bacteria, level = df6$Bacteria)
df6b <- df6 %>% select(Bacteria, Phage, Plasmid, Phage_Plasmid, Other, Nspacers) %>% 
  gather(Group, Proportion, Phage, Plasmid, Phage_Plasmid, Other)
getPalette = c("darkgreen","darkred","#944DFF","#10ACFF","#E90FF0","#AA0083","#780043","#9F9F9F")

label1 <- paste0(df6b$Bacteria[seq(1,6)], "\n(", df6b$Nspacers[seq(1,6)], ")")
ff <- ggplot(df6b) +
  geom_bar(aes(x = forcats::fct_rev(Bacteria), y = Proportion, fill = Group), position="fill", stat="identity", width = 0.8) +
  coord_flip() + scale_fill_manual("Type of genes", values = getPalette) + theme_minimal() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L), breaks = seq(0, 1, .1)) +
  scale_x_discrete(labels = rev(forcats::fct_rev(label1))) +
  ylab("Proportion of gene types") + xlab("Species") +
  theme(axis.text.y = element_text(face = "italic", hjust = 1, color = rev(name_colors)), text = element_text(size = 20), 
        legend.position="top", legend.title = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust=1)) + guides(fill = guide_legend(reverse = TRUE))

# Figure
out <- ggarrange(fa, fb, fc, fd, ff, fe,
                 labels = c("a)", "b)", "c)", "d)", "e)", "f)"), 
                 font.label = list(size = 20, color = "black", face = "bold"), nrow = 2, ncol = 3)

tiff("Figures/fig1.tiff", width = 24, height = 14, units = "in", compression = "lzw", res = 300)
print(out)
dev.off()



