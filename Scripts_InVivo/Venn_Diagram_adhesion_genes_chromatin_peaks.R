#!/usr/bin/Rscript

#Robin W. Yeo
#02_22_2021

#This script generates Venn diagrams to visualize the overlap of:
#i) shared cell adhesion genes with nearby chromatin peaks that open in young qNSCs (vs old qNSCs) and in old aNSCs (vs young aNSCs)
#ii) shared cell adhesion genes with nearby chromatin peaks that open in old aNSCs and old NPCs

rm(list=ls())
library(VennDiagram)
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")

qNSC <- read.table(file="Original_Data/Differential_Peaks/Differential_ATAC_peaks_adhesion_GO_qNSC_Y_vs_O.txt", header=TRUE, sep="\t", quote="")
aNSC <- read.table(file="Original_Data/Differential_Peaks/Differential_ATAC_peaks_adhesion_GO_aNSC_Y_vs_O.txt", header=TRUE, sep="\t", quote="")

qNSC_YOUNG <- qNSC[qNSC$Fold>0,]
aNSC_OLD <- aNSC[aNSC$Fold<0,]

qNSC_YOUNG_genes <- as.data.frame(unique(qNSC_YOUNG$SYMBOL))
aNSC_OLD_genes <- as.data.frame(unique(aNSC_OLD$SYMBOL))

venn.diagram(
  x = list(unique(qNSC_YOUNG$SYMBOL), unique(aNSC_OLD$SYMBOL)),
  category.names = c("qNSC","aNSC"),
  filename = 'Output_Figs/venn_diagram_adhesion_SYMBOL_Y_qNSC_O_aNSC.jpg',
  output=TRUE,
  
  lwd = 2,
  fill = c("dodgerblue","red2"),
  
  cex = 1.6,
  fontface = "bold",
  fontfamily = "sans"
)


NPC <- read.table(file="Original_Data/Differential_Peaks/NPC_Old_Differential_Sites_withFC_FDR_0.05.narrowPeak", header=TRUE, sep="\t", quote="")
GO <- read.table(file="Original_Data/GO_Cell_Adhesion_GO0007155.txt",sep="\t",header=TRUE)
NPC <- NPC[NPC$SYMBOL %in% GO$Gene.name,]
NPC_OLD <- NPC[NPC$Fold<0,]

venn.diagram(
  x = list(unique(aNSC_OLD$SYMBOL), unique(NPC_OLD$SYMBOL)),
  category.names = c("aNSC","NPC"),
  filename = 'Output_Figs/venn_diagram_adhesion_aNSC_NPC_SYMBOL.jpg',
  output=TRUE,
  
  lwd = 2,
  fill = c("red2","purple1"),
  
  cex = 1.6,
  fontface = "bold",
  fontfamily = "sans"
)
