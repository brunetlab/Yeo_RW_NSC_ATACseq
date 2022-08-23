#!/usr/bin/Rscript

#Robin W. Yeo
#01/24/2019

#This script reads in *.csv file containing significance and fold-change values for GO terms (found in Extended Data Tables) and plots them using ggplot2.

rm(list=ls())
library(ggplot2)

setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
theme <- theme(aspect.ratio =3/2, panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black",size=15),axis.text.y=element_text(colour="black",size=15),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

##Outputs young vs. old barplots with GO enrichments for qNSc and aNSC peaks
GO <- read.csv("Original_Data/GO_data_summary_Young_vs_Old_2022.csv")
qNSC <- GO[which(GO$Condition=="qNSC"),]
aNSC <- GO[which(GO$Condition=="aNSC"),]

pdf("Output_Figs/GO/invivo_qNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=qNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-6,5,1)) + geom_bar(stat="identity", fill=(qNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="qNSC GO Enrichment", y="", x="")
dev.off()

pdf("Output_Figs/GO/invivo_aNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=aNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-4,6,1)) + geom_bar(stat="identity", fill=(aNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="aNSC GO Enrichment", y="", x="")
dev.off()


##Outputs a barplot with GO enrichments for the top 1000 genes driving the positive and negative PC2 directions from Y_O_Q_A PCA
GO <- read.csv("Original_Data/GO_data_summary_Q_A_PC2.csv")

pdf("Output_Figs/GO/PC2_Y_O_Q_A_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=GO, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-16,6,2)) + geom_bar(stat="identity", fill=(GO$Color_Code), width=0.85) + theme + coord_flip() +labs(title="PC2 GO Enrichment", y="", x="")
dev.off()


##Outputs a barplot with IPA enrichments for old aNSC Canonical Pathways
IPA <- read.csv("Original_Data/IPA_data_summary_Old_aNSC.csv")

pdf("Output_Figs/GO/IPA_old_aNSC_canonical_pathways_plot.pdf", width=10, height=15)
p <- ggplot(data=IPA, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(0,5,1)) + geom_bar(stat="identity", fill=(IPA$Color_Code), width=0.85) + theme + coord_flip() +labs(title="Old aNSC IPA Enrichment", y="", x="")
dev.off()


##Outputs a barplot with GO enrichments for NPC peaks
NPC <- read.csv("Original_Data/GO_data_summary_NPC.csv")

pdf("Output_Figs/GO/invivo_NPC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=NPC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-5,5,1)) + geom_bar(stat="identity", fill=(NPC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="NPC GO", y="", x="")
dev.off()


##Outputs a barplot with GO enrichments for endothelial peaks
Endo <- read.csv("Original_Data/GO_data_summary_Endo.csv")

pdf("Output_Figs/GO/invivo_Endo_Y_vs_O_GO_plot.pdf", width=15, height=10)
p <- ggplot(data=Endo, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-7,7,1)) + geom_bar(stat="identity", fill=(Endo$Color_Code), width=0.85) + theme + coord_flip() +labs(title="Endo GO", y="", x="")
dev.off()


##Outputs a barplot with GO enrichments for astrocyte peaks
Ast <- read.csv("Original_Data/GO_data_summary_Ast.csv")

pdf("Output_Figs/GO/invivo_Ast_Y_vs_O_GO_plot.pdf", width=15, height=10)
p <- ggplot(data=Ast, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-5,5,1)) + geom_bar(stat="identity", fill=(Ast$Color_Code), width=0.85) + theme + coord_flip() +labs(title="Ast GO", y="", x="")
dev.off()

