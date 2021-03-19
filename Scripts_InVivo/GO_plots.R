#!/usr/bin/Rscript

#01/24/2019

library(ggplot2)


rm(list=ls())
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")


theme<-theme(aspect.ratio =3/2, panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black",size=15),axis.text.y=element_text(colour="black",size=15),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


##Outputs a barplot with GO enrichments for young qNSc vs. aNSC peaks
GO <- read.csv("Original_Data/GO_data_summary_qNSC_vs_aNSC.csv")

pdf("Output_Figs/qNSC_vs_aNSC_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=GO, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-6,10,1)) + geom_bar(stat="identity", fill=(GO$Color_Code), width=0.85) + theme + coord_flip() +labs(title="qNSC vs. aNSC GO Enrichment", y="", x="")
dev.off()



##Outputs young vs. old barplots with GO enrichments for qNSc and aNSC peaks
GO <- read.csv("Original_Data/GO_data_summary_Young_vs_Old.csv")
qNSC <- GO[which(GO$Condition=="qNSC"),]
aNSC <- GO[which(GO$Condition=="aNSC"),]

pdf("Output_Figs/qNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=qNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-6,5,1)) + geom_bar(stat="identity", fill=(qNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="qNSC GO Enrichment", y="", x="")
dev.off()

pdf("Output_Figs/aNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=aNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-4,6,1)) + geom_bar(stat="identity", fill=(aNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="aNSC GO Enrichment", y="", x="")
dev.off()


##Outputs a barplot with GO enrichments for the top 1000 genes driving the positive and negative PC2 directions from Y_O_Q_A PCA
GO <- read.csv("Original_Data/GO_data_summary_Q_A_PC2.csv")

pdf("Output_Figs/PC2_Y_O_Q_A_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=GO, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-16,6,2)) + geom_bar(stat="identity", fill=(GO$Color_Code), width=0.85) + theme + coord_flip() +labs(title="PC2 GO Enrichment", y="", x="")
dev.off()


##Outputs a barplot with IPA enrichments for old aNSC Canonical Pathways
IPA <- read.csv("Original_Data/IPA_data_summary_Old_aNSC.csv")

pdf("Output_Figs/IPA_old_aNSC_canonical_pathways_plot.pdf", width=10, height=15)
p <- ggplot(data=IPA, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(0,5,1)) + geom_bar(stat="identity", fill=(IPA$Color_Code), width=0.85) + theme + coord_flip() +labs(title="Old aNSC IPA Enrichment", y="", x="")
dev.off()

