#!/usr/bin/Rscript

#Robin W. Yeo
#04/27/2020

#This script reads in *.csv file containing significance and fold-change values for GO terms (found in Extended Data Tables) and plots them using ggplot2.

rm(list=ls())
library(ggplot2)
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
theme<-theme(aspect.ratio =3/2, panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black",size=15),axis.text.y=element_text(colour="black",size=15),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


##Outputs young vs. old barplots with GO enrichments for qNSc and aNSC peaks
GO <- read.csv("Original_Data/GO_data_summary_Young_vs_Old_InVitro.csv")
qNSC <- GO[which(GO$Condition=="qNSC"),]
aNSC <- GO[which(GO$Condition=="aNSC"),]

pdf("Output_Figs/GO/invitro_qNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=qNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-5,7,1)) + geom_bar(stat="identity", fill=(qNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="qNSC GO Enrichment", y="", x="")
dev.off()

pdf("Output_Figs/GO/invitro_aNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=aNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-5,4,1)) + geom_bar(stat="identity", fill=(aNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="aNSC GO Enrichment", y="", x="")
dev.off()

