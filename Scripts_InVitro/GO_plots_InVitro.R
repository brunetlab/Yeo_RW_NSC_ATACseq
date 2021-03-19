#!/usr/bin/Rscript

#04/27/2020

library(ggplot2)


rm(list=ls())
setwd("~/Dropbox/RWY_ATAC_Code_Checking/")

theme<-theme(aspect.ratio =3/2, panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black",size=15),axis.text.y=element_text(colour="black",size=15),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


##Outputs a barplot with GO enrichments for in vitro young qNSc vs. aNSC peaks
GO <- read.csv("IN_VITRO/Output_Data_InVitro/GO_data_summary_qNSC_vs_aNSC_InVitro.csv")

pdf("IN_VITRO/Output_Figs_InVitro/qNSC_vs_aNSC_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=GO, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-8,4,1)) + geom_bar(stat="identity", fill=(GO$Color_Code), width=0.85) + theme + coord_flip() +labs(title="qNSC vs. aNSC GO Enrichment", y="", x="")
dev.off()



##Outputs young vs. old barplots with GO enrichments for qNSc and aNSC peaks
GO <- read.csv("IN_VITRO/Output_Data_InVitro/GO_data_summary_Young_vs_Old_InVitro.csv")
qNSC <- GO[which(GO$Condition=="qNSC"),]
aNSC <- GO[which(GO$Condition=="aNSC"),]

pdf("IN_VITRO/Output_Figs_InVitro/qNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=qNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-5,7,1)) + geom_bar(stat="identity", fill=(qNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="qNSC GO Enrichment", y="", x="")
dev.off()

pdf("IN_VITRO/Output_Figs_InVitro/aNSC_Y_vs_O_GO_plot.pdf", width=10, height=10)
p <- ggplot(data=aNSC, aes(x=reorder(Term,(Order)), y=log_p_val))
p + scale_y_continuous(breaks=seq(-5,4,1)) + geom_bar(stat="identity", fill=(aNSC$Color_Code), width=0.85) + theme + coord_flip() +labs(title="aNSC GO Enrichment", y="", x="")
dev.off()

