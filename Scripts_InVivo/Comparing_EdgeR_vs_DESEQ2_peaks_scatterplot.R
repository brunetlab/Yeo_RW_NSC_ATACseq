#!/usr/bin/Rscript

#Robin W. Yeo
#07/19/2022

#During submission of the original paper, significant changes were made to Diffbind changing the original scripts used to call differential ATAC-seq peaks.
#1) we use Diffbind v3 EdgeR (with v2 settings: allPeaks$config$design <- FALSE) to generate new differential peak sets for qNSCs and aNSCs and check their overlap with the original peaks from Diffbind v2 EdgeR
#2) we use Diffbind v3 DESEQ2 (with v2 settings: allPeaks$config$design <- FALSE) to generate alternative differential peak sets for qNSCs and aNSCs and visualize how the original differential peaks' FDRs correlate with the p-values output from DESEQ2

##### V3 DIFFBIND EDGER WITH V2 CONFIGURATIONS
rm(list=ls())
library(DiffBind)
library(ChIPseeker)
library(parallel)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

for (cell_type in c("qNSC","aNSC")) {
  setwd("/Users/robinyeo/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
  allPeaks <- dba(sampleSheet=paste("Original_Data/",cell_type,"_peaks.csv",sep=""))
  print(allPeaks) 
  
  allPeaks$config$design <- FALSE
  
  allPeaks = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, summits=FALSE, filter=0, minCount=1, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)
  allPeaks = dba.contrast(allPeaks, categories=DBA_CONDITION, minMembers=2)
  allPeaks.analyze = dba.analyze(allPeaks, bParallel=TRUE, method=DBA_EDGER)
  
  rep_all <- dba.report(allPeaks.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
  Report <- annotatePeak(makeGRangesFromDataFrame(rep_all,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  
  Young <- rep_all[rep_all$Fold>0 & rep_all$FDR<0.05, ]
  Old <- rep_all[rep_all$Fold<0 & rep_all$FDR<0.05, ]
  
  Young <- rep_all[rep_all$Fold>0 & rep_all$'p-value'<0.01, ]
  Old <- rep_all[rep_all$Fold<0 & rep_all$'p-value'<0.01, ]
  
  write.table(as.data.frame(Report),file=paste("Output_Data/Diffbind_v3_comparison/", cell_type, "_Diffbind_EDGER_v3_v2_settings.csv",sep=""), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

##### V3 DIFFBIND DESEQ2 WITH V2 CONFIGURATIONS
rm(list=ls())
library(DiffBind)
library(ChIPseeker)
library(parallel)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


for (cell_type in c("qNSC","aNSC")) {
  setwd("/Users/robinyeo/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
  allPeaks <- dba(sampleSheet=paste("Original_Data/",cell_type,"_peaks.csv",sep=""))
  print(allPeaks) 
  
  allPeaks$config$design <- FALSE
  
  allPeaks = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, summits=FALSE, filter=0, minCount=1, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)
  allPeaks = dba.contrast(allPeaks, categories=DBA_CONDITION, minMembers=2)
  allPeaks.analyze = dba.analyze(allPeaks, bParallel=TRUE, method=DBA_DESEQ2)
  
  rep_all <- dba.report(allPeaks.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_DESEQ2, DataType=DBA_DATA_FRAME)
  Report <- annotatePeak(makeGRangesFromDataFrame(rep_all,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  
  Young <- rep_all[rep_all$Fold>0 & rep_all$FDR<0.05, ]
  Old <- rep_all[rep_all$Fold<0 & rep_all$FDR<0.05, ]
  
  Young <- rep_all[rep_all$Fold>0 & rep_all$'p-value'<0.01, ]
  Old <- rep_all[rep_all$Fold<0 & rep_all$'p-value'<0.01, ]
  
  write.table(as.data.frame(Report),file=paste("Output_Data/Diffbind_v3_comparison/", cell_type, "_Diffbind_DESEQ2_v3_v2_settings.csv",sep=""), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
}



##### COMPARING Y_QNSC ORIGNAL DAPs TO V3 EDGER & V3 DESEQ2 DAPS (WITH V2 SETTINGS)
rm(list=ls())
library(DiffBind)
library(ChIPseeker)
library(parallel)
library(ggplot2)
library(ggpubr)
library("VennDiagram")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#Before loading the new Diffbind v3 outputs, I had to manually trim columns 30+ (with Excel)
orig_EDGER_qNSC <- read.table("Original_Data/Differential_Peaks/qNSC_Young_Differential_Sites_withFC_FDR_0.05.narrowPeak", header=TRUE, sep="\t", quote="")
v3_EDGER_qNSC <- read.csv("Output_Data/Diffbind_v3_comparison/qNSC_Diffbind_EDGER_v3_v2_settings.csv")
v3_DESEQ2_qNSC <- read.csv("Output_Data/Diffbind_v3_comparison/qNSC_Diffbind_DESEQ2_v3_v2_settings.csv")

names(orig_EDGER_qNSC)[10:11] <- c("p.value.orig", "FDR.orig")
names(v3_EDGER_qNSC)[10:11] <- c("p.value.edgerv3", "FDR.edgerv3")
names(v3_DESEQ2_qNSC)[10:11] <- c("p.value.deseq2v3", "FDR.deseq2v3")

orig_EDGER_qNSC_Y <- orig_EDGER_qNSC[orig_EDGER_qNSC$Fold>0,c("seqnames","start","end","FDR.orig")]
v3_EDGER_qNSC_Y <- v3_EDGER_qNSC[v3_EDGER_qNSC$Fold>0,c("seqnames","start","end","FDR.edgerv3")]
v3_DESEQ2_qNSC_Y <- v3_DESEQ2_qNSC[v3_DESEQ2_qNSC$Fold>0,c("seqnames","start","end","FDR.deseq2v3")]

orig_EDGER_qNSC_Y_DAPs <- orig_EDGER_qNSC_Y[orig_EDGER_qNSC_Y$FDR.orig<=0.05,] #6670
v3_EDGER_qNSC_Y_DAPs <- v3_EDGER_qNSC_Y[v3_EDGER_qNSC_Y$FDR.edgerv3<=0.05,] #7221

#intersection between original peaks and v3 diffbind (w/ v2 settings) peaks
comparison <- merge(orig_EDGER_qNSC_Y_DAPs,v3_EDGER_qNSC_Y_DAPs, by=c("seqnames","start","end")) #5796
grid.newpage()
draw.pairwise.venn(area1=6670, area2=7221,cross.area=5796,category=c("Orig EdgeR","V3 EdgeR"),fill=c("Blue","Lightblue"))
dev.print(device=jpeg,filename="Output_Figs/venn_orig_v3_edger_Y_qNSC_peaks.png",width=480,height=300)

theme <- theme(legend.position="none", panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.text=element_text(colour="black", size=rel(2)))

#looking at original DAPs and what their FDRs look like in the v3 EdgeR method
comparison <- merge(orig_EDGER_qNSC_Y_DAPs,v3_EDGER_qNSC, by=c("seqnames","start","end")) #6583
ggplot(comparison, aes(x=-log(FDR.orig),y=-log(FDR.edgerv3))) + geom_point() + theme + xlab("-log(FDR) original peaks") + ylab("-log(FDR) v3 EdgeR peaks") + ggtitle("Comparison of Young qNSC FDRs between original and v3 EdgeR sets") + geom_smooth() + stat_cor()
ggsave("Output_Figs/correlation_orig_v3_edger_Y_qNSC_peaks.png")

#looking at original DAPs and what their FDRs look like in the v3 DESEQ2 method
comparison <- merge(orig_EDGER_qNSC_Y_DAPs,v3_DESEQ2_qNSC, by=c("seqnames","start","end")) #6583
ggplot(comparison, aes(x=-log(FDR.orig),y=-log(p.value.deseq2v3))) + geom_point() + theme + xlab("-log(FDR) original peaks") + ylab("-log(p.value) v3 DESEQ2 peaks") + ggtitle("Comparison of Young qNSC FDRs between original and v3 DESEQ2 sets") + geom_smooth() + stat_cor()
ggsave("Output_Figs/correlation_orig_v3_deseq2_Y_qNSC_peaks.png")


##### COMPARING O_ANSC ORIGNAL DAPs TO V3 EDGER & V3 DESEQ2 DAPS (WITH V2 SETTINGS)
rm(list=ls())
library(DiffBind)
library(ChIPseeker)
library(parallel)
library(ggplot2)
library(ggpubr)
library("VennDiagram")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#Before loading the new Diffbind v3 outputs, I had to manually trim columns 30+ (with Excel)
orig_EDGER_aNSC <- read.table("Original_Data/Differential_Peaks/aNSC_Old_Differential_Sites_withFC_FDR_0.05.narrowPeak", header=TRUE, sep="\t", quote="")
v3_EDGER_aNSC <- read.csv("Output_Data/Diffbind_v3_comparison/aNSC_Diffbind_EDGER_v3_v2_settings.csv")
v3_DESEQ2_aNSC <- read.csv("Output_Data/Diffbind_v3_comparison/aNSC_Diffbind_DESEQ2_v3_v2_settings.csv")

names(orig_EDGER_aNSC)[10:11] <- c("p.value.orig", "FDR.orig")
names(v3_EDGER_aNSC)[10:11] <- c("p.value.edgerv3", "FDR.edgerv3")
names(v3_DESEQ2_aNSC)[10:11] <- c("p.value.deseq2v3", "FDR.deseq2v3")

orig_EDGER_aNSC_O <- orig_EDGER_aNSC[orig_EDGER_aNSC$Fold<0,c("seqnames","start","end","FDR.orig")]
v3_EDGER_aNSC_O <- v3_EDGER_aNSC[v3_EDGER_aNSC$Fold<0,c("seqnames","start","end","FDR.edgerv3")]
v3_DESEQ2_aNSC_O <- v3_DESEQ2_aNSC[v3_DESEQ2_aNSC$Fold<0,c("seqnames","start","end","FDR.deseq2v3")]

orig_EDGER_aNSC_O_DAPs <- orig_EDGER_aNSC_O[orig_EDGER_aNSC_O$FDR.orig<=0.05,] #2013
v3_EDGER_aNSC_O_DAPs <- v3_EDGER_aNSC_O[v3_EDGER_aNSC_O$FDR.edgerv3<=0.05,] #3112

#intersection between original peaks and v3 diffbind (w/ v2 settings) peaks
comparison <- merge(orig_EDGER_aNSC_O_DAPs,v3_EDGER_aNSC_O_DAPs, by=c("seqnames","start","end")) #1850
grid.newpage()
draw.pairwise.venn(area1=2013, area2=3112,cross.area=1850,category=c("Orig EdgeR","V3 EdgeR"),fill=c("Red","Orange"))
dev.print(device=jpeg,filename="Output_Figs/venn_orig_v3_edger_O_aNSC_peaks.png",width=480,height=300)

theme <- theme(legend.position="none", panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.text=element_text(colour="black", size=rel(2)))

#looking at original DAPs and what their FDRs look like in the v3 EdgeR method
comparison <- merge(orig_EDGER_aNSC_O_DAPs,v3_EDGER_aNSC, by=c("seqnames","start","end")) #1984
ggplot(comparison, aes(x=-log(FDR.orig),y=-log(FDR.edgerv3))) + geom_point() + theme + xlab("-log(FDR) original peaks") + ylab("-log(FDR) v3 EdgeR peaks") + ggtitle("Comparison of Young qNSC FDRs between original and v3 EdgeR sets") + geom_smooth() + stat_cor()
ggsave("Output_Figs/correlation_orig_v3_edger_O_aNSC_peaks.png")

#looking at original DAPs and what their FDRs look like in the v3 DESEQ2 method
comparison <- merge(orig_EDGER_aNSC_O_DAPs,v3_DESEQ2_aNSC, by=c("seqnames","start","end")) #1984
ggplot(comparison, aes(x=-log(FDR.orig),y=-log(p.value.deseq2v3))) + geom_point() + theme + xlab("-log(FDR) original peaks") + ylab("-log(p.value) v3 DESEQ2 peaks") + ggtitle("Comparison of Young qNSC FDRs between original and v3 DESEQ2 sets") + geom_smooth() + stat_cor()
ggsave("Output_Figs/correlation_orig_v3_deseq2_O_aNSC_peaks.png")





