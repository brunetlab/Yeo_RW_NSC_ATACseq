#!/usr/bin/Rscript

#Robin Yeo
#01/31/2020

#This script uses Diffbind to generate a consensus peak matrix for Young/Old qNSC/aNSC filled with the accessibility of the pooled sub-sampled (to 30M) reads for each of Y_qNSC, O_qNSC, Y_aNSC, O_aNSC



source("https://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("devtools")
biocLite("BiocParallel")

#Load DiffBind library

library("sva")
library("devtools")
library("BiocParallel")
library("limma")

#For list2df (https://www.rdocumentation.org/packages/qdapTools/versions/1.3.3/topics/list2df)
if (!require("pacman")) install.packages("pacman")
pacman::p_load_gh("trinker/qdapTools")

library(devtools)
library(easyGgplot2)

library(GenomicFeatures)
library(GenomicRanges)

library(org.Mm.eg.db)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)





rm(list=ls())
library(DiffBind)
library(ChIPseeker)
library(data.table)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/Dropbox/RWY_ATAC_Code_Checking/")

#Load in Data based on a tab-delemited .csv file that contains information about paths to 
#final.bam files & peak files

allPeaks <- dba(sampleSheet="Original_Data/Pooled_subsampled_Y_O_Q_A.csv")

#Prints data summary of allpeaks: number of peaks in each peakset, as well as the total number of unique peaks after merging overlapping ones
print(allPeaks) 

# 4 Samples, 42312 sites in matrix (87796 total):
#   ID Tissue        Factor Condition Caller Intervals
# 1 Y_qNSC   qNSC  Non-Dividing     Young    bed     70705
# 2 Y_aNSC   aNSC Proliferative     Young    bed     74840
# 3 O_qNSC   qNSC  Non-Dividing       Old    bed     69525
# 4 O_aNSC   aNSC Proliferative       Old    bed     88383


# #Calculate a binding matrix with scores based on read counts for every sample (DBA_SCORE_READS)
#Note: this is the RAW COUNT MATRIX since I'm already uploading sub-sampled bam files (to 30M reads each)
allPeaks.count = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

Y_O_qNSC_aNSC_consensus <- as.data.frame(dba.peakset(allPeaks.count, bRetrieve=TRUE))
Y_O_qNSC_aNSC_consensus <- makeGRangesFromDataFrame(Y_O_qNSC_aNSC_consensus, keep.extra.columns = TRUE)
Y_O_qNSC_aNSC_consensus <- as.data.frame(annotatePeak(Y_O_qNSC_aNSC_consensus, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db"))
write.table(Y_O_qNSC_aNSC_consensus, file="Output_Data/Y_O_Q_A_countMatrix_Pooled_30M_annotated.txt", sep="\t", row.names=FALSE) #87,796 peaks


Y_O_qNSC_aNSC_consensus_3UTR <- Y_O_qNSC_aNSC_consensus[Y_O_qNSC_aNSC_consensus$annotation == "3' UTR", ] #1668 peaks
Y_O_qNSC_aNSC_consensus_distal <- Y_O_qNSC_aNSC_consensus[Y_O_qNSC_aNSC_consensus$annotation == "Distal Intergenic", ] #31,660 peaks
Y_O_qNSC_aNSC_consensus_exon <- Y_O_qNSC_aNSC_consensus[Y_O_qNSC_aNSC_consensus$annotation %like% "Exon", ] #3904 peaks
Y_O_qNSC_aNSC_consensus_intron <- Y_O_qNSC_aNSC_consensus[Y_O_qNSC_aNSC_consensus$annotation %like% "Intron", ] #28,571 peaks
Y_O_qNSC_aNSC_consensus_promoter <- Y_O_qNSC_aNSC_consensus[Y_O_qNSC_aNSC_consensus$annotation%like% "Promoter", ] #20,633 peaks

write.table(Y_O_qNSC_aNSC_consensus_3UTR, file="Output_Data/Y_O_qNSC_aNSC_countMatrix_Pooled_30M_Annotated_3UTR.txt", quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)
write.table(Y_O_qNSC_aNSC_consensus_distal, file="Output_Data/Y_O_qNSC_aNSC_countMatrix_Pooled_30M_Annotated_distal.txt", quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)
write.table(Y_O_qNSC_aNSC_consensus_exon, file="Output_Data/Y_O_qNSC_aNSC_countMatrix_Pooled_30M_Annotated_exon.txt", quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)
write.table(Y_O_qNSC_aNSC_consensus_intron, file="Output_Data/Y_O_qNSC_aNSC_countMatrix_Pooled_30M_Annotated_intron.txt", quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)
write.table(Y_O_qNSC_aNSC_consensus_promoter, file="Output_Data/Y_O_qNSC_aNSC_countMatrix_Pooled_30M_Annotated_promoter.txt", quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)

