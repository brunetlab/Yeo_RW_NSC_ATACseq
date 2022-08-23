#!/usr/bin/Rscript

#Robin W. Yeo
#07/30/21

#The purpose of this script is to generate a count matrix of Young/Old qNSC/aNSC ATAC-seq peaks that overlap with H3K27ac and p300 ChIP-seq marks from Martynoga et al., 2013
#The overlapping peaks used as input here were generated in the scripts "intersect_H3K27ac_p300_ChIP_ATAC_Martynoga_2013.sh"

##---------------------------------------------------------------------------------------
#Generating Young/Old qNSC/aNSC consensus count matrix of ATAC/H3k27ac/p300 peaks
##---------------------------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")

allPeaks <- dba(sampleSheet="Original_Data/Y_O_qNSC_aNSC_peaks_H3K27ac_p300.csv")

#Prints data summary of allpeaks: number of peaks in each peakset, as well as the total number of unique peaks after merging overlapping ones
print(allPeaks)

# 11 Samples, 5640 sites in matrix:
#   ID Tissue Factor Condition Replicate Intervals
# 1  Y_qNSC_1   qNSC Y_qNSC     Young         1      3728
# 2  Y_aNSC_1   aNSC Y_aNSC     Young         1      2758
# 3  O_qNSC_1   qNSC O_qNSC       Old         1      3728
# 4  O_aNSC_1   aNSC O_aNSC       Old         1      2758
# 5  Y_qNSC_2   qNSC Y_qNSC     Young         2      3728
# 6  O_qNSC_2   qNSC O_qNSC       Old         2      3728
# 7  O_aNSC_2   aNSC O_aNSC       Old         2      2758
# 8  Y_qNSC_3   qNSC Y_qNSC     Young         3      3728
# 9  Y_aNSC_3   aNSC Y_aNSC     Young         3      2758
# 10 O_qNSC_3   qNSC O_qNSC       Old         3      3728
# 11 O_aNSC_3   aNSC O_aNSC       Old         3      2758

#Calculate a binding matrix with scores based on read counts for every sample
allPeaks.count.Raw = dba.count(allPeaks, summits=FALSE, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

#This converts allPeaks.count.Raw to a dataframe for writing out tables
countMatrix.Raw <- dba.peakset(allPeaks.count.Raw, bRetrieve=TRUE, DataType = DBA_DATA_FRAME)

#Exporting the count matrix to a txt file for downstream use
write.table(countMatrix.Raw, file="Output_Data/Y_O_qNSC_aNSC_countMatrix_H3K27ac_p300_Raw_2022.txt", sep="\t", row.names=FALSE)

#Exporting only the peaks (not sample count matrix values)
write.table(countMatrix.Raw[,c(1:3)], file="Output_Data/Y_O_qNSC_aNSC_peakset_H3K27ac_p300_2022.txt", sep="\t", row.names=FALSE)

