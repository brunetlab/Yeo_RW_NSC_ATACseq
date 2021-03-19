#!/usr/bin/Rscript

#Robin Yeo
#04/27/2020

#This script uses Diffbind to generate a consensus peakset for the in vitro Y+O Q+A libraries.
#It then generates a raw count matrix (no normalization) and exports it for downstream analysis.

#This script then combines all of the Y+O Q+A libraries from both in vivo and in vitro experiments and generates a consensus peakset and count matrix for downstream PCA+heatmaps.

rm(list=ls())
library(data.table)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/Dropbox/RWY_ATAC_Code_Checking/")


#Load in all of the replicates using a tab-delemited .csv file that contains information about file paths.

allPeaks <- dba(sampleSheet="IN_VITRO/Original_Data_InVitro/Y_O_qNSC_aNSC_InVitro_peaks.csv") 

#Prints data summary of allpeaks: number of peaks in each peakset, as well as the total number of unique peaks after merging overlapping ones
print(allPeaks) 

# 10 Samples, 121497 sites in matrix:
#   ID Tissue Factor Condition Replicate Caller Intervals
# 1  Y_aNSC_2   aNSC Y_aNSC     Young         2    bed    142735
# 2  Y_aNSC_3   aNSC Y_aNSC     Young         3    bed    142735
# 3  Y_aNSC_4   aNSC Y_aNSC     Young         4    bed    142735
# 4  O_aNSC_1   aNSC O_aNSC       Old         1    bed    150677
# 5  O_aNSC_2   aNSC O_aNSC       Old         2    bed    150677
# 6  O_aNSC_4   aNSC O_aNSC       Old         4    bed    150677
# 7  Y_qNSC_1   qNSC Y_qNSC     Young         1    bed     95078
# 8  Y_qNSC_2   qNSC Y_qNSC     Young         2    bed     95078
# 9  O_qNSC_2   qNSC O_qNSC       Old         2    bed     87458
# 10 O_qNSC_3   qNSC O_qNSC       Old         3    bed     87458

#Calculate a binding matrix with scores based on read counts for every sample
allPeaks.count.Raw = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

#This converts allPeaks.count.Raw to a dataframe for writing out tables
countMatrix.Raw <- dba.peakset(allPeaks.count.Raw, bRetrieve=TRUE, DataType = DBA_DATA_FRAME)

#Exporting the count matrix to a txt file for downstream use
write.table(countMatrix.Raw, file="IN_VITRO/Output_Data_InVitro/invitro_consensus_countMatrix_Raw.txt", sep="\t", row.names=FALSE)

#Exporting only the peaks (not sample count matrix values)
write.table(countMatrix.Raw[,c(1:3)], file="IN_VITRO/Output_Data_InVitro/invitro_consensus_peakset.txt", sep="\t", row.names=FALSE)



##---------------------------------------------------------------------------------------
#Generating Young/Old qNSC/aNSC consensus count matrix for combined in vivo+in vitro experiments
##---------------------------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/Dropbox/RWY_ATAC_Code_Checking/")


#Load in all of the replicates using a tab-delemited .csv file that contains information about file paths.
allPeaks <- dba(sampleSheet="IN_VITRO/Original_Data_InVitro/Combined_InVitro_InVivo_peaks.csv") 

#Prints data summary of allpeaks: number of peaks in each peakset, as well as the total number of unique peaks after merging overlapping ones
print(allPeaks) 

# 21 Samples, 156963 sites in matrix:
#   ID  Tissue Factor Condition Replicate Caller Intervals
# 1   Y_qNSC_1_invivo  InVivo Y_qNSC     Young         1    bed     70705
# 2   Y_aNSC_1_invivo  InVivo Y_aNSC     Young         1    bed     74840
# 3   O_qNSC_1_invivo  InVivo O_qNSC       Old         1    bed     69525
# 4   O_aNSC_1_invivo  InVivo O_aNSC       Old         1    bed     88383
# 5   Y_qNSC_2_invivo  InVivo Y_qNSC     Young         2    bed     70705
# 6   O_qNSC_2_invivo  InVivo O_qNSC       Old         2    bed     69525
# 7   O_aNSC_2_invivo  InVivo O_aNSC       Old         2    bed     88383
# 8   Y_qNSC_3_invivo  InVivo Y_qNSC     Young         3    bed     70705
# 9   Y_aNSC_3_invivo  InVivo Y_aNSC     Young         3    bed     74840
# 10  O_qNSC_3_invivo  InVivo O_qNSC       Old         3    bed     69525
# 11  O_aNSC_3_invivo  InVivo O_aNSC       Old         3    bed     88383
# 12 Y_aNSC_2_invitro InVitro Y_aNSC     Young         2    bed    142735
# 13 Y_aNSC_3_invitro InVitro Y_aNSC     Young         3    bed    142735
# 14 Y_aNSC_4_invitro InVitro Y_aNSC     Young         4    bed    142735
# 15 O_aNSC_1_invitro InVitro O_aNSC       Old         1    bed    150677
# 16 O_aNSC_2_invitro InVitro O_aNSC       Old         2    bed    150677
# 17 O_aNSC_4_invitro InVitro O_aNSC       Old         4    bed    150677
# 18 Y_qNSC_1_invitro InVitro Y_qNSC     Young         1    bed     95078
# 19 Y_qNSC_2_invitro InVitro Y_qNSC     Young         2    bed     95078
# 20 O_qNSC_2_invitro InVitro O_qNSC       Old         2    bed     87458
# 21 O_qNSC_3_invitro InVitro O_qNSC       Old         3    bed     87458


#Calculate a binding matrix with scores based on read counts for every sample
allPeaks.count.Raw = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

#This converts allPeaks.count.Raw to a dataframe for writing out tables
countMatrix.Raw <- dba.peakset(allPeaks.count.Raw, bRetrieve=TRUE, DataType = DBA_DATA_FRAME)

#Exporting the count matrix to a txt file for downstream use
write.table(countMatrix.Raw, file="IN_VITRO/Output_Data_InVitro/Combined_InVitro_InVivo_countMatrix_Raw.txt", sep="\t", row.names=FALSE)

#Exporting only the peaks (not sample count matrix values)
write.table(countMatrix.Raw[,c(1:3)], file="IN_VITRO/Output_Data_InVitro/Combined_InVitro_InVivo_peakset.txt", sep="\t", row.names=FALSE)


