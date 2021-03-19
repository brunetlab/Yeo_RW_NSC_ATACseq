#!/usr/bin/Rscript

#Robin Yeo
#04/24/2020

#This script is to use Diffbind to generate lists of differentially accessible peaks that are annotated with nearby genes (for downstream GO/KEGG enrichment).
#The first part of the script computes differential peaks for young vs. old based on the Y+O consensus peakset
#The second part of the script computes the Q vs A differential peaks using the Y+O Q+A consensus peakset

rm(list=ls())
library(DiffBind)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#This for loop outputs a list of chromatin peaks that are differentially open in both the young and old conditions
#The differential peaks are used downstream for motif enrichment and the annotated genes are used downstream for GO/KEGG pathway enrichment.


cell_types <- c("qNSC", "aNSC")
for (cell_type in cell_types) {
threshold=0.05
setwd("~/Dropbox/RWY_ATAC_Code_Checking/")
allPeaks <- dba(sampleSheet=paste("IN_VITRO/Original_Data_InVitro/",cell_type,"_InVitro_peaks.csv",sep=""))
print(allPeaks) 

allPeaks.count = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)
uncorrected.contrast = dba.contrast(allPeaks.count, categories=DBA_CONDITION, minMembers=2)
allPeaks.contrast.analyze = dba.analyze(uncorrected.contrast, bCorPlot=FALSE, bParallel=TRUE, bTagwise=FALSE, bFullLibrarySize=TRUE, bReduceObjects=FALSE,method=DBA_EDGER) 

#Outputting the full list of peaks from the consensus peakset with associated FDRs and FCs from differential peak calling
rep_all <- dba.report(allPeaks.contrast.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report_all <- annotatePeak(makeGRangesFromDataFrame(rep_all,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(as.data.frame(Report_all),file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/All_Diff_peaks/",cell_type, "_Differential_Sites_withFC.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#Outputting the list of peaks that are differentially enriched in young or old conditions with FDR < 0.05
rep_thresh <- dba.report(allPeaks.contrast.analyze, contrast=1, th=threshold, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report <- annotatePeak(makeGRangesFromDataFrame(rep_thresh,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
Report <- as.data.frame(Report)
indexes = Report$Fold>0
Young_Enriched = Report[indexes,]
Old_Enriched = Report[!indexes,]

write.table(Young_Enriched, file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/",cell_type, "_Young_InVitro_Differential_Sites_withFC_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(Old_Enriched, file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/",cell_type, "_Old_InVitro_Differential_Sites_withFC_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#RWY's outputs from print(allPeaks)

# 4 Samples, 76324 sites in matrix:
#   ID Tissue       Factor Condition Replicate Caller Intervals
# 1 Y_qNSC_1   qNSC Non-Dividing     Young         1    bed     95078
# 2 Y_qNSC_2   qNSC Non-Dividing     Young         2    bed     95078
# 3 O_qNSC_2   qNSC Non-Dividing       Old         2    bed     87458
# 4 O_qNSC_3   qNSC Non-Dividing       Old         3    bed     87458
# 
# 6 Samples, 107967 sites in matrix:
#   ID Tissue        Factor Condition Replicate Caller Intervals
# 1 Y_aNSC_2   aNSC Proliferative     Young         2    bed    142735
# 2 Y_aNSC_3   aNSC Proliferative     Young         3    bed    142735
# 3 Y_aNSC_4   aNSC Proliferative     Young         4    bed    142735
# 4 O_aNSC_1   aNSC Proliferative       Old         1    bed    150677
# 5 O_aNSC_2   aNSC Proliferative       Old         2    bed    150677
# 6 O_aNSC_4   aNSC Proliferative       Old         4    bed    150677

###################################################################################################
### The following outputs differential peaks for qNSC vs. aNSC in either young or old conditions
###################################################################################################

############
#NOTE: THE THRESHOLD HERE HAS BEEN SET TO 0.0001 (INSTEAD OF 0.05 AS ABOVE AND IN THE IN VIVO ANALYSES)
############


rm(list=ls())
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
threshold=0.0001
setwd("~/Dropbox/RWY_ATAC_Code_Checking/")
allPeaks <- dba(sampleSheet="IN_VITRO/Original_Data_InVitro/Y_O_qNSC_aNSC_InVitro_peaks.csv")
print(allPeaks) 

# 10 Samples, 121497 sites in matrix:
#   ID Tissue        Factor Condition Replicate Caller Intervals
# 1  Y_aNSC_2   aNSC Proliferative     Young         2    bed    142735
# 2  Y_aNSC_3   aNSC Proliferative     Young         3    bed    142735
# 3  Y_aNSC_4   aNSC Proliferative     Young         4    bed    142735
# 4  O_aNSC_1   aNSC Proliferative       Old         1    bed    150677
# 5  O_aNSC_2   aNSC Proliferative       Old         2    bed    150677
# 6  O_aNSC_4   aNSC Proliferative       Old         4    bed    150677
# 7  Y_qNSC_1   qNSC  Non-Dividing     Young         1    bed     95078
# 8  Y_qNSC_2   qNSC  Non-Dividing     Young         2    bed     95078
# 9  O_qNSC_2   qNSC  Non-Dividing       Old         2    bed     87458
# 10 O_qNSC_3   qNSC  Non-Dividing       Old         3    bed     87458

allPeaks.count = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

#This calculates differential peaks during differentiation of YOUNG cell types
uncorrected.contrast = dba.contrast(allPeaks.count, group1=allPeaks.count$masks$Y_qNSC,group2=allPeaks.count$masks$Y_aNSC, name1="Y_qNSC",name2="Y_aNSC", minMembers=2)
Y.allPeaks.contrast.analyze = dba.analyze(uncorrected.contrast, bCorPlot=FALSE, bParallel=TRUE, bTagwise=FALSE, bFullLibrarySize=TRUE, bReduceObjects=FALSE,method=DBA_EDGER) 

rep_all_Y <- dba.report(Y.allPeaks.contrast.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report_all_Y <- annotatePeak(makeGRangesFromDataFrame(rep_all_Y,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(as.data.frame(Report_all_Y),file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/All_Diff_peaks/Young_qNSC_aNSC_Differential_Sites_withFC.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

rep_Y <- dba.report(Y.allPeaks.contrast.analyze, contrast=1, th=threshold, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report <- annotatePeak(makeGRangesFromDataFrame(rep_Y,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
Report <- as.data.frame(Report)

indexes = Report$Fold>0
Q_Enriched = Report[indexes,]
A_Enriched = Report[!indexes,]

write.table(Q_Enriched, file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/qNSC_Young_Q_A_InVitro_Differential_Report_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(A_Enriched, file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/aNSC_Young_Q_A_InVitro_Differential_Report_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


#This calculates differential peaks during differentiation of OLD cell types
uncorrected.contrast = dba.contrast(allPeaks.count, group1=allPeaks.count$masks$O_qNSC,group2=allPeaks.count$masks$O_aNSC, name1="O_qNSC",name2="O_aNSC", minMembers=2)
O.allPeaks.contrast.analyze = dba.analyze(uncorrected.contrast, bCorPlot=FALSE, bParallel=TRUE, bTagwise=FALSE, bFullLibrarySize=TRUE, bReduceObjects=FALSE,method=DBA_EDGER) 

rep_all_O <- dba.report(O.allPeaks.contrast.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report_all_O <- annotatePeak(makeGRangesFromDataFrame(rep_all_O,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(as.data.frame(Report_all_O),file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/All_Diff_peaks/Old_qNSC_aNSC_Differential_Sites_withFC.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

rep_O <- dba.report(O.allPeaks.contrast.analyze, contrast=1, th=threshold, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report <- annotatePeak(makeGRangesFromDataFrame(rep_O,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
Report <- as.data.frame(Report)

indexes = Report$Fold>0
Q_Enriched = Report[indexes,]
A_Enriched = Report[!indexes,]

write.table(Q_Enriched, file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/qNSC_Old_Q_A_InVitro_Differential_Report_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(A_Enriched, file=paste("IN_VITRO/Output_Data_InVitro/Differential_Peaks/aNSC_Old_Q_A_InVitro_Differential_Report_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
