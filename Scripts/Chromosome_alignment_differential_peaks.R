#08_01_2018
#Bioconductor version 3.6

#This script will visualize the differentially accessible peaks between young and old qNSCs and aNSCs respectively aligned along the mouse chromosomes.

rm(list=ls())
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

### VISUALIZE THE DIFFERENTIALLY ACCESSIBLE PEAKS ALONG CHROMOSOME TRACKS
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")


pdf("Output_Figs/qNSC_Young_Differential_covplot.pdf")
covplot(readPeakFile("Output_Data/Differential_Peaks/qNSC_Young_Differential_Sites_withFC_FDR_0.05.narrowPeak"))
dev.off()

pdf("Output_Figs/qNSC_Old_Differential_covplot.pdf")
covplot(readPeakFile("Output_Data/Differential_Peaks/qNSC_Old_Differential_Sites_withFC_FDR_0.05.narrowPeak"))
dev.off()

pdf("Output_Figs/aNSC_Young_Differential_covplot.pdf")
covplot(readPeakFile("Output_Data/Differential_Peaks/aNSC_Young_Differential_Sites_withFC_FDR_0.05.narrowPeak"))
dev.off()

pdf("Output_Figs/aNSC_Old_Differential_covplot.pdf")
covplot(readPeakFile("Output_Data/Differential_Peaks/aNSC_Old_Differential_Sites_withFC_FDR_0.05.narrowPeak"))
dev.off()
