#!/usr/bin/Rscript

#Robin W. Yeo
#07/19/2022

#This script generates heatmaps of chromatin accessibility values at differentially accessible peaks within cell adhesion/migration pathways that change with age in qNSCs and aNSCs respectively.

library(pheatmap)
library(RColorBrewer)

##########################################################################################
#Plotting heatmap for differentially accessible peaks within GO Cell Adhesion category that change with age
##########################################################################################

rm(list=ls())
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
res <- read.table("Original_Data/Differential_Peaks/All_Diff_peaks/qNSC_Differential_Sites_withFC.narrowPeak", header=TRUE,quote="",sep="\t")
res_diff <- res[res$FDR<=0.05,]

#To subset to only differential peaks in cell adhesion pathway
GO <- read.table(file="Original_Data/GO_Cell_Adhesion_GO0007155.txt",sep="\t",header=TRUE)
res_diff_GO <- res_diff[res_diff$SYMBOL %in% GO$Gene.name,]
write.table(res_diff_GO, file="Output_Data/Differential_ATAC_peaks_adhesion_GO_qNSC_Y_vs_O.txt",sep="\t", row.names=F,quote=FALSE)
res_diff_GO <- res_diff_GO[,c("Y_qNSC_1","Y_qNSC_2","Y_qNSC_3","O_qNSC_1","O_qNSC_2", "O_qNSC_3")]

breaksList <- seq (-2,2,by=0.04)

pdf(paste("Output_Figs/Adhesion_heatmaps/qNSC_Cell_Adhesion.pdf"))
pheatmap(res_diff_GO,
         scale = "row",
         display_numbers = FALSE,
         #color = colorRampPalette(brewer.pal(n = 6, name = "YlOrRd"))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         treeheight_row = 0,
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaksList)),
         breaks=breaksList,
         show_rownames = FALSE,
         show_colnames = FALSE,
         border_color = NA)

dev.off()

rm(list=ls()) 
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
res <- read.table("Original_Data/Differential_Peaks/All_Diff_peaks/aNSC_Differential_Sites_withFC.narrowPeak", header=TRUE,quote="",sep="\t")
res_diff <- res[res$FDR<=0.05,]

#To subset to only differential peaks in cell adhesion pathway
GO <- read.table(file="Original_Data/GO_Cell_Adhesion_GO0007155.txt",sep="\t",header=TRUE)
res_diff_GO <- res_diff[res_diff$SYMBOL %in% GO$Gene.name,]
write.table(res_diff_GO, file="Output_Data/Differential_ATAC_peaks_adhesion_GO_aNSC_Y_vs_O.txt",sep="\t", row.names=F,quote=FALSE)
res_diff_GO <- res_diff_GO[,c("Y_aNSC_1","Y_aNSC_3","O_aNSC_1","O_aNSC_2", "O_aNSC_3")]

breaksList <- seq (-2,2,by=0.04)

pdf(paste("Output_Figs/Adhesion_heatmaps/aNSC_Cell_Adhesion.pdf"))
pheatmap(res_diff_GO,
         scale = "row",
         display_numbers = FALSE,
         #color = colorRampPalette(brewer.pal(n = 6, name = "YlOrRd"))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         treeheight_row = 0,
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaksList)),
         breaks=breaksList,
         show_rownames = FALSE,
         show_colnames = FALSE,
         border_color = NA)

dev.off()



##########################################################################################
#Plotting heatmap for differentially accessible peaks within GO negative Regulation of Cell Migration category that change with age
##########################################################################################

rm(list=ls())
library(pheatmap)
library(RColorBrewer)
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
res <- read.table("Original_Data/Differential_Peaks/All_Diff_peaks/qNSC_Differential_Sites_withFC.narrowPeak", header=TRUE,quote="",sep="\t")
res_diff <- res[res$FDR<=0.05,]

#To subset to only differential peaks in negative regulation cell migration pathway
GO <- read.table(file="Original_Data/negative_regulation_cell_migration_0030336.txt",sep="\t",header=TRUE)
res_diff_GO <- res_diff[res_diff$SYMBOL %in% GO$Symbol,]
write.table(res_diff_GO, file="Output_Data/Differential_ATAC_peaks_neg_regulation_migration_GO_qNSC_Y_vs_O.txt",sep="\t", row.names=F,quote=FALSE)
res_diff_GO <- res_diff_GO[,c("Y_qNSC_1","Y_qNSC_2","Y_qNSC_3","O_qNSC_1","O_qNSC_2", "O_qNSC_3")]

breaksList <- seq (-2,2,by=0.04)

pdf(paste("Output_Figs/Adhesion_heatmaps/qNSC_Negative_Regulation_Cell_Migration.pdf"))
pheatmap(res_diff_GO,
         scale = "row",
         display_numbers = FALSE,
         #color = colorRampPalette(brewer.pal(n = 6, name = "YlOrRd"))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         treeheight_row = 0,
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaksList)),
         breaks=breaksList,
         show_rownames = FALSE,
         show_colnames = FALSE,
         border_color = NA)
dev.off()


rm(list=ls()) 
library(pheatmap)
library(RColorBrewer)
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
res <- read.table("Original_Data/Differential_Peaks/All_Diff_peaks/aNSC_Differential_Sites_withFC.narrowPeak", header=TRUE,quote="",sep="\t")
res_diff <- res[res$FDR<=0.05,]

#To subset to only differential peaks in negative regulation cell migration pathway
GO <- read.table(file="Original_Data/negative_regulation_cell_migration_0030336.txt",sep="\t",header=TRUE)
res_diff_GO <- res_diff[res_diff$SYMBOL %in% GO$Symbol,]
write.table(res_diff_GO, file="Output_Data/Differential_ATAC_peaks_neg_cell_migration_GO_aNSC_Y_vs_O.txt",sep="\t", row.names=F,quote=FALSE)

res_diff_GO <- res_diff_GO[,c("Y_aNSC_1","Y_aNSC_3","O_aNSC_1","O_aNSC_2", "O_aNSC_3")]

breaksList <- seq (-2,2,by=0.04)

pdf(paste("Output_Figs/Adhesion_heatmaps/aNSC_Negative_Regulation_Cell_Migration.pdf"))
pheatmap(res_diff_GO,
         scale = "row",
         display_numbers = FALSE,
         #color = colorRampPalette(brewer.pal(n = 6, name = "YlOrRd"))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         treeheight_row = 0,
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaksList)),
         breaks=breaksList,
         show_rownames = FALSE,
         show_colnames = FALSE,
         border_color = NA)
dev.off()


