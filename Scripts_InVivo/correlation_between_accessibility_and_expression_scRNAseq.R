#!/usr/bin/Rscript

#Robin W. Yeo
#07/19/2022

#In this script, I want to visualize the correlation between promoter accessibility and gene expression genome-wide in the 5 cell types of interest in the young/old SVZ.
#I labelled ATAC-seq peaks from young/old SVZ populations as promoters using ChIPSeeker and then assess the correlation of promoter accessibility (VST-normalized) with gene expression (from scRNA-seq data).

rm(list=ls())
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(Seurat)
library(ggplot2)
library(ggpubr)
library(data.table)
library("Hmisc")
library(pheatmap)

setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")

#This reads in a VST_normalized count matrix of 24,887 ATAC-seq peaks identified as "Promoter" by ChIPSeeker from all young/old samples (outputted by "Generate_Consensus_Count_Matrix_Diffbind.R")
promoter_peaks <- readPeakFile("Original_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_consensus_countMatrix_VST.txt", header=TRUE)

promoter_peaks.count <- as.data.frame(annotatePeak(promoter_peaks, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db"))
promoter_peaks.count$SYMBOL <- tolower(promoter_peaks.count$SYMBOL)
genelist <- tolower(promoter_peaks.count$SYMBOL)

load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

# Subset data to signature genes
colnames(d) <- tolower(colnames(d))
expr_data <- d[, colnames(d) %in% genelist]
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype")]
expr_data <- cbind(meta, expr_data)

# Reorder factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
expr_data$celltype <- factor(expr_data$celltype,  levels=CELLS, ordered=T)
expr_data$age <- factor(expr_data$age, levels=c("y", "o"), ordered=T)

#Calculating mean expression value of each gene across all young cells [Note: first 3 cols of expr_data are c("age","replicate","celltype")]
Y_all_cells <- expr_data[expr_data$age=="y",]
Y_all_cells <- as.data.frame(colMeans(Y_all_cells[,c(4:ncol(Y_all_cells))]))
Y_all_cells <- tibble::rownames_to_column(Y_all_cells, "SYMBOL")
names(Y_all_cells)[2] <- "Y_all_cells_sc"

#Calculating mean expression value of each gene across all old cells [Note: first 3 cols of expr_data are c("age","replicate","celltype")]
O_all_cells <- expr_data[expr_data$age=="o",]
O_all_cells <- as.data.frame(colMeans(O_all_cells[,c(4:ncol(O_all_cells))]))
O_all_cells <- tibble::rownames_to_column(O_all_cells, "SYMBOL")
names(O_all_cells)[2] <- "O_all_cells_sc"

#Initializing new dataframe containing 14,503 genes to fill in with average expression values by young/old cell type
expression_sc <- as.data.frame(Y_all_cells$SYMBOL)
names(expression_sc)[1] <- "SYMBOL"

for (cell_type in CELLS) {
  temp <- expr_data[expr_data$celltype==cell_type & expr_data$age=="y",]
  temp <- as.data.frame(colMeans(temp[,c(4:ncol(temp))])) 
  names(temp)[1] <- paste0("Y_",cell_type,"_sc")
  expression_sc <- cbind(expression_sc,temp)
}

for (cell_type in CELLS) {
  temp <- expr_data[expr_data$celltype==cell_type & expr_data$age=="o",]
  temp <- as.data.frame(colMeans(temp[,c(4:ncol(temp))]))
  names(temp)[1] <- paste0("O_",cell_type,"_sc")
  expression_sc <- cbind(expression_sc,temp)
}

expression_sc <- cbind(expression_sc,Y_all_cells[,2])
expression_sc <- cbind(expression_sc,O_all_cells[,2])

#Merging both dataframes based on gene IDs so that I can compare promoter accessibility to expression levels
merged_df <- merge(promoter_peaks.count,expression_sc,by="SYMBOL")
names(merged_df)
#[1] "SYMBOL"                "seqnames"              "start"                 "end"                   "width"                 "strand"                "Y_Endo_1"             
#[8] "Y_qNSC_1"              "Y_aNSC_1"              "Y_NPC_1"               "O_Endo_1"              "O_Ast_1"               "O_qNSC_1"              "O_aNSC_1"             
#[15] "O_NPC_1"               "Y_Endo_2"              "Y_Ast_2"               "Y_qNSC_2"              "Y_NPC_2"               "O_Endo_2"              "O_Ast_2"              
#[22] "O_qNSC_2"              "O_aNSC_2"              "Y_Ast_3"               "Y_qNSC_3"              "Y_aNSC_3"              "Y_NPC_3"               "O_Ast_3"              
#[29] "O_qNSC_3"              "O_aNSC_3"              "O_NPC_3"               "annotation"            "geneChr"               "geneStart"             "geneEnd"              
#[36] "geneLength"            "geneStrand"            "geneId"                "transcriptId"          "distanceToTSS"         "ENSEMBL"               "GENENAME"             
#[43] "Y_Astrocytes_qNSCs_sc" "Y_aNSCs_NPCs_sc"       "Y_Neuroblasts_sc"      "Y_Neurons_sc"          "Y_OPC_sc"              "Y_Oligodendrocytes_sc" "Y_Endothelial_sc"     
#[50] "Y_Mural_cells_sc"      "Y_Microglia_sc"        "Y_Macrophages_sc"      "Y_T_cells_sc"          "O_Astrocytes_qNSCs_sc" "O_aNSCs_NPCs_sc"       "O_Neuroblasts_sc"     
#[57] "O_Neurons_sc"          "O_OPC_sc"              "O_Oligodendrocytes_sc" "O_Endothelial_sc"      "O_Mural_cells_sc"      "O_Microglia_sc"        "O_Macrophages_sc"     
#[64] "O_T_cells_sc"          "Y_all_cells[, 2]"      "O_all_cells[, 2]"   

#Taking the mean of VST-normalized ATAC-seq accessibility values for each age/celltype
labels <- list("Y_Ast_ATAC", "Y_qNSC_ATAC", "Y_aNSC_ATAC", "Y_NPC_ATAC", "Y_Endo_ATAC", "O_Ast_ATAC", "O_qNSC_ATAC", "O_aNSC_ATAC", "O_NPC_ATAC", "O_Endo_ATAC")
merged_df$Y_Ast_ATAC <- rowMeans(cbind(merged_df$Y_Ast_2,merged_df$Y_Ast_3))
merged_df$Y_qNSC_ATAC <- rowMeans(cbind(merged_df$Y_qNSC_1,merged_df$Y_qNSC_2,merged_df$Y_qNSC_3))
merged_df$Y_aNSC_ATAC <- rowMeans(cbind(merged_df$Y_aNSC_1,merged_df$Y_aNSC_3))
merged_df$Y_NPC_ATAC <- rowMeans(cbind(merged_df$Y_NPC_1,merged_df$Y_NPC_2,merged_df$Y_NPC_3))
merged_df$Y_Endo_ATAC <- rowMeans(cbind(merged_df$Y_Endo_1,merged_df$Y_Endo_2))
merged_df$O_Ast_ATAC <- rowMeans(cbind(merged_df$O_Ast_1,merged_df$O_Ast_2,merged_df$O_Ast_3))
merged_df$O_qNSC_ATAC <- rowMeans(cbind(merged_df$O_qNSC_1,merged_df$O_qNSC_2,merged_df$O_qNSC_3))
merged_df$O_aNSC_ATAC <- rowMeans(cbind(merged_df$O_aNSC_1,merged_df$O_aNSC_2,merged_df$O_aNSC_3))
merged_df$O_NPC_ATAC <- rowMeans(cbind(merged_df$O_NPC_1,merged_df$O_NPC_3))
merged_df$O_Endo_ATAC <- rowMeans(cbind(merged_df$O_Endo_1,merged_df$O_Endo_2))


## --------------------------------------------------------------------------
#Plotting deciled boxplots to correlate promoter accessibility with gene expression using Dulken scRNA-seq for all 5 Y/O cell types sorted from the SVZ
## --------------------------------------------------------------------------
library(OneR)

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black",size=10),axis.text.y=element_text(colour="black",size=15),
             axis.title.x=element_text(colour="black",size=15,vjust=-1),axis.title.y=element_text(colour="black",size=15),
             axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

deciles <- c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")

#--------------------------------------------------------------------------
Y_Endo_bin <- bin(merged_df$Y_Endo_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=Y_Endo_bin,y=merged_df$Y_Endothelial_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="hotpink") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/Y_Endo_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

O_Endo_bin <- bin(merged_df$O_Endo_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=O_Endo_bin,y=merged_df$O_Endothelial_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="deeppink") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/O_Endo_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)
#--------------------------------------------------------------------------
Y_Ast_bin <- bin(merged_df$Y_Ast_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=Y_Ast_bin,y=merged_df$Y_Astrocytes_qNSCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="lightgreen") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/Y_Ast_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

O_Ast_bin <- bin(merged_df$O_Ast_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=O_Ast_bin,y=merged_df$O_Astrocytes_qNSCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="darkgreen") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/O_Ast_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)
#--------------------------------------------------------------------------
Y_qNSC_bin <- bin(merged_df$Y_qNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=Y_qNSC_bin,y=merged_df$Y_Astrocytes_qNSCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="lightskyblue") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/Y_qNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

O_qNSC_bin <- bin(merged_df$O_qNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=O_qNSC_bin,y=merged_df$O_Astrocytes_qNSCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="dodgerblue") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/O_qNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)
#--------------------------------------------------------------------------
Y_aNSC_bin <- bin(merged_df$Y_aNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=Y_aNSC_bin,y=merged_df$Y_aNSCs_NPCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="orange") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/Y_aNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

O_aNSC_bin <- bin(merged_df$O_aNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=O_aNSC_bin,y=merged_df$O_aNSCs_NPCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="red2") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) +
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/O_aNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)
#--------------------------------------------------------------------------
Y_NPC_bin <- bin(merged_df$Y_NPC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=Y_NPC_bin,y=merged_df$Y_aNSCs_NPCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="plum") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/Y_NPC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

O_NPC_bin <- bin(merged_df$O_NPC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=O_NPC_bin,y=merged_df$O_aNSCs_NPCs_sc))
gg + #geom_jitter(aes(color='blue'),alpha=0.3) +
  geom_boxplot(fill="purple1") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  coord_cartesian(ylim=c(0, 1)) + 
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/ATAC_scRNAseq_correlation/O_NPC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)
#--------------------------------------------------------------------------

# R version 4.1.0 (2021-05-18)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.6.6
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] OneR_2.2                                  org.Mm.eg.db_3.13.0                       Hmisc_4.7-0                               Formula_1.2-4                            
# [5] survival_3.3-1                            lattice_0.20-45                           data.table_1.14.2                         ggpubr_0.4.0                             
# [9] SeuratObject_4.0.4                        Seurat_4.1.0                              TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0 GenomicFeatures_1.44.2                   
# [13] AnnotationDbi_1.54.1                      pheatmap_1.0.12                           ChIPseeker_1.28.3                         ggplot2_3.3.5                            
# [17] DESeq2_1.32.0                             SummarizedExperiment_1.22.0               Biobase_2.52.0                            MatrixGenerics_1.4.3                     
# [21] matrixStats_0.62.0                        GenomicRanges_1.44.0                      GenomeInfoDb_1.28.4                       IRanges_2.26.0                           
# [25] S4Vectors_0.30.2                          BiocGenerics_0.40.0      