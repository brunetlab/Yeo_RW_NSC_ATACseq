#!/usr/bin/Rscript

#Robin W. Yeo
#07/30/21

#In this script, I generate PCAs from the consensus count matrix of young qNSC/aNSC ATAC-seq peaks that overlap publicly available H3K27ac and p300 ChIPseq peaks from cultured NSCs (Martynoga et al., 2013)
#The consensus count matrix used as input for PCA visualization here was generated with the script "Generate_Consensus_Count_Matrix_Diffbind_enhancer_annotation.R"

#----------------------------------------
#YOUNG/OLD QNSC/ANSC READS IN NSC ATAC PEAKS THAT OVERLAP WITH H3K27AC AND P300
#----------------------------------------

rm(list=ls())
library("DESeq2")
library("ggplot2")
library(ChIPseeker)
library("pheatmap")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
set.seed(42)

setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")
consensus_count_matrix <- read.table("RWY_Output_Data/Y_O_qNSC_aNSC_countMatrix_H3K27ac_p300_Raw_2022.txt", header=TRUE, sep="\t")

sampleTable<-  colnames(consensus_count_matrix[,c(4:14)])
sampleTable<- as.data.frame(sampleTable)
names(sampleTable)[1]<- "Sample"
sampleTable$Group <- c("Y_qNSC", "Y_aNSC", "O_qNSC","O_aNSC", "Y_qNSC", "O_qNSC","O_aNSC", "Y_qNSC", "Y_aNSC", "O_qNSC","O_aNSC")
str(sampleTable)

dds <- DESeqDataSetFromMatrix(countData=consensus_count_matrix[,c(4:14)], sampleTable, ~Group)

dds <- DESeq(dds)

vst <- varianceStabilizingTransformation(dds, blind =F)

vstMat <- assay(vst)

########################################################
##CONSENSUS

pca.vst <- prcomp(t(vstMat), scale = F)
summary(pca.vst)

#format for ggplot
pca.vst.df <- as.data.frame(pca.vst$x)
str(pca.vst.df)
pca.vst.df$Age <- c("Young","Young","Old","Old","Young","Old","Old","Young","Young","Old","Old")
pca.vst.df$Cell_Type <- c("qNSC","aNSC","qNSC","aNSC","qNSC","qNSC","aNSC","qNSC","aNSC","qNSC","aNSC")
pca.vst.df$ID <- c("Y_qNSC","Y_aNSC","O_qNSC","O_aNSC","Y_qNSC","O_qNSC","O_aNSC","Y_qNSC","Y_aNSC","O_qNSC","O_aNSC")

#Y_qNSC lightskyblue
#O_qNSC dodgerblue
#Y_aNSC orange
#O_aNSC red2

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text=element_text(colour="black", size=rel(2)),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
             legend.title = element_text(color = "black", size = 40),legend.text = element_text(color = "black", size = 30))

pdf("Output_Figs/Y_O_Q_A_H3K27ac_p300_PC1_PC2_2022.pdf",height=11.8,width=15)
ggplot(pca.vst.df, aes(PC1,PC2, color=pca.vst.df$ID)) +
  geom_point(size=12, shape=18)+
  scale_color_manual(values=c("lightskyblue","orange","dodgerblue","red2"),breaks=c("Y_qNSC","Y_aNSC","O_qNSC","O_aNSC"),labels=c("Young qNSC","Young aNSC","Old qNSC","Old aNSC"), name="")+
  ggtitle("PCA - qNSC and aNSCs") +
  xlab("PC1 (50.50%)")+
  ylab("PC2 (10.68%)")+
  theme
dev.off()

pdf("Output_Figs/Y_O_Q_A_H3K27ac_p300_PC1_PC3_2022.pdf",height=11.8,width=15)
ggplot(pca.vst.df, aes(PC1,PC3, color=pca.vst.df$ID)) +
  geom_point(size=12, shape=18)+
  scale_color_manual(values=c("lightskyblue","orange","dodgerblue","red2"),breaks=c("Y_qNSC","Y_aNSC","O_qNSC","O_aNSC"),labels=c("Young qNSC","Young aNSC","Old qNSC","Old aNSC"), name="")+
  ggtitle("PCA - qNSC and aNSCs") +
  xlab("PC1 (50.50%)")+
  ylab("PC3 (7.70%)")+
  theme
dev.off()

