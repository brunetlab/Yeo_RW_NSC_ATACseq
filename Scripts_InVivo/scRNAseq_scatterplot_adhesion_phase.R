#!/usr/bin/Rscript

#Robin W. Yeo
#07/19/2022

#This script uses scRNAseq data from young/old SVZ (Dulken et al.) to visualize the relationship between adhesion signature and cell cycle phase in young/old qNSCs, aNSCs, and neuroblasts.

rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Dropbox/RWY_ATAC_Code_Checking/")

# Load Data
load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

#This outputs violin plots for young+old qnsc/ast and ansc/npc expression levels for genes within the Cell Adhesion GO category
filename <- "GO_Cell_Adhesion_GO0007155"
GO <- read.table(file=paste0("Original_Data/",filename,".txt"), row.names=NULL, sep="\t", header=TRUE,fill=TRUE)
gene_list_GO <- data.frame(lapply(unique(GO$Gene.name), as.character), stringsAsFactors=FALSE)
genelist <- tolower(c(gene_list_GO))

# Subset data to signature genes
colnames(d) <- tolower(colnames(d))
adhesion_data <- d[, colnames(d) %in% genelist]
meta <- d[, colnames(d) %in% c("s.score", "g2m.score", "phase","age", "replicate", "celltype")]
adhesion_data$adhesion_response <- rowSums(adhesion_data)
adhesion_data <- cbind(meta, adhesion_data)

# Reorder factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
adhesion_data$celltype <- factor(adhesion_data$celltype,  levels=CELLS, ordered=T)
adhesion_data$age <- factor(adhesion_data$age, levels=c("y", "o"), ordered=T)

# Plot parameters
ageColors <- c("orange", "red2")
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"))

p <- ggplot(data=dplyr::filter(adhesion_data, celltype == "aNSCs_NPCs"), aes(x=phase, y=adhesion_response)) +
  geom_violin(aes(fill=age), trim=T,scale="width", draw_quantiles = c(.5)) +
  scale_fill_manual(values=ageColors) +
  facet_wrap(~celltype, scales = "free", nrow =2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(legend.position="none") +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  #stat_compare_means(method="wilcox.test")
p

ggsave("Output_Figs/aNSC_NPC_violinplot_GO_Cell_Adhesion_by_phase_all_ages.pdf",p, height=6, width=6)


ageColors <- c("lightskyblue", "dodgerblue")

p <- ggplot(data=dplyr::filter(adhesion_data, celltype == "Astrocytes_qNSCs"), aes(x=phase, y=adhesion_response)) +
  geom_violin(aes(fill=age), trim=T,scale="width", draw_quantiles = c(.5)) +
  scale_fill_manual(values=ageColors) +
  facet_wrap(~celltype, scales = "free", nrow =2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(legend.position="none") +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  #stat_compare_means(method="wilcox.test")
p

ggsave("Output_Figs/Ast_qNSC_violinplot_GO_Cell_Adhesion_by_phase_all_ages.pdf",p, height=6, width=6)


ageColors <- c("mediumorchid1", "purple4")

p <- ggplot(data=dplyr::filter(adhesion_data, celltype == "Neuroblasts"), aes(x=phase, y=adhesion_response)) +
  geom_violin(aes(fill=age), trim=T,scale="width", draw_quantiles = c(.5)) +
  scale_fill_manual(values=ageColors) +
  facet_wrap(~celltype, scales = "free", nrow =2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(legend.position="none") +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#stat_compare_means(method="wilcox.test")
p

ggsave("Output_Figs/Neuroblasts_violinplot_GO_Cell_Adhesion_by_phase_all_ages.pdf",p, height=6, width=6)



###############Testing for phase/adhesion coupling

qNSCs <- dplyr::filter(adhesion_data, celltype == "Astrocytes_qNSCs") #1597 qNSCs/Astrocytes
Y_qNSCs <- dplyr::filter(qNSCs, age == "y") #1117
O_qNSCs <- dplyr::filter(qNSCs, age == "o") #480

aNSCs <- dplyr::filter(adhesion_data, celltype == "aNSCs_NPCs") #866 aNSCs/NPCs
Y_aNSCs <- dplyr::filter(aNSCs, age == "y") #784
O_aNSCs <- dplyr::filter(aNSCs, age == "o") #82

Neuroblasts <- dplyr::filter(adhesion_data, celltype == "Neuroblasts") #1014 Neuroblasts
Y_neuro <- dplyr::filter(Neuroblasts, age == "y") #868
O_neuro <- dplyr::filter(Neuroblasts, age == "o") #146

png("Output_Figs/Y_qNSC_barplot_cell_cycle_phase.png", width = 360, height = 480)
barplot(prop.table(table(Y_qNSCs$phase)), main="Y_qNSC", col="lightskyblue",ylim=c(0,1), cex.axis=2.5, cex.names=2.5)
dev.off()

png("Output_Figs/O_qNSC_barplot_cell_cycle_phase.png", width = 360, height = 480)
barplot(prop.table(table(O_qNSCs$phase)), main="O_qNSC", col="dodgerblue",ylim=c(0,1), cex.axis=2.5, cex.names=2.5)
dev.off()

png("Output_Figs/Y_aNSC_barplot_cell_cycle_phase.png", width = 360, height = 480)
barplot(prop.table(table(Y_aNSCs$phase)), main="Y_aNSC", col="orange",ylim=c(0,1), cex.axis=2.5, cex.names=2.5)
dev.off()

png("Output_Figs/O_aNSC_barplot_cell_cycle_phase.png", width = 360, height = 480)
barplot(prop.table(table(O_aNSCs$phase)), main="O_aNSC", col="red2",ylim=c(0,1), cex.axis=2.5, cex.names=2.5)
dev.off()

png("Output_Figs/Y_neuro_barplot_cell_cycle_phase.png", width = 360, height = 480)
barplot(prop.table(table(Y_neuro$phase)), main="Y_neuro", col="mediumorchid1",ylim=c(0,1), cex.axis=2.5, cex.names=2.5)
dev.off()

png("Output_Figs/O_neuro_barplot_cell_cycle_phase.png", width = 360, height = 480)
barplot(prop.table(table(O_neuro$phase)), main="O_neuro", col="purple4",ylim=c(0,1), cex.axis=2.5, cex.names=2.5)
dev.off()


barplot(prop.table(table(Y_qNSCs$phase)), main="Y_qNSC", col="dodgerblue",ylim=c(0,1), cex.axis=2.5, cex.names=2.5)
barplot(prop.table(table(Y_aNSCs$phase)), main="Y_aNSC", col="red2",ylim=c(0,0.6), cex.axis=2.5, cex.names=2.5)
barplot(prop.table(table(Y_neuro$phase)), main="Y_neuro", col="purple4",ylim=c(0,0.6), cex.axis=2.5, cex.names=2.5)


young_NSC_lineage <- rbind(Y_qNSCs, Y_aNSCs, Y_neuro)
old_NSC_lineage <- rbind(O_qNSCs, O_aNSCs, O_neuro)
whole_NSC_lineage <- rbind(young_NSC_lineage, old_NSC_lineage)

g <- ggplot(data=young_NSC_lineage, aes(x=s.score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,250)
g <- g+ggplot(data=old_NSC_lineage, aes(x=s.score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,250)
ggsave("Output_Figs/scatterplot_GO_Cell_Adhesion_by_phase_all_ages.pdf",g, height=6, width=12)

g <- ggplot(data=young_NSC_lineage, aes(x=s.score, y=adhesion_response, color=celltype)) + geom_density_2d() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,250)
g <- g+ggplot(data=old_NSC_lineage, aes(x=s.score, y=adhesion_response, color=celltype)) + geom_density_2d() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,250)
ggsave("Output_Figs/contourplot_GO_Cell_Adhesion_by_phase_all_ages.pdf",g, height=6, width=12)


g <- ggplot(data=whole_NSC_lineage, aes(x=s.score, y=adhesion_response, color=celltype, shape=age)) + geom_point(size=2.5) + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,250)
ggsave("Output_Figs/scatterplot_GO_Cell_Adhesion_s_score_combined.pdf",g, height=6, width=12)






#Anne asked me to generate separate plots with young (downsampled) and old of each cell type on their own graph

Y_qNSCs_480 <- sample_n(Y_qNSCs, 480)
Y_aNSCs_82 <- sample_n(Y_aNSCs, 82)
Y_neuro_146 <- sample_n(Y_neuro, 146)


g<-ggplot(data=rbind(Y_qNSCs_480, O_qNSCs), aes(x=s.score, y=adhesion_response, color=celltype, shape=age)) + geom_point() + scale_color_manual(values=c("dodgerblue")) + scale_shape_manual(values=c(19,1))+ ylim(0,250) + xlim(-0.2,1)
ggsave("Output_Figs/scatterplot_GO_Cell_Adhesion_s_score_downampled_qNSCs.pdf",g, height=6, width=6)

g<-ggplot(data=rbind(Y_aNSCs_82, O_aNSCs), aes(x=s.score, y=adhesion_response, color=celltype, shape=age)) + geom_point() + scale_color_manual(values=c("red2")) + scale_shape_manual(values=c(19,1))+ ylim(0,250) + xlim(-0.2,1)
ggsave("Output_Figs/scatterplot_GO_Cell_Adhesion_s_score_downampled_aNSCs.pdf",g, height=6, width=6)

g<-ggplot(data=rbind(Y_neuro_146, O_neuro), aes(x=s.score, y=adhesion_response, color=celltype, shape=age)) + geom_point() + scale_color_manual(values=c("purple4")) + scale_shape_manual(values=c(19,1))+ ylim(0,250) + xlim(-0.2,1)
ggsave("Output_Figs/scatterplot_GO_Cell_Adhesion_s_score_downampled_neuro.pdf",g, height=6, width=6)


theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text=element_text(colour="black", size=rel(2)),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position = "none")

g <- ggplot(data=rbind(Y_qNSCs_480,Y_aNSCs_82,Y_neuro_146), aes(x=s.score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,250) + xlim(-0.2,1) + theme
ggsave("Output_Figs/scatterplot_GO_Cell_Adhesion_by_phase_all_ages_downsampled_young.pdf",g, height=6, width=18)

g <- ggplot(data=old_NSC_lineage, aes(x=s.score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,250) + xlim(-0.2,1) + theme
ggsave("Output_Figs/scatterplot_GO_Cell_Adhesion_by_phase_all_ages_downsampled_old.pdf",g, height=6, width=18)





######################################################################################################################
rm(list=ls())
library(Seurat)
library(ggplot2)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Dropbox/RWY_ATAC_Code_Checking/")

# Load Data
load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

#This outputs violin plots for young+old qnsc/ast and ansc/npc expression levels for genes within the Cell Adhesion GO category
filename <- "GO_Cell_Adhesion_GO0007155"
GO <- read.table(file=paste0("Original_Data/",filename,".txt"), row.names=NULL, sep="\t", header=TRUE,fill=TRUE)
gene_list_GO <- data.frame(lapply(unique(GO$Gene.name), as.character), stringsAsFactors=FALSE)
genelist <- tolower(c(gene_list_GO))

# Subset data to signature genes
colnames(d) <- tolower(colnames(d))
adhesion_data <- d[, colnames(d) %in% genelist]
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype")]
adhesion_data$adhesion_response <- rowSums(adhesion_data)
adhesion_data <- cbind(meta, adhesion_data)

# Reorder factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
adhesion_data$celltype <- factor(adhesion_data$celltype,  levels=CELLS, ordered=T)
adhesion_data$age <- factor(adhesion_data$age, levels=c("y", "o"), ordered=T)

p <- ggplot(data=dplyr::filter(adhesion_data, celltype %in% c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts"), age == "y"), aes(x=age, y=adhesion_response, fill=celltype)) +
  geom_violin(trim=T, draw_quantiles = c(.5)) +
  scale_fill_manual(values=c("dodgerblue", "red2", "purple4")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

ggsave("Output_Figs/Young_qNSC_aNSC_NB_violinplot_GO_Cell_Adhesion.pdf",p, height=6, width=8)

