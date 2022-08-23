#04/08/2019

#Adapted script from:
# Fig3 Interferon Gamma Response Signature Violins
# Written by Matthew Buckley

#The purpose of this script is to visualize gene expression changes within Ben Dulken's scRNA-seq data based on Gene Ontology pathways (downloaded from http://www.informatics.jax.org)


#### PART 1: Compare GO pathways in Young vs. Old

rm(list=ls())
library(Seurat)
library(ggplot2)
library(ggpubr)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")

# Load Data
load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

filenames <- list.files("Original_Data/GO_lists")
filenames <- gsub("\\..*","",filenames)


for (filename in filenames) {
  GO <- read.table(file=paste0("Original_Data/GO_lists/",filename,".txt"), row.names=NULL, sep="\t", header=TRUE,fill=TRUE)
  gene_list_GO <- data.frame(lapply(unique(GO$Symbol), as.character), stringsAsFactors=FALSE)
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
  AGE_CELLTYPES <- c("y_Astrocytes_qNSCs", "o_Astrocytes_qNSCs", "y_aNSCs_NPCs", "o_aNSCs_NPCs", "y_Neuroblasts", "o_Neuroblasts")
  adhesion_data$age_celltype <- paste0(adhesion_data$age,"_",adhesion_data$celltype)
  adhesion_data$age_celltype <- factor(adhesion_data$age_celltype, levels=AGE_CELLTYPES, ordered=T)
  
  # Plot parameters
  ageColors <- c("orange", "red2")
  cellColors <- c("lightskyblue", "orange","plum")
  agecelltypeColors <- c("lightskyblue", "dodgerblue", "orange", "red2", "plum", "purple1")
  s=0.4
  a1=0.1
  
  p <- ggplot(data=adhesion_data[adhesion_data$celltype == "Astrocytes_qNSCs" | adhesion_data$celltype == "aNSCs_NPCs" | adhesion_data$celltype == "Neuroblasts",],
              aes(x=age, y=adhesion_response)) +
    geom_violin(aes(fill=age_celltype), trim=T,scale="width", draw_quantiles = c(.5)) +
    geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.35, h = 0)) +
    scale_fill_manual(values=agecelltypeColors) +
    facet_wrap(~celltype, scales = "free", nrow =1) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 7)) +
    theme(legend.position="none") +
    theme_bw() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    stat_compare_means(method="wilcox.test")
  p
  ggsave(paste0("Output_Figs/Violin_plots/",filename,"_violin_plot.pdf"),p, height=6, width=15)
}
