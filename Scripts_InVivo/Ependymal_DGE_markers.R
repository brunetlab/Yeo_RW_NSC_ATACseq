# DGE analysis for ependymal vs NSCs cell markers in old mice for Yeo et al. paper
# load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(Matrix.utils)
library(glmnet)
library(caret) 
library(dplyr)
library(ggplot2)
library(Metrics)
library(stringr)
library(DESeq2)

#setwd("/labs/abrunet1/Eric/SummerRotation-GAGs")

# Make stacked violin plots with top genes
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# Load seurat object, contains all data
svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")
#dbn_preds <- readRDS("../SvzClockV5/transfer/DBN2019/data/bootstrap_pseudocell_15.rds")


#=============== Ependymal cell markers in old mice ===============

# Take old (>18 months) subset
old_svz <- subset(x = svz, subset = Age >= 18)
#old_svz <- subset(svz,subset.name="Age",low.threshold=18)

# Find markers between ependymal and aNSC-NPC
markers <- FindMarkers(old_svz, ident.1 = "Ependymal", ident.2 = "aNSC_NPC" ,test.use="negbinom", assay="RNA")
head(markers)
write.csv(markers,"ependymal_vs_aNSC_NPC_markers.csv")

# Find markers between ependymal and Astrocyte-qNSC
markers2 <- FindMarkers(old_svz, ident.1 = "Ependymal", ident.2 = "Astrocyte_qNSC" ,test.use="negbinom", assay="RNA")
head(markers2)
write.csv(markers2,"ependymal_vs_Astrocyte_qNSC_markers.csv")

# Find markers between ependymal and ALL other celltypes
markers3 <- FindMarkers(old_svz, ident.1 ="Ependymal", ident.2 = NULL, test.use="negbinom", assay="RNA")
head(markers3)
write.csv(markers3,"ependymal_vs_ALL_markers.csv")

# Make some UMAP plots with top genes
features <- c("Rarres2","Clu")
FeaturePlot(old_svz, features = features,reduction = "umap", label=TRUE)

features <- c("Pltp","Ascc1")
FeaturePlot(old_svz, features = features,reduction = "umap", label=TRUE)

# Violin plots of top 5 markers
StackedVlnPlot(obj = old_svz, features = rownames(markers)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers2)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers3)[1:5])



#=============== Ependymal cell markers in young mice (visualized in old) ===============

# Take young (<4 months) subset
young_svz <- subset(x = svz, subset = Age <= 5)

# Find markers between ependymal and aNSC-NPC
markers <- FindMarkers(young_svz, ident.1 = "Ependymal", ident.2 = "aNSC_NPC" ,test.use="negbinom", assay="RNA")
head(markers)
write.csv(markers,"young_ependymal_vs_aNSC_NPC_markers.csv")

# Find markers between ependymal and Astrocyte-qNSC
markers2 <- FindMarkers(young_svz, ident.1 = "Ependymal", ident.2 = "Astrocyte_qNSC" ,test.use="negbinom", assay="RNA")
head(markers2)
write.csv(markers2,"young_ependymal_vs_Astrocyte_qNSC_markers.csv")

# Find markers between ependymal and ALL other celltypes
markers3 <- FindMarkers(young_svz, ident.1 ="Ependymal", ident.2 = NULL, test.use="negbinom", assay="RNA")
head(markers3)
write.csv(markers3,"young_ependymal_vs_ALL_markers.csv")

# Make some UMAP plots with top genes
features <- c("Chchd10","Ascc1")
FeaturePlot(young_svz, features = features,reduction = "umap", label=TRUE)
FeaturePlot(old_svz, features = features,reduction = "umap", label=TRUE)

# Violin plots of top 5 markers
StackedVlnPlot(obj = young_svz, features = rownames(markers)[1:5])
StackedVlnPlot(obj = young_svz, features = rownames(markers2)[1:5])
StackedVlnPlot(obj = young_svz, features = rownames(markers3)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers2)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers3)[1:5])


#=============== Ependymal cell markers in all mice (visualized in old) ===============

# Find markers between ependymal and aNSC-NPC
markers <- FindMarkers(svz, ident.1 = "Ependymal", ident.2 = "aNSC_NPC" ,test.use="negbinom", assay="RNA")
head(markers)
write.csv(markers,"all_ependymal_vs_aNSC_NPC_markers.csv")

# Find markers between ependymal and Astrocyte-qNSC
markers2 <- FindMarkers(svz, ident.1 = "Ependymal", ident.2 = "Astrocyte_qNSC" ,test.use="negbinom", assay="RNA")
head(markers2)
write.csv(markers2,"all_ependymal_vs_Astrocyte_qNSC_markers.csv")

# Find markers between ependymal and ALL other celltypes
markers3 <- FindMarkers(svz, ident.1 ="Ependymal", ident.2 = NULL, test.use="negbinom", assay="RNA")
head(markers3)
write.csv(markers3,"all_ependymal_vs_ALL_markers.csv")

# Make some UMAP plots with top genes
features <- c("Rarres2","Clu")
FeaturePlot(young_svz, features = features,reduction = "umap", label=TRUE)
FeaturePlot(old_svz, features = features,reduction = "umap", label=TRUE)

# Violin plots of top 5 markers
StackedVlnPlot(obj = young_svz, features = rownames(markers)[1:5])
StackedVlnPlot(obj = young_svz, features = rownames(markers2)[1:5])
StackedVlnPlot(obj = young_svz, features = rownames(markers3)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers2)[1:5])
StackedVlnPlot(obj = old_svz, features = rownames(markers3)[1:5])



#=============== Visualizing aNSC cell markers (GFAP, CD133, EFGR) ===============
aNSC_markers <- c("Gfap", "Prom1", "Egfr", "Mki67")
StackedVlnPlot(obj = young_svz, features = aNSC_markers)
StackedVlnPlot(obj = old_svz, features = aNSC_markers)

VlnPlot(obj = young_svz, features=aNSC_markers, idents=c("Astrocyte_qNSC", "Ependymal"))
VlnPlot(obj = old_svz, features=aNSC_markers, idents=c("Astrocyte_qNSC", "Ependymal"))

#=============== Reactive astrocyte markers ===============
# From Luo et al., 2006: S100\beta, GFAP <-- both shared with aNSCs
# From Luo et al., 2008: CD24, possibly \beta-catenin <-- shared with ependymal
reactive_astrocyte_markers <- c("S100b","Gfap","Cd24a","Ctnnb1")
StackedVlnPlot(obj = young_svz, features = reactive_astrocyte_markers)
StackedVlnPlot(obj = old_svz, features = reactive_astrocyte_markers)

# How does % of cells with all reactive markers change with age?
reactive_svz <- subset(svz, S100b > 0 &  Gfap > 0 & Cd24a > 0 & Ctnnb1 > 0)
print(dim(reactive_svz))
print(reactive_svz[[]]$Age)
print(reactive_svz[[]]$Celltype.LowRes)

reactive_cells <- rownames(reactive_svz@meta.data)

#### FeaturePlot of cells with all markers
#svz[["reactive_marker"]] <- ifelse(svz[["RNA"]]["S100b"]>0 & svz[["RNA"]]["Gfap"]>0 & svz[["RNA"]]["Cd24a"]>0 & svz[["RNA"]]["Ctnnb1"]>0, "reactive", "non-reactive")
DimPlot(svz, cells.highlight = reactive_cells, cols.highlight = "red", cols = "gray") +
  scale_color_manual(labels = c("Non-reactive", "Reactive"), values = c("grey", "red"))
ggsave("plots/umap_reactive_all_markers.pdf", width = 7.5, height = 6)