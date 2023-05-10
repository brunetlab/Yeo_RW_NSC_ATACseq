library(Seurat)
library(ggplot2)
library(stringr)
#library(MAST)

# Load Data
load("data/Dulken_svz_celltypes_2019-01-31.rda")
svz <- UpdateSeuratObject(object = svz)
svz <- SetIdent(svz, value="Age")
sub_svz <- subset(x = svz, subset = Celltype %in% c("aNSCs_NPCs", "Astrocytes_qNSCs", "Neuroblasts"))

sub_svz <- RunPCA(sub_svz, verbose=F)
sub_svz <- FindNeighbors(sub_svz, dims = 1:21)
sub_svz <- RunUMAP(sub_svz, dims=1:21)

DimPlot(sub_svz, group.by="Age", pt.size = .6)
ggsave("plots/umap_nsc_lineage_dulken.pdf", width = 7.5, height = 6)
DimPlot(sub_svz, group.by="Celltype", pt.size = .6, cols=c("dodgerblue", "purple4", "red2"))
ggsave("plots/umap_nsc_lineage_dulken_celltype.pdf", width = 7.5, height = 6)