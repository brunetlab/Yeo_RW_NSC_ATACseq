library(Seurat)
library(tidyverse)
library(ggthemes)
library(cowplot)

Celltypes <- c("Oligodendro", "Microglia", "Endothelial",
               "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

#======================================== Shared Genes

svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")
meta <- svz[[]] # Metadata. Age, Celltype, and Celltype.LowRes most important
harmony <- Reductions(svz, slot="umap_har")@cell.embeddings # 2D batch effect correction reduction
colnames(harmony) <- c("umap_har_1", "umap_har_2") 
umap <- Reductions(svz, slot="umap")@cell.embeddings # Non corrected UMAP reduction
d <- tbl_df(cbind(umap, harmony, meta))

DefaultAssay(svz) <- "RNA"
exps <- svz[["RNA"]]@data


#========================================
# Bootstrap Cell Chrono Celltype Specific Specific Genes
celltypeDF <- readRDS("../../01_clocks/05_visualize/data/models_upset_data.rds")
generalG <- rownames(celltypeDF[rowSums(celltypeDF) >  1, ])

# Other exploratory figures
genes <- c("Alcam", "Ctnnd2", "Lsamp", "Itgb8", "Ntm")


exps2 <- exps[rownames(exps) %in% genes, ]
exprs_t <- t(as.matrix(exps2))
meta_counts <- cbind(d, exprs_t)
mclong <- meta_counts %>% pivot_longer(c(36:length(colnames(meta_counts))), names_to = "gene", values_to = "expression")
mclong$Celltype.LowRes <- factor(mclong$Celltype.LowRes, levels = Celltypes)


mclong2 <- mclong %>%
  select(hash.ID, Age, orig.ident, Celltype.LowRes, gene, expression) %>%
  group_by(hash.ID, Age, orig.ident, Celltype.LowRes, gene) %>%
  filter(Celltype.LowRes %in% c("Astrocyte_qNSC")) %>%
  filter(gene %in% c("Alcam", "Ctnnd2", "Itgb8")) %>%
  summarize(Mean = mean(expression))


ggplot(mclong2, aes(x = Age, y = Mean, color = gene)) +
  facet_wrap(Celltype.LowRes~., scales= "free_y") +
  geom_point() +
  geom_smooth(method = "loess", aes(fill = gene)) +
  theme_cowplot() +
  scale_color_few() +
  scale_fill_few()

ggsave("plots/other_trajectories.pdf", width=3.53, height=3.53)

# Other exploratory figures

exps2 <- exps[rownames(exps) %in% genes, ]
exprs_t <- t(as.matrix(exps2))
meta_counts <- cbind(d, exprs_t)
mclong <- meta_counts %>% pivot_longer(c(36:length(colnames(meta_counts))), names_to = "gene", values_to = "expression")
mclong$Celltype.LowRes <- factor(mclong$Celltype.LowRes, levels = Celltypes)


mclong2 <- mclong %>%
  select(hash.ID, Age, orig.ident, Celltype.LowRes, gene, expression) %>%
  group_by(hash.ID, Age, orig.ident, Celltype.LowRes, gene) %>%
  filter(Celltype.LowRes %in% c("aNSC_NPC")) %>%
  filter(gene %in% c("Alcam",  "Lsamp",  "Ntm")) %>%
  summarize(Mean = mean(expression))


ggplot(mclong2, aes(x = Age, y = Mean, color = gene)) +
  facet_wrap(Celltype.LowRes~., scales= "free_y") +
  geom_point() +
  geom_smooth(method = "loess", aes(fill = gene)) +
  theme_cowplot() +
  scale_color_few() +
  scale_fill_few()

ggsave("plots/other_trajectories2.pdf", width=3.53, height=3.53)



