#04/08/2019

#Adapted script from Matthew Buckley

#The purpose of this script is to see gene expression changes within Ben Dulken's scRNA-seq data based on cell-adhesion hits from the ATAC data

rm(list=ls())
library(Seurat)
library(ggplot2)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")

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
adhesion_data$gene_sum <- rowSums(adhesion_data)
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

p <- ggplot(data=dplyr::filter(adhesion_data, celltype == "aNSCs_NPCs"), aes(x=age, y=gene_sum)) +
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

ggsave("Output_Figs/aNSC_NPC_violinplot_GO_Cell_Adhesion_2020.pdf",p, height=6, width=6)


ageColors <- c("lightskyblue", "dodgerblue")

p <- ggplot(data=dplyr::filter(adhesion_data, celltype == "Astrocytes_qNSCs"), aes(x=age, y=gene_sum)) +
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

ggsave("Output_Figs/Ast_qNSC_violinplot_GO_Cell_Adhesion_2020.pdf",p, height=6, width=6)


######################################################################################################################
rm(list=ls())
library(Seurat)
library(ggplot2)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")

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
adhesion_data$gene_sum <- rowSums(adhesion_data)
adhesion_data <- cbind(meta, adhesion_data)

# Reorder factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
adhesion_data$celltype <- factor(adhesion_data$celltype,  levels=CELLS, ordered=T)
adhesion_data$age <- factor(adhesion_data$age, levels=c("y", "o"), ordered=T)

p <- ggplot(data=dplyr::filter(adhesion_data, celltype %in% c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts"), age == "y"), aes(x=age, y=gene_sum, fill=celltype)) +
  geom_violin(trim=T, draw_quantiles = c(.5)) +
  scale_fill_manual(values=c("dodgerblue", "red2", "purple4")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

ggsave("Output_Figs/Young_qNSC_aNSC_NB_violinplot_GO_Cell_Adhesion.pdf",p, height=6, width=8)

# R version 3.5.2 (2018-12-20)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_3.1.5                             ChIPseeker_1.18.0                        ReactomePA_1.26.0                       
# [4] clusterProfiler_3.10.1                   ggplot2_3.3.2                            org.Mm.eg.db_3.7.0                      
# [7] TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4 GenomicFeatures_1.34.8                   AnnotationDbi_1.44.0                    
# [10] Biobase_2.42.0                           GenomicRanges_1.34.0                     GenomeInfoDb_1.18.2                     
# [13] IRanges_2.16.0                           S4Vectors_0.20.1                         BiocGenerics_0.28.0                     
# [16] RColorBrewer_1.1-2                       pheatmap_1.0.12                          calibrate_1.7.7                         
# [19] MASS_7.3-51.6                           
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.1.8                         fastmatch_1.1-0                         plyr_1.8.6                             
# [4] igraph_1.2.5                            lazyeval_0.2.2                          splines_3.5.2                          
# [7] listenv_0.8.0                           BiocParallel_1.16.6                     gridBase_0.4-7                         
# [10] urltools_1.7.3                          digest_0.6.25                           htmltools_0.5.0                        
# [13] GOSemSim_2.8.0                          viridis_0.5.1                           GO.db_3.7.0                            
# [16] gdata_2.18.0                            magrittr_1.5                            checkmate_2.0.0                        
# [19] memoise_1.1.0                           cluster_2.1.0                           ROCR_1.0-7                             
# [22] globals_0.12.5                          Biostrings_2.50.2                       graphlayouts_0.7.0                     
# [25] matrixStats_0.56.0                      enrichplot_1.2.0                        prettyunits_1.1.1                      
# [28] colorspace_1.4-1                        blob_1.2.1                              rappdirs_0.3.1                         
# [31] ggrepel_0.8.2                           dplyr_1.0.0                             crayon_1.3.4                           
# [34] RCurl_1.98-1.2                          jsonlite_1.6.1                          TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
# [37] graph_1.60.0                            zoo_1.8-8                               survival_3.2-3                         
# [40] ape_5.3                                 glue_1.4.1                              polyclip_1.10-0                        
# [43] gtable_0.3.0                            zlibbioc_1.28.0                         XVector_0.22.0                         
# [46] leiden_0.3.3                            UpSetR_1.4.0                            DelayedArray_0.8.0                     
# [49] graphite_1.28.2                         future.apply_1.5.0                      scales_1.1.1                           
# [52] DOSE_3.8.2                              DBI_1.1.0                               Rcpp_1.0.4.6                           
# [55] plotrix_3.7-8                           viridisLite_0.3.0                       progress_1.2.2                         
# [58] reticulate_1.16                         gridGraphics_0.5-0                      rsvd_1.0.3                             
# [61] bit_1.1-15.2                            reactome.db_1.66.0                      europepmc_0.4                          
# [64] tsne_0.1-3                              htmlwidgets_1.5.1                       httr_1.4.1                             
# [67] fgsea_1.8.0                             gplots_3.0.3                            ellipsis_0.3.1                         
# [70] ica_1.0-2                               pkgconfig_2.0.3                         XML_3.99-0.3                           
# [73] farver_2.0.3                            uwot_0.1.8                              ggplotify_0.0.5                        
# [76] tidyselect_1.1.0                        rlang_0.4.6                             reshape2_1.4.4                         
# [79] munsell_0.5.0                           tools_3.5.2                             generics_0.0.2                         
# [82] RSQLite_2.2.0                           ggridges_0.5.2                          stringr_1.4.0                          
# [85] yaml_2.2.1                              bit64_0.9-7                             fitdistrplus_1.1-1                     
# [88] tidygraph_1.2.0                         caTools_1.17.1.2                        purrr_0.3.4                            
# [91] RANN_2.6.1                              ggraph_2.0.3                            pbapply_1.4-2                          
# [94] future_1.17.0                           nlme_3.1-145                            DO.db_2.9                              
# [97] xml2_1.3.2                              biomaRt_2.38.0                          compiler_3.5.2                         
# [100] rstudioapi_0.11                         png_0.1-7                               plotly_4.9.2.1                         
# [103] tibble_3.0.1                            tweenr_1.0.1                            stringi_1.4.6                          
# [106] lattice_0.20-41                         Matrix_1.2-18                           vctrs_0.3.1                            
# [109] pillar_1.4.4                            lifecycle_0.2.0                         BiocManager_1.30.10                    
# [112] lmtest_0.9-37                           triebeard_0.3.0                         RcppAnnoy_0.0.16                       
# [115] irlba_2.3.3                             data.table_1.12.8                       cowplot_1.0.0                          
# [118] bitops_1.0-6                            patchwork_1.0.1                         rtracklayer_1.42.2                     
# [121] qvalue_2.14.1                           R6_2.4.1                                KernSmooth_2.23-16                     
# [124] gridExtra_2.3                           codetools_0.2-16                        boot_1.3-25                            
# [127] gtools_3.8.2                            SummarizedExperiment_1.12.0             withr_2.2.0                            
# [130] sctransform_0.2.1                       GenomicAlignments_1.18.1                Rsamtools_1.34.1                       
# [133] GenomeInfoDbData_1.2.0                  hms_0.5.3                               grid_3.5.2                             
# [136] tidyr_1.1.0                             rvcheck_0.1.8                           Rtsne_0.15                             
# [139] ggforce_0.3.2    