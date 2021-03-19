# Figure 3A Marker tSNE
# Matthew Buckley

rm(list=ls())
library(Seurat)
library(ggplot2)
library(tidyverse)


# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")


# Load Data
load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata
features <- colnames(d)

GO <- read.delim(file="Original_Data/GO_Cell_Adhesion_GO0007155.txt", header=TRUE)
genelist <- tolower(GO$Gene.name)


colnames(d) <- tolower(colnames(d))
adhesion_data <- d[, colnames(d) %in% genelist]
adhesion_data$gene_sum <- rowSums(adhesion_data)
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype", "tsne_1", "tsne_2")]
d2 <- cbind(meta, adhesion_data)
d2$age <- factor(d2$age, levels = c("y", "o"), ordered=TRUE)

low_col <- "grey"
high_col <- "darkred"
size_range <- c(.8, 1.2)
alpha_range <- c(0.10, .90)

# Saturate color at 99.5% percentile of response
threshold <- quantile(d2$gene_sum, .99)
adh_sig_cap <- d2$gene_sum
adh_sig_cap[adh_sig_cap > threshold] <- threshold
d2$adh_sig_cap <- adh_sig_cap

# Plot IFN gamma response
q <- ggplot(data = filter(d2, age == "y"), aes(x = tsne_1, y = tsne_2,
			size = adh_sig_cap, alpha = adh_sig_cap, color = adh_sig_cap))
q <- q + geom_point() + scale_size(range = size_range) + scale_alpha(range = alpha_range)
q <- q + scale_colour_gradient(low = low_col, high = high_col)
q <- q + theme(legend.position="none", panel.background = element_blank(),axis.line = element_line(colour = "black"))
q

ggsave("Output_Figs/YOUNG_tSNE_GO_Cell_Adhesion_GO0007155.png", q, height = 7, width = 6)


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
#   [1] forcats_0.5.0                            stringr_1.4.0                            dplyr_1.0.0                             
# [4] purrr_0.3.4                              readr_1.3.1                              tidyr_1.1.0                             
# [7] tibble_3.0.1                             tidyverse_1.3.0                          Seurat_3.1.5                            
# [10] ChIPseeker_1.18.0                        ReactomePA_1.26.0                        clusterProfiler_3.10.1                  
# [13] ggplot2_3.3.2                            org.Mm.eg.db_3.7.0                       TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4
# [16] GenomicFeatures_1.34.8                   AnnotationDbi_1.44.0                     Biobase_2.42.0                          
# [19] GenomicRanges_1.34.0                     GenomeInfoDb_1.18.2                      IRanges_2.16.0                          
# [22] S4Vectors_0.20.1                         BiocGenerics_0.28.0                      RColorBrewer_1.1-2                      
# [25] pheatmap_1.0.12                          calibrate_1.7.7                          MASS_7.3-51.6                           
# 
# loaded via a namespace (and not attached):
#   [1] reticulate_1.16                         tidyselect_1.1.0                        RSQLite_2.2.0                          
# [4] htmlwidgets_1.5.1                       grid_3.5.2                              BiocParallel_1.16.6                    
# [7] Rtsne_0.15                              munsell_0.5.0                           codetools_0.2-16                       
# [10] ica_1.0-2                               future_1.17.0                           withr_2.2.0                            
# [13] colorspace_1.4-1                        GOSemSim_2.8.0                          rstudioapi_0.11                        
# [16] ROCR_1.0-7                              DOSE_3.8.2                              listenv_0.8.0                          
# [19] urltools_1.7.3                          GenomeInfoDbData_1.2.0                  polyclip_1.10-0                        
# [22] bit64_0.9-7                             farver_2.0.3                            vctrs_0.3.1                            
# [25] generics_0.0.2                          R6_2.4.1                                graphlayouts_0.7.0                     
# [28] rsvd_1.0.3                              bitops_1.0-6                            fgsea_1.8.0                            
# [31] gridGraphics_0.5-0                      DelayedArray_0.8.0                      assertthat_0.2.1                       
# [34] scales_1.1.1                            ggraph_2.0.3                            enrichplot_1.2.0                       
# [37] gtable_0.3.0                            globals_0.12.5                          tidygraph_1.2.0                        
# [40] rlang_0.4.6                             splines_3.5.2                           rtracklayer_1.42.2                     
# [43] lazyeval_0.2.2                          broom_0.5.6                             europepmc_0.4                          
# [46] checkmate_2.0.0                         modelr_0.1.8                            BiocManager_1.30.10                    
# [49] yaml_2.2.1                              reshape2_1.4.4                          backports_1.1.8                        
# [52] qvalue_2.14.1                           tools_3.5.2                             ggplotify_0.0.5                        
# [55] gridBase_0.4-7                          ellipsis_0.3.1                          gplots_3.0.3                           
# [58] ggridges_0.5.2                          Rcpp_1.0.4.6                            plyr_1.8.6                             
# [61] progress_1.2.2                          zlibbioc_1.28.0                         RCurl_1.98-1.2                         
# [64] prettyunits_1.1.1                       pbapply_1.4-2                           viridis_0.5.1                          
# [67] cowplot_1.0.0                           zoo_1.8-8                               haven_2.3.1                            
# [70] SummarizedExperiment_1.12.0             ggrepel_0.8.2                           cluster_2.1.0                          
# [73] fs_1.4.1                                magrittr_1.5                            data.table_1.12.8                      
# [76] DO.db_2.9                               reprex_0.3.0                            triebeard_0.3.0                        
# [79] lmtest_0.9-37                           RANN_2.6.1                              reactome.db_1.66.0                     
# [82] fitdistrplus_1.1-1                      matrixStats_0.56.0                      hms_0.5.3                              
# [85] patchwork_1.0.1                         XML_3.99-0.3                            readxl_1.3.1                           
# [88] gridExtra_2.3                           compiler_3.5.2                          biomaRt_2.38.0                         
# [91] KernSmooth_2.23-16                      crayon_1.3.4                            htmltools_0.5.0                        
# [94] lubridate_1.7.9                         DBI_1.1.0                               tweenr_1.0.1                           
# [97] dbplyr_1.4.4                            rappdirs_0.3.1                          boot_1.3-25                            
# [100] Matrix_1.2-18                           cli_2.0.2                               gdata_2.18.0                           
# [103] igraph_1.2.5                            pkgconfig_2.0.3                         TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
# [106] rvcheck_0.1.8                           GenomicAlignments_1.18.1                plotly_4.9.2.1                         
# [109] xml2_1.3.2                              XVector_0.22.0                          rvest_0.3.5                            
# [112] digest_0.6.25                           sctransform_0.2.1                       RcppAnnoy_0.0.16                       
# [115] tsne_0.1-3                              graph_1.60.0                            Biostrings_2.50.2                      
# [118] cellranger_1.1.0                        leiden_0.3.3                            fastmatch_1.1-0                        
# [121] uwot_0.1.8                              Rsamtools_1.34.1                        gtools_3.8.2                           
# [124] graphite_1.28.2                         lifecycle_0.2.0                         nlme_3.1-145                           
# [127] jsonlite_1.6.1                          fansi_0.4.1                             viridisLite_0.3.0                      
# [130] pillar_1.4.4                            lattice_0.20-41                         httr_1.4.1                             
# [133] plotrix_3.7-8                           survival_3.2-3                          GO.db_3.7.0                            
# [136] glue_1.4.1                              UpSetR_1.4.0                            png_0.1-7                              
# [139] bit_1.1-15.2                            ggforce_0.3.2                           stringi_1.4.6                          
# [142] blob_1.2.1                              caTools_1.17.1.2                        memoise_1.1.0                          
# [145] irlba_2.3.3                             future.apply_1.5.0                      ape_5.3    
