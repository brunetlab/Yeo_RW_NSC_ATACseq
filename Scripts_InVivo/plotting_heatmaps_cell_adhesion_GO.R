
#install.packages("calibrate")
#install.packages("pheatmap")


library(pheatmap)
library(RColorBrewer)


##########################################################################################
#Plotting heatmaps for various cell adhesion/migration GO pathways that were enriched in qNSCs (vs. aNSCs)
##########################################################################################

rm(list=ls())
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")
res <- read.table("Output_Data/Differential_Peaks/All_Diff_peaks/Young_qNSC_aNSC_Differential_Sites_withFC.narrowPeak", header=TRUE,quote="",sep="\t")
res_diff <- res[res$FDR<0.05,]


GO_lists <- list.files("Original_Data/GO_adhesion_gene_lists", all.files=FALSE)

#To subset to only differential peaks in cell adhesion pathway
for (pathway in GO_lists){

GO <- read.table(file=paste("Original_Data/GO_adhesion_gene_lists/",pathway,sep=""),sep="\t",header=TRUE)
res_diff_GO <- res_diff[res_diff$SYMBOL %in% GO$Symbol,]

res_diff_GO <- res_diff_GO[,c("Y_qNSC_1","Y_qNSC_2","Y_qNSC_3","Y_aNSC_1","Y_aNSC_3")]
breaksList <- seq (-2,2,by=0.04)

scaled_height <- nrow(res_diff_GO)/100
print(scaled_height)

pdf(paste("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/Output_Figs/Adhesion_heatmaps/",pathway,".pdf"), width=10, height = scaled_height)
pheatmap(res_diff_GO,
         scale = "row",
         display_numbers = FALSE,
         #color = colorRampPalette(brewer.pal(n = 6, name = "YlOrRd"))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         treeheight_row = 0,
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         show_rownames = FALSE,
         show_colnames = FALSE,
         border_color = NA)
dev.off()

}

#write.table(res_diff_GO, file="Output_Data/Differential_ATAC_peaks_adhesion_GO_young_Q_vs_A.txt",sep="\t", row.names=F,quote=FALSE)


##########################################################################################
#Plotting heatmap for differentially accessible peaks within GO Cell Adhesion category that change with age
##########################################################################################

rm(list=ls())
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")
res <- read.table("Output_Data/Differential_Peaks/All_Diff_peaks/qNSC_Differential_Sites_withFC.narrowPeak", header=TRUE,quote="",sep="\t")
res_diff <- res[res$FDR<=0.05,]

#To subset to only differential peaks in cell adhesion pathway
GO <- read.table(file="Original_Data/GO_Cell_Adhesion_GO0007155.txt",sep="\t",header=TRUE)
res_diff_GO <- res_diff[res_diff$SYMBOL %in% GO$Gene.name,]

res_diff_GO <- res_diff_GO[,c("Y_qNSC_1","Y_qNSC_2","Y_qNSC_3","O_qNSC_1","O_qNSC_2", "O_qNSC_3")]

breaksList <- seq (-2,2,by=0.04)

pdf(paste("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/Output_Figs/Adhesion_heatmaps/qNSC_Cell_Adhesion.pdf"))
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
#write.table(res_diff_GO, file="Output_Data/Differential_ATAC_peaks_adhesion_GO_qNSC_Y_vs_O.txt",sep="\t", row.names=F,quote=FALSE)


rm(list=ls()) 
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")
res <- read.table("Output_Data/Differential_Peaks/All_Diff_peaks/aNSC_Differential_Sites_withFC.narrowPeak", header=TRUE,quote="",sep="\t")
res_diff <- res[res$FDR<=0.05,]

#To subset to only differential peaks in cell adhesion pathway
GO <- read.table(file="Original_Data/GO_Cell_Adhesion_GO0007155.txt",sep="\t",header=TRUE)
res_diff_GO <- res_diff[res_diff$SYMBOL %in% GO$Gene.name,]

res_diff_GO <- res_diff_GO[,c("Y_aNSC_1","Y_aNSC_3","O_aNSC_1","O_aNSC_2", "O_aNSC_3")]

breaksList <- seq (-2,2,by=0.04)

pdf(paste("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/Output_Figs/Adhesion_heatmaps/aNSC_Cell_Adhesion.pdf"))
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
#write.table(res_diff_GO, file="Output_Data/Differential_ATAC_peaks_adhesion_GO_aNSC_Y_vs_O.txt",sep="\t", row.names=F,quote=FALSE)

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
