
#install.packages("calibrate")
#install.packages("pheatmap")

rm(list=ls())

library(calibrate)
library(pheatmap)
library(RColorBrewer)

library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(ChIPseeker)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")

res <- read.table("Output_Data/Differential_Peaks/All_Diff_peaks/Young_qNSC_aNSC_Differential_Sites_withFC.narrowPeak", header=TRUE, sep="\t", quote="")
res_diff <- res[res$FDR<0.05,]
res_diff_Q <- res_diff[res_diff$Fold>0,]
res_diff_A <- res_diff[res_diff$Fold<0,]


pdf("Output_Figs/Y_qNSC_aNSC_Volcano_plot.pdf")
# Make a basic volcano plot
#NOTE: I am plotting "-Fold" here since Fold>0 corresponds to qNSCs and I want aNSCs to be on the right of the plot
with(res, plot(-Fold, -log10(FDR), pch=20, main="Volcano plot", xlim=c(-6,6), ylim=c(0,7), grid()))
with(subset(res, FDR<.05 & Fold>0), points(-Fold, -log10(FDR), pch=20, col="dodgerblue")) #qNSCs (Fold > 0)
with(subset(res, FDR<.05 & Fold<0), points(-Fold, -log10(FDR), pch=20, col="red2")) #aNSCs (Fold < 0)
dev.off()

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
#   [1] ChIPseeker_1.18.0                        ReactomePA_1.26.0                        clusterProfiler_3.10.1                  
# [4] ggplot2_3.3.2                            org.Mm.eg.db_3.7.0                       TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4
# [7] GenomicFeatures_1.34.8                   AnnotationDbi_1.44.0                     Biobase_2.42.0                          
# [10] GenomicRanges_1.34.0                     GenomeInfoDb_1.18.2                      IRanges_2.16.0                          
# [13] S4Vectors_0.20.1                         BiocGenerics_0.28.0                      RColorBrewer_1.1-2                      
# [16] pheatmap_1.0.12                          calibrate_1.7.7                          MASS_7.3-51.6                           
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.8.0                             colorspace_1.4-1                        ellipsis_0.3.1                         
# [4] ggridges_0.5.2                          qvalue_2.14.1                           XVector_0.22.0                         
# [7] rstudioapi_0.11                         farver_2.0.3                            urltools_1.7.3                         
# [10] graphlayouts_0.7.0                      ggrepel_0.8.2                           bit64_0.9-7                            
# [13] xml2_1.3.2                              splines_3.5.2                           GOSemSim_2.8.0                         
# [16] polyclip_1.10-0                         jsonlite_1.6.1                          Rsamtools_1.34.1                       
# [19] gridBase_0.4-7                          GO.db_3.7.0                             graph_1.60.0                           
# [22] ggforce_0.3.2                           graphite_1.28.2                         BiocManager_1.30.10                    
# [25] compiler_3.5.2                          httr_1.4.1                              rvcheck_0.1.8                          
# [28] backports_1.1.8                         Matrix_1.2-18                           TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
# [31] tweenr_1.0.1                            prettyunits_1.1.1                       tools_3.5.2                            
# [34] igraph_1.2.5                            gtable_0.3.0                            glue_1.4.1                             
# [37] GenomeInfoDbData_1.2.0                  reshape2_1.4.4                          DO.db_2.9                              
# [40] dplyr_1.0.0                             rappdirs_0.3.1                          fastmatch_1.1-0                        
# [43] Rcpp_1.0.4.6                            enrichplot_1.2.0                        vctrs_0.3.1                            
# [46] Biostrings_2.50.2                       gdata_2.18.0                            rtracklayer_1.42.2                     
# [49] ggraph_2.0.3                            stringr_1.4.0                           lifecycle_0.2.0                        
# [52] gtools_3.8.2                            XML_3.99-0.3                            DOSE_3.8.2                             
# [55] europepmc_0.4                           zlibbioc_1.28.0                         scales_1.1.1                           
# [58] tidygraph_1.2.0                         reactome.db_1.66.0                      hms_0.5.3                              
# [61] SummarizedExperiment_1.12.0             yaml_2.2.1                              memoise_1.1.0                          
# [64] gridExtra_2.3                           UpSetR_1.4.0                            biomaRt_2.38.0                         
# [67] triebeard_0.3.0                         stringi_1.4.6                           RSQLite_2.2.0                          
# [70] plotrix_3.7-8                           checkmate_2.0.0                         caTools_1.17.1.2                       
# [73] boot_1.3-25                             BiocParallel_1.16.6                     rlang_0.4.6                            
# [76] pkgconfig_2.0.3                         matrixStats_0.56.0                      bitops_1.0-6                           
# [79] lattice_0.20-41                         purrr_0.3.4                             GenomicAlignments_1.18.1               
# [82] cowplot_1.0.0                           bit_1.1-15.2                            tidyselect_1.1.0                       
# [85] plyr_1.8.6                              magrittr_1.5                            R6_2.4.1                               
# [88] gplots_3.0.3                            generics_0.0.2                          DelayedArray_0.8.0                     
# [91] DBI_1.1.0                               pillar_1.4.4                            withr_2.2.0                            
# [94] RCurl_1.98-1.2                          tibble_3.0.1                            crayon_1.3.4                           
# [97] KernSmooth_2.23-16                      viridis_0.5.1                           progress_1.2.2                         
# [100] grid_3.5.2                              data.table_1.12.8                       blob_1.2.1                             
# [103] digest_0.6.25                           tidyr_1.1.0                             gridGraphics_0.5-0                     
# [106] munsell_0.5.0                           viridisLite_0.3.0                       ggplotify_0.0.5  