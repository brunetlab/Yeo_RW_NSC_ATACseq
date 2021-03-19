#06_18_2018
#R version 3.5.2 (2018-12-20)

#The purpose of this script is to use ChIPseeker (v1.18.0) to annotate ATAC-seq peaks and integrate peaks with RNA-seq data.
#Specifically, I want to see if promoter accessibility is correlated with expression.

#This script reads in the VST-normalized count matrix of all 24,887 promoter peaks from the consensus peakset (note that for some genes there are multiple promoter peaks)
#and associates each consensus peak with a gene then plots the gene expression vs accessibility

#There are 24,887 consensus promoter peaks, 15,329 genes in the expression matrix and 19,156 gene/promoter pairs (note that some genes have multiple promoters)

#Note: approx ~100-200 points per plot are truncated due to the y-axis range of (5-15) I've chosen
#This appears as a Warning message e.g. "Removed 166 rows containing non-finite values (stat_boxplot)."


rm(list=ls())
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")

#This reads in a generated peakset of all promoters (outputted by "Generate_Consensus_Count_Matrix_Diffbind.R") and calculates a count matrix for the cell_type reads that fall within promoter peaks
promoter_peaks <- readPeakFile("Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_consensus_countMatrix_VST.txt", header=TRUE)

promoter_peaks.count <- as.data.frame(annotatePeak(promoter_peaks, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db"))

#This reads in the averaged RNA-seq values for each cell type
expression_values <- read.table("Original_Data/Leeman_RNAseq_average_expression_values.txt")
colnames(expression_values)[colnames(expression_values)=="Gene_Name"] <- "SYMBOL"

#Merging both dataframes based on gene IDs so that I can compare promoter accessibility to expression levels
merged_df <- merge(promoter_peaks.count,expression_values,by="SYMBOL")

labels <- list("Y_Ast_ATAC", "Y_qNSC_ATAC", "Y_aNSC_ATAC", "Y_NPC_ATAC", "Y_Endo_ATAC", "O_Ast_ATAC", "O_qNSC_ATAC", "O_aNSC_ATAC", "O_NPC_ATAC", "O_Endo_ATAC")
#names(high_accessibility_promoters) <- labels
merged_df$Y_Ast_ATAC <- rowMeans(cbind(merged_df$Y_Ast_2,merged_df$Y_Ast_3))
merged_df$Y_qNSC_ATAC <- rowMeans(cbind(merged_df$Y_qNSC_1,merged_df$Y_qNSC_2,merged_df$Y_qNSC_3))
merged_df$Y_aNSC_ATAC <- rowMeans(cbind(merged_df$Y_aNSC_1,merged_df$Y_aNSC_3))
merged_df$Y_NPC_ATAC <- rowMeans(cbind(merged_df$Y_NPC_1,merged_df$Y_NPC_2,merged_df$Y_NPC_3))
merged_df$Y_Endo_ATAC <- rowMeans(cbind(merged_df$Y_Endo_1,merged_df$Y_Endo_2))

merged_df$O_Ast_ATAC <- rowMeans(cbind(merged_df$O_Ast_1,merged_df$O_Ast_2,merged_df$O_Ast_3))
merged_df$O_qNSC_ATAC <- rowMeans(cbind(merged_df$O_qNSC_1,merged_df$O_qNSC_2,merged_df$O_qNSC_3))
merged_df$O_aNSC_ATAC <- rowMeans(cbind(merged_df$O_aNSC_1,merged_df$O_aNSC_2,merged_df$O_aNSC_3))
merged_df$O_NPC_ATAC <- rowMeans(cbind(merged_df$O_NPC_1,merged_df$O_NPC_3))
merged_df$O_Endo_ATAC <- rowMeans(cbind(merged_df$O_Endo_1,merged_df$O_Endo_2))

## --------------------------------------------------------------------------
#Plotting deciled boxplots to correlate promoter accessibility with gene expression (from Dena's bulk RNA-seq data)

library(OneR)

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black",size=10),axis.text.y=element_text(colour="black",size=15),
             axis.title.x=element_text(colour="black",size=15,vjust=-1),axis.title.y=element_text(colour="black",size=15),
             axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

deciles <- c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")

#Y_qNSC lightskyblue
#O_qNSC dodgerblue
#Y_aNSC orange
#O_aNSC red2

Y_qNSC_bin <- bin(merged_df$Y_qNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=Y_qNSC_bin,y=merged_df$Y_qNSC))
gg + 
  geom_boxplot(fill="lightskyblue") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  ylim(5,15) +
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/Y_qNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

O_qNSC_bin <- bin(merged_df$O_qNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=O_qNSC_bin,y=merged_df$O_qNSC))
gg + 
  geom_boxplot(fill="dodgerblue") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  ylim(5,15) +
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/O_qNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

Y_aNSC_bin <- bin(merged_df$Y_aNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=Y_aNSC_bin,y=merged_df$Y_aNSC))
gg + 
  geom_boxplot(fill="orange") + 
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  ylim(5,15) +
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/Y_aNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)

O_aNSC_bin <- bin(merged_df$O_aNSC_ATAC, nbins = 10, labels = NULL, method = c("content"), na.omit = TRUE)
gg <- ggplot(data = merged_df, mapping = aes(x=O_aNSC_bin,y=merged_df$O_aNSC))
gg + 
  geom_boxplot(fill="red2") +
  labs(y="Gene expression", x='Promoter Accessibility Percentiles') +
  ylim(5,15) +
  guides(color=FALSE) +
  scale_x_discrete(labels= deciles) + 
  theme
ggsave(file="Output_Figs/O_aNSC_consensus_promoters_accessibility_vs_expression_BOXPLOT.png", width =11, height=11)


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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12                          ggplot2_3.3.2                            DESeq2_1.22.2                            org.Mm.eg.db_3.7.0                       TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4
# [6] GenomicFeatures_1.34.8                   AnnotationDbi_1.44.0                     ChIPseeker_1.18.0                        DiffBind_2.10.0                          SummarizedExperiment_1.12.0             
# [11] DelayedArray_0.8.0                       BiocParallel_1.16.6                      matrixStats_0.56.0                       Biobase_2.42.0                           GenomicRanges_1.34.0                    
# [16] GenomeInfoDb_1.18.2                      IRanges_2.16.0                           S4Vectors_0.20.1                         BiocGenerics_0.28.0                      data.table_1.12.8                       
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.1.8                         GOstats_2.48.0                          Hmisc_4.3-1                             fastmatch_1.1-0                         plyr_1.8.6                             
# [6] igraph_1.2.5                            GSEABase_1.44.0                         splines_3.5.2                           BatchJobs_1.8                           gridBase_0.4-7                         
# [11] urltools_1.7.3                          amap_0.8-16                             digest_0.6.25                           htmltools_0.5.0                         GOSemSim_2.8.0                         
# [16] viridis_0.5.1                           GO.db_3.7.0                             gdata_2.18.0                            magrittr_1.5                            checkmate_2.0.0                        
# [21] memoise_1.1.0                           BBmisc_1.11                             cluster_2.1.0                           limma_3.38.3                            Biostrings_2.50.2                      
# [26] annotate_1.60.1                         graphlayouts_0.7.0                      systemPipeR_1.16.1                      enrichplot_1.2.0                        prettyunits_1.1.1                      
# [31] colorspace_1.4-1                        blob_1.2.1                              ggrepel_0.8.2                           xfun_0.15                               dplyr_1.0.0                            
# [36] crayon_1.3.4                            RCurl_1.98-1.2                          jsonlite_1.6.1                          TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 graph_1.60.0                           
# [41] genefilter_1.64.0                       brew_1.0-6                              survival_3.2-3                          sendmailR_1.2-1                         glue_1.4.1                             
# [46] polyclip_1.10-0                         gtable_0.3.0                            zlibbioc_1.28.0                         XVector_0.22.0                          UpSetR_1.4.0                           
# [51] Rgraphviz_2.26.0                        scales_1.1.1                            DOSE_3.8.2                              DBI_1.1.0                               edgeR_3.24.3                           
# [56] Rcpp_1.0.4.6                            plotrix_3.7-8                           htmlTable_2.0.0                         viridisLite_0.3.0                       xtable_1.8-4                           
# [61] progress_1.2.2                          gridGraphics_0.5-0                      foreign_0.8-76                          bit_1.1-15.2                            europepmc_0.4                          
# [66] Formula_1.2-3                           AnnotationForge_1.24.0                  htmlwidgets_1.5.1                       httr_1.4.1                              fgsea_1.8.0                            
# [71] gplots_3.0.3                            RColorBrewer_1.1-2                      acepack_1.4.1                           ellipsis_0.3.1                          pkgconfig_2.0.3                        
# [76] XML_3.99-0.3                            farver_2.0.3                            nnet_7.3-14                             locfit_1.5-9.4                          labeling_0.3                           
# [81] ggplotify_0.0.5                         tidyselect_1.1.0                        rlang_0.4.6                             reshape2_1.4.4                          munsell_0.5.0                          
# [86] tools_3.5.2                             generics_0.0.2                          RSQLite_2.2.0                           ggridges_0.5.2                          stringr_1.4.0                          
# [91] yaml_2.2.1                              knitr_1.29                              bit64_0.9-7                             tidygraph_1.2.0                         caTools_1.17.1.2                       
# [96] purrr_0.3.4                             ggraph_2.0.3                            RBGL_1.58.2                             DO.db_2.9                               xml2_1.3.2                             
# [101] biomaRt_2.38.0                          compiler_3.5.2                          rstudioapi_0.11                         geneplotter_1.60.0                      tibble_3.0.1                           
# [106] tweenr_1.0.1                            stringi_1.4.6                           lattice_0.20-41                         Matrix_1.2-18                           vctrs_0.3.1                            
# [111] pillar_1.4.4                            lifecycle_0.2.0                         BiocManager_1.30.10                     triebeard_0.3.0                         cowplot_1.0.0                          
# [116] bitops_1.0-6                            rtracklayer_1.42.2                      qvalue_2.14.1                           R6_2.4.1                                latticeExtra_0.6-28                    
# [121] hwriter_1.3.2                           ShortRead_1.40.0                        KernSmooth_2.23-16                      gridExtra_2.3                           boot_1.3-25                            
# [126] MASS_7.3-51.6                           gtools_3.8.2                            Category_2.48.1                         rjson_0.2.20                            withr_2.2.0                            
# [131] GenomicAlignments_1.18.1                Rsamtools_1.34.1                        GenomeInfoDbData_1.2.0                  hms_0.5.3                               rpart_4.1-15                           
# [136] grid_3.5.2                              tidyr_1.1.0                             rvcheck_0.1.8                           ggforce_0.3.2                           base64enc_0.1-3
