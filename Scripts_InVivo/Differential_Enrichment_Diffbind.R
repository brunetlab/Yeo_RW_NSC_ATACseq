#!/usr/bin/Rscript

#Robin Yeo
#07/24/2018

#This script is to use Diffbind to generate lists of differentially accessible peaks (young vs. old) for each of the 5 cell types that are annotated with nearby genes (for downstream GO/KEGG enrichment).

rm(list=ls())
library(DiffBind)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


#This for loop outputs a list of chromatin peaks that are differentially open in both the young and old conditions for each of the 5 sorted cell types
#The differential peaks are used downstream for motif enrichment and the annotated genes are used downstream for GO/KEGG pathway enrichment.
cell_types <- c("Endo", "Ast","qNSC", "aNSC", "NPC")
for (cell_type in cell_types) {
threshold=0.05
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")

allPeaks <- dba(sampleSheet=paste("Original_Data/",cell_type,"_peaks.csv",sep=""))
print(allPeaks) 

allPeaks.count = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)
uncorrected.contrast = dba.contrast(allPeaks.count, categories=DBA_CONDITION, minMembers=2)
allPeaks.contrast.analyze = dba.analyze(uncorrected.contrast, bCorPlot=FALSE, bParallel=TRUE, bTagwise=FALSE, bFullLibrarySize=TRUE, bReduceObjects=FALSE,method=DBA_EDGER) 

#Outputting the full list of peaks from the consensus peakset with associated FDRs and FCs from differential peak calling
rep_all <- dba.report(allPeaks.contrast.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report_all <- annotatePeak(makeGRangesFromDataFrame(rep_all,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(as.data.frame(Report_all),file=paste("Output_Data/Differential_Peaks/All_Diff_peaks/",cell_type, "_Differential_Sites_withFC.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#Outputting the list of peaks that are differentially enriched in young or old conditions with FDR < 0.05
rep_thresh <- dba.report(allPeaks.contrast.analyze, contrast=1, th=threshold, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report <- annotatePeak(makeGRangesFromDataFrame(rep_thresh,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
Report <- as.data.frame(Report)
indexes = Report$Fold>0
Young_Enriched = Report[indexes,]
Old_Enriched = Report[!indexes,]

write.table(Young_Enriched, file=paste("Output_Data/Differential_Peaks/",cell_type, "_Young_Differential_Sites_withFC_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(Old_Enriched, file=paste("Output_Data/Differential_Peaks/",cell_type, "_Old_Differential_Sites_withFC_FDR_",threshold,".narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#RWY's outputs from print(allPeaks)

# 4 Samples, 33550 sites in matrix:
#   ID      Tissue       Factor Condition Replicate Caller Intervals
# 1 Y_Endo_1 Endothelial Non-Dividing     Young         1    bed     28280
# 2 Y_Endo_2 Endothelial Non-Dividing     Young         2    bed     28280
# 3 O_Endo_1 Endothelial Non-Dividing       Old         1    bed     20237
# 4 O_Endo_2 Endothelial Non-Dividing       Old         2    bed     20237
# 
# 5 Samples, 66961 sites in matrix:
#   ID    Tissue       Factor Condition Replicate Caller Intervals
# 1 Y_Ast_2 Astrocyte Non-Dividing     Young         2    bed     71458
# 2 Y_Ast_3 Astrocyte Non-Dividing     Young         3    bed     71458
# 3 O_Ast_1 Astrocyte Non-Dividing       Old         1    bed     68576
# 4 O_Ast_2 Astrocyte Non-Dividing       Old         2    bed     68576
# 5 O_Ast_3 Astrocyte Non-Dividing       Old         3    bed     68576
# 
# 6 Samples, 64193 sites in matrix:
#   ID Tissue       Factor Condition Replicate Caller Intervals
# 1 Y_qNSC_1   qNSC Non-Dividing     Young         1    bed     70705
# 2 Y_qNSC_2   qNSC Non-Dividing     Young         2    bed     70705
# 3 Y_qNSC_3   qNSC Non-Dividing     Young         3    bed     70705
# 4 O_qNSC_1   qNSC Non-Dividing       Old         1    bed     69525
# 5 O_qNSC_2   qNSC Non-Dividing       Old         2    bed     69525
# 6 O_qNSC_3   qNSC Non-Dividing       Old         3    bed     69525
# 
# 5 Samples, 58801 sites in matrix:
#   ID Tissue        Factor Condition Replicate Caller Intervals
# 1 Y_aNSC_1   aNSC Proliferative     Young         1    bed     74840
# 2 Y_aNSC_3   aNSC Proliferative     Young         3    bed     74840
# 3 O_aNSC_1   aNSC Proliferative       Old         1    bed     88383
# 4 O_aNSC_2   aNSC Proliferative       Old         2    bed     88383
# 5 O_aNSC_3   aNSC Proliferative       Old         3    bed     88383
# 
# 5 Samples, 55421 sites in matrix:
#   ID Tissue        Factor Condition Replicate Caller Intervals
# 1 Y_NPC_1    NPC Proliferative     Young         1    bed     91406
# 2 Y_NPC_2    NPC Proliferative     Young         2    bed     91406
# 3 Y_NPC_3    NPC Proliferative     Young         3    bed     91406
# 4 O_NPC_1    NPC Proliferative       Old         1    bed     69025
# 5 O_NPC_3    NPC Proliferative       Old         3    bed     69025


###################################################################################################
### The following outputs differential peaks for qNSC vs. aNSC in either young or old conditions
###################################################################################################

rm(list=ls())
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
threshold=0.05
setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking")
allPeaks <- dba(sampleSheet="Original_Data/Y_O_qNSC_aNSC_peaks.csv")
print(allPeaks) 

# 11 Samples, 87796 sites in matrix:
#   ID Tissue Factor Condition Replicate Caller Intervals
# 1  Y_qNSC_1   qNSC Y_qNSC     Young         1    bed     70705
# 2  Y_aNSC_1   aNSC Y_aNSC     Young         1    bed     74840
# 3  O_qNSC_1   qNSC O_qNSC       Old         1    bed     69525
# 4  O_aNSC_1   aNSC O_aNSC       Old         1    bed     88383
# 5  Y_qNSC_2   qNSC Y_qNSC     Young         2    bed     70705
# 6  O_qNSC_2   qNSC O_qNSC       Old         2    bed     69525
# 7  O_aNSC_2   aNSC O_aNSC       Old         2    bed     88383
# 8  Y_qNSC_3   qNSC Y_qNSC     Young         3    bed     70705
# 9  Y_aNSC_3   aNSC Y_aNSC     Young         3    bed     74840
# 10 O_qNSC_3   qNSC O_qNSC       Old         3    bed     69525
# 11 O_aNSC_3   aNSC O_aNSC       Old         3    bed     88383

allPeaks.count = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

#This calculates differential peaks during differentiation of YOUNG cell types
uncorrected.contrast = dba.contrast(allPeaks.count, group1=allPeaks.count$masks$Y_qNSC,group2=allPeaks.count$masks$Y_aNSC, name1="Y_qNSC",name2="Y_aNSC", minMembers=2)
Y.allPeaks.contrast.analyze = dba.analyze(uncorrected.contrast, bCorPlot=FALSE, bParallel=TRUE, bTagwise=FALSE, bFullLibrarySize=TRUE, bReduceObjects=FALSE,method=DBA_EDGER) 

rep_all_Y <- dba.report(Y.allPeaks.contrast.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report_all_Y <- annotatePeak(makeGRangesFromDataFrame(rep_all_Y,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(as.data.frame(Report_all_Y),file=paste("Output_Data/Differential_Peaks/All_Diff_peaks/Young_qNSC_aNSC_Differential_Sites_withFC.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

rep_Y <- dba.report(Y.allPeaks.contrast.analyze, contrast=1, th=threshold, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report <- annotatePeak(makeGRangesFromDataFrame(rep_Y,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
Report <- as.data.frame(Report)

indexes = Report$Fold>0
Q_Enriched = Report[indexes,]
A_Enriched = Report[!indexes,]

write.table(Q_Enriched, file=paste("Output_Data/Differential_Peaks/qNSC_Young_Q_A_Differential_Report_FDR_0.05.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(A_Enriched, file=paste("Output_Data/Differential_Peaks/aNSC_Young_Q_A_Differential_Report_FDR_0.05.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


#This calculates differential peaks during differentiation of OLD cell types
uncorrected.contrast = dba.contrast(allPeaks.count, group1=allPeaks.count$masks$O_qNSC,group2=allPeaks.count$masks$O_aNSC, name1="O_qNSC",name2="O_aNSC", minMembers=2)
O.allPeaks.contrast.analyze = dba.analyze(uncorrected.contrast, bCorPlot=FALSE, bParallel=TRUE, bTagwise=FALSE, bFullLibrarySize=TRUE, bReduceObjects=FALSE,method=DBA_EDGER) 

rep_all_O <- dba.report(O.allPeaks.contrast.analyze, contrast=1, th=1, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report_all_O <- annotatePeak(makeGRangesFromDataFrame(rep_all_O,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(as.data.frame(Report_all_O),file=paste("Output_Data/Differential_Peaks/All_Diff_peaks/Old_qNSC_aNSC_Differential_Sites_withFC.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

rep_O <- dba.report(O.allPeaks.contrast.analyze, contrast=1, th=threshold, bCount=TRUE, method=DBA_EDGER, DataType=DBA_DATA_FRAME)
Report <- annotatePeak(makeGRangesFromDataFrame(rep_O,keep.extra.columns=TRUE), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
Report <- as.data.frame(Report)

indexes = Report$Fold>0
Q_Enriched = Report[indexes,]
A_Enriched = Report[!indexes,]

write.table(Q_Enriched, file=paste("Output_Data/Differential_Peaks/qNSC_Old_Q_A_Differential_Report_FDR_0.05.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(A_Enriched, file=paste("Output_Data/Differential_Peaks/aNSC_Old_Q_A_Differential_Report_FDR_0.05.narrowPeak",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


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