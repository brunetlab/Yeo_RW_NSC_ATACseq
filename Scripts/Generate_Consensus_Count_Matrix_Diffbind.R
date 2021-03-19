#!/usr/bin/Rscript

#Robin Yeo
#04/02/2020

#This script uses Diffbind to generate a consensus peakset using all cell types (young + old).
#It then generates a raw count matrix (no normalization) and exports it for downstream analysis.
#This script also uses ChIPSeeker in order to generate consensus peaksets consisting of exclusively distal, intronic, and promoter peaks.

#This script then outputs the NSC consensus peakset and count matrix consisting of only Young and Old qNSC and aNSC peaks.

#This script additionally outputs a set of genome-wide promoters using ChIPSeeker (not related to my ATAC-sea peaks) for other downstream analyses


rm(list=ls())
library(data.table)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")

#Load in all of the replicates using a tab-delimited .csv file that contains information about file paths.
#Here the csv file contains all 25 high quality libraries (2-3 per condition) using the replicated peaks (i.e. young qNSC replicates would all share the same peak file path)
#from the Kundaje pipeline and individual bam files for each library (i.e. each library has a unique bam file path).

allPeaks <- dba(sampleSheet="Original_Data/All_Peaks_HQ_Replicated_Peaks_Individual.csv") 

#Prints data summary of allpeaks: number of peaks in each peakset, as well as the total number of unique peaks after merging overlapping ones
print(allPeaks) 

# 25 Samples, 141970 sites in matrix:
#   ID      Tissue        Factor Condition Replicate Caller Intervals
# 1  Y_Endo_1 Endothelial  Non-Dividing     Young         1    bed     28280
# 2  Y_qNSC_1        qNSC  Non-Dividing     Young         1    bed     70705
# 3  Y_aNSC_1        aNSC Proliferative     Young         1    bed     74840
# 4   Y_NPC_1         NPC Proliferative     Young         1    bed     91406
# 5  O_Endo_1 Endothelial  Non-Dividing       Old         1    bed     20237
# 6   O_Ast_1   Astrocyte  Non-Dividing       Old         1    bed     68576
# 7  O_qNSC_1        qNSC  Non-Dividing       Old         1    bed     69525
# 8  O_aNSC_1        aNSC Proliferative       Old         1    bed     88383
# 9   O_NPC_1         NPC Proliferative       Old         1    bed     69025
# 10 Y_Endo_2 Endothelial  Non-Dividing     Young         2    bed     28280
# 11  Y_Ast_2   Astrocyte  Non-Dividing     Young         2    bed     71458
# 12 Y_qNSC_2        qNSC  Non-Dividing     Young         2    bed     70705
# 13  Y_NPC_2         NPC Proliferative     Young         2    bed     91406
# 14 O_Endo_2 Endothelial  Non-Dividing       Old         2    bed     20237
# 15  O_Ast_2   Astrocyte  Non-Dividing       Old         2    bed     68576
# 16 O_qNSC_2        qNSC  Non-Dividing       Old         2    bed     69525
# 17 O_aNSC_2        aNSC Proliferative       Old         2    bed     88383
# 18  Y_Ast_3   Astrocyte  Non-Dividing     Young         3    bed     71458
# 19 Y_qNSC_3        qNSC  Non-Dividing     Young         3    bed     70705
# 20 Y_aNSC_3        aNSC Proliferative     Young         3    bed     74840
# 21  Y_NPC_3         NPC Proliferative     Young         3    bed     91406
# 22  O_Ast_3   Astrocyte  Non-Dividing       Old         3    bed     68576
# 23 O_qNSC_3        qNSC  Non-Dividing       Old         3    bed     69525
# 24 O_aNSC_3        aNSC Proliferative       Old         3    bed     88383
# 25  O_NPC_3         NPC Proliferative       Old         3    bed     69025

#Calculate a binding matrix with scores based on read counts for every sample
allPeaks.count.Raw = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

#This converts allPeaks.count.Raw to a dataframe for writing out tables
countMatrix.Raw <- dba.peakset(allPeaks.count.Raw, bRetrieve=TRUE, DataType = DBA_DATA_FRAME)

#Exporting the count matrix to a txt file for downstream use
write.table(countMatrix.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/consensus_countMatrix_Raw.txt", sep="\t", row.names=FALSE)

#Exporting only the peaks (not sample count matrix values)
write.table(countMatrix.Raw[,c(1:3)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/consensus_peakset.txt", sep="\t", row.names=FALSE)

# Annotating each peak in the global consensus count matrix with a genomic identity to output count matrices subsetted by promoter to correlate with gene expression downstream
countMatrix.Raw <- dba.peakset(allPeaks.count.Raw, bRetrieve=TRUE)
countMatrix.Raw <- as.data.frame(annotatePeak(countMatrix.Raw, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db"))

#Exporting the annotated count matrix to a txt file for downstream use
write.table(countMatrix.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/consensus_countMatrix_Raw_Annotated.txt", sep="\t", row.names=FALSE)

#Outputting an annotated count matrix, an unannotated count matrix, and the peakset based on global consensus promoter peaks
countMatrix.promoter.Raw <- countMatrix.Raw[countMatrix.Raw$annotation %like% "Promoter",]
write.table(countMatrix.promoter.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_consensus_countMatrix_Raw_Annotated.txt", sep="\t", row.names=FALSE)
write.table(countMatrix.promoter.Raw[,c(1:3,6:30)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_consensus_countMatrix_Raw.txt", sep="\t", row.names=FALSE) #removing rows 4 and 5 which are "width" and "strand"
write.table(countMatrix.promoter.Raw[,c(1:3)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_consensus_peakset.txt", sep="\t", row.names=FALSE)





##---------------------------------------------------------------------------------------
#Generating Young/Old qNSC/aNSC consensus count matrix
##---------------------------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")


#Load in all of the replicates using a tab-delimited .csv file that contains information about file paths.
#Here the csv file contains all 25 high quality libraries (2-3 per condition) using the replicated peaks (i.e. young qNSC replicates would all share the same peak file path)
#from the Kundaje pipeline and individual bam files for each library (i.e. each library has a unique bam file path).

allPeaks <- dba(sampleSheet="Original_Data/Y_O_qNSC_aNSC_peaks.csv") 

#Prints data summary of allpeaks: number of peaks in each peakset, as well as the total number of unique peaks after merging overlapping ones
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


#Calculate a binding matrix with scores based on read counts for every sample
allPeaks.count.Raw = dba.count(allPeaks, minOverlap=0, bParallel=TRUE, score=DBA_SCORE_READS, bRemoveDuplicates=FALSE)

#This converts allPeaks.count.Raw to a dataframe for writing out tables
countMatrix.Raw <- dba.peakset(allPeaks.count.Raw, bRetrieve=TRUE, DataType = DBA_DATA_FRAME)

#Exporting the count matrix to a txt file for downstream use
write.table(countMatrix.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/Y_O_qNSC_aNSC_countMatrix_Raw.txt", sep="\t", row.names=FALSE)

#Exporting only the peaks (not sample count matrix values)
write.table(countMatrix.Raw[,c(1:3)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/Y_O_qNSC_aNSC_peakset.txt", sep="\t", row.names=FALSE)

# Annotating each peak in the NSC peakset with a genomic identity to output count matrices subsetted by promoter, intron, distal region, etc.
# For each genomic subset, I am outputting an annotated count matrix, an unannotated count matrix, and the peakset
countMatrix.Raw <- dba.peakset(allPeaks.count.Raw, bRetrieve=TRUE)
countMatrix.Raw <- as.data.frame(annotatePeak(countMatrix.Raw, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db"))

#Exporting the annotated count matrix to a txt file for downstream use
write.table(countMatrix.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/Y_O_qNSC_aNSC_countMatrix_Raw_Annotated.txt", sep="\t", row.names=FALSE)

#Promoters: 20633 peaks
countMatrix.promoter.Raw <- countMatrix.Raw[countMatrix.Raw$annotation %like% "Promoter",]
write.table(countMatrix.promoter.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_Y_O_qNSC_aNSC_countMatrix_Raw_Annotated.txt", sep="\t", row.names=FALSE)
write.table(countMatrix.promoter.Raw[,c(1:3,6:16)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_Y_O_qNSC_aNSC_countMatrix_Raw.txt", sep="\t", row.names=FALSE) #removing rows 4 and 5 which are "width" and "strand"
write.table(countMatrix.promoter.Raw[,c(1:3)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/promoter_Y_O_qNSC_aNSC_peakset.txt", sep="\t", row.names=FALSE)

#Introns: 28571 peaks
countMatrix.intronic.Raw <- countMatrix.Raw[countMatrix.Raw$annotation %like% "Intron",]
write.table(countMatrix.intronic.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/intronic_Y_O_qNSC_aNSC_countMatrix_Raw_Annotated.txt", sep="\t", row.names=FALSE)
write.table(countMatrix.intronic.Raw[,c(1:3,6:16)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/intronic_Y_O_qNSC_aNSC_countMatrix_Raw.txt", sep="\t", row.names=FALSE) #removing rows 4 and 5 which are "width" and "strand"
write.table(countMatrix.intronic.Raw[,c(1:3)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/intronic_Y_O_qNSC_aNSC_peakset.txt", sep="\t", row.names=FALSE)

#Distal Intergenic: 31660 peaks
countMatrix.distal.Raw <- countMatrix.Raw[countMatrix.Raw$annotation == "Distal Intergenic",]
write.table(countMatrix.distal.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/distal_Y_O_qNSC_aNSC_countMatrix_Raw_Annotated.txt", sep="\t", row.names=FALSE)
write.table(countMatrix.distal.Raw[,c(1:3,6:16)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/distal_Y_O_qNSC_aNSC_countMatrix_Raw.txt", sep="\t", row.names=FALSE) #removing rows 4 and 5 which are "width" and "strand"
write.table(countMatrix.distal.Raw[,c(1:3)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/distal_Y_O_qNSC_aNSC_peakset.txt", sep="\t", row.names=FALSE)

#Distal Intergenic + Intronic: 60231 peaks
countMatrix.DI.Raw <- countMatrix.Raw[countMatrix.Raw$annotation == "Distal Intergenic" | countMatrix.Raw$annotation %like% "Intron",]
write.table(countMatrix.DI.Raw, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/DI_Y_O_qNSC_aNSC_countMatrix_Raw_Annotated.txt", sep="\t", row.names=FALSE)
write.table(countMatrix.DI.Raw[,c(1:3,6:16)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/DI_Y_O_qNSC_aNSC_countMatrix_Raw.txt", sep="\t", row.names=FALSE) #removing rows 4 and 5 which are "width" and "strand"
write.table(countMatrix.DI.Raw[,c(1:3)], file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/DI_Y_O_qNSC_aNSC_peakset.txt", sep="\t", row.names=FALSE)


##---------------------------------------------------------------------------------------
#GENERATING COUNT MATRICES FOR ALL GENOME-WIDE PROMOTERS (NOTE: THIS IS NOT RELATED TO THE ATAC-SEQ DATA)
#This outputs two peaksets consisting of 24,244 promoter peaks (using either +/- 1kb or -3000bp/+100bp)
##---------------------------------------------------------------------------------------

promoters_1k <- as.data.frame(getPromoters(TxDb=txdb, upstream=1000, downstream=1000))
promoters_3k <- as.data.frame(getPromoters(TxDb=txdb, upstream=3000, downstream=100))
write.table(promoters_1k, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/genome_wide_promoter_peaks_min1k_plus1k.txt", sep="\t", row.names=FALSE)
write.table(promoters_3k, file="Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/genome_wide_promoter_peaks_min3k_plus100bp.txt", sep="\t", row.names=FALSE)


# RWY SESSION_INFO()
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
#   [1] TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4 GenomicFeatures_1.34.8                   AnnotationDbi_1.44.0                     ChIPseeker_1.18.0                        DiffBind_2.10.0                         
# [6] SummarizedExperiment_1.12.0              DelayedArray_0.8.0                       BiocParallel_1.16.6                      matrixStats_0.56.0                       Biobase_2.42.0                          
# [11] GenomicRanges_1.34.0                     GenomeInfoDb_1.18.2                      IRanges_2.16.0                           S4Vectors_0.20.1                         BiocGenerics_0.28.0                     
# [16] data.table_1.12.8                       
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.1.8                         GOstats_2.48.0                          fastmatch_1.1-0                         plyr_1.8.6                              igraph_1.2.5                           
# [6] GSEABase_1.44.0                         splines_3.5.2                           BatchJobs_1.8                           gridBase_0.4-7                          ggplot2_3.3.2                          
# [11] urltools_1.7.3                          amap_0.8-16                             digest_0.6.25                           GOSemSim_2.8.0                          viridis_0.5.1                          
# [16] GO.db_3.7.0                             gdata_2.18.0                            magrittr_1.5                            checkmate_2.0.0                         memoise_1.1.0                          
# [21] BBmisc_1.11                             limma_3.38.3                            Biostrings_2.50.2                       annotate_1.60.1                         graphlayouts_0.7.0                     
# [26] systemPipeR_1.16.1                      enrichplot_1.2.0                        prettyunits_1.1.1                       colorspace_1.4-1                        blob_1.2.1                             
# [31] ggrepel_0.8.2                           dplyr_1.0.0                             crayon_1.3.4                            RCurl_1.98-1.2                          jsonlite_1.6.1                         
# [36] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 graph_1.60.0                            genefilter_1.64.0                       brew_1.0-6                              survival_3.2-3                         
# [41] sendmailR_1.2-1                         glue_1.4.1                              polyclip_1.10-0                         gtable_0.3.0                            zlibbioc_1.28.0                        
# [46] XVector_0.22.0                          UpSetR_1.4.0                            Rgraphviz_2.26.0                        scales_1.1.1                            DOSE_3.8.2                             
# [51] pheatmap_1.0.12                         DBI_1.1.0                               edgeR_3.24.3                            Rcpp_1.0.4.6                            plotrix_3.7-8                          
# [56] viridisLite_0.3.0                       xtable_1.8-4                            progress_1.2.2                          gridGraphics_0.5-0                      bit_1.1-15.2                           
# [61] europepmc_0.4                           AnnotationForge_1.24.0                  httr_1.4.1                              fgsea_1.8.0                             gplots_3.0.3                           
# [66] RColorBrewer_1.1-2                      ellipsis_0.3.1                          pkgconfig_2.0.3                         XML_3.99-0.3                            farver_2.0.3                           
# [71] locfit_1.5-9.4                          ggplotify_0.0.5                         tidyselect_1.1.0                        rlang_0.4.6                             reshape2_1.4.4                         
# [76] munsell_0.5.0                           tools_3.5.2                             generics_0.0.2                          RSQLite_2.2.0                           ggridges_0.5.2                         
# [81] stringr_1.4.0                           yaml_2.2.1                              bit64_0.9-7                             tidygraph_1.2.0                         caTools_1.17.1.2                       
# [86] purrr_0.3.4                             ggraph_2.0.3                            RBGL_1.58.2                             DO.db_2.9                               xml2_1.3.2                             
# [91] biomaRt_2.38.0                          compiler_3.5.2                          rstudioapi_0.11                         tibble_3.0.1                            tweenr_1.0.1                           
# [96] stringi_1.4.6                           lattice_0.20-41                         Matrix_1.2-18                           vctrs_0.3.1                             pillar_1.4.4                           
# [101] lifecycle_0.2.0                         BiocManager_1.30.10                     triebeard_0.3.0                         cowplot_1.0.0                           bitops_1.0-6                           
# [106] rtracklayer_1.42.2                      qvalue_2.14.1                           R6_2.4.1                                latticeExtra_0.6-28                     hwriter_1.3.2                          
# [111] ShortRead_1.40.0                        KernSmooth_2.23-16                      gridExtra_2.3                           boot_1.3-25                             MASS_7.3-51.6                          
# [116] gtools_3.8.2                            Category_2.48.1                         rjson_0.2.20                            GenomicAlignments_1.18.1                Rsamtools_1.34.1                       
# [121] GenomeInfoDbData_1.2.0                  hms_0.5.3                               grid_3.5.2                              tidyr_1.1.0                             rvcheck_0.1.8                          
# [126] ggforce_0.3.2                           base64enc_0.1-3                        