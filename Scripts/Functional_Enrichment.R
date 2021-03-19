#07_25_2018
#Bioconductor version 3.6

#The purpose of this script is to implement ChIPseeker (v1.18.0) to annotate ATAC-seq peaks to determine the proportion of distal vs. promoter vs. introns etc...
#I output the distribution of different genomic elements for each pooled peakset then I manually create a table to plot it for Y/O Q/A

rm(list=ls())
library(ChIPseeker)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


setwd("~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/")

All_Cell_Types_Names <-list("Y_Astrocyte","Y_qNSC","Y_aNSC","Y_NPC","O_Endo","O_Astrocyte","O_aNSC","O_NPC")

for (i in c(1:8))
{
a <- readPeakFile(paste("Original_Data/FINAL_PEAK_FILES/", All_Cell_Types_Names[i], "_ppr.naive_overlap.filt.narrowPeak",sep=''))
a <- annotatePeak(a, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(a@annoStat, file=paste("Output_Data/Functional_Enrichment/", All_Cell_Types_Names[i], "_functional_enrichment.txt",sep=''), sep="\t", quote=FALSE)
}

#O_qNSC and Y_endo do not have the same peak file names (since these came from Optimal_Set) so I did these manually
a <- readPeakFile("Original_Data/FINAL_PEAK_FILES/O_qNSC_rep1-rep3.naive_overlap.filt.narrowPeak")
a <- annotatePeak(a, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(a@annoStat, file="Output_Data/Functional_Enrichment/O_qNSC_functional_enrichment.txt", sep="\t", quote=FALSE)

a <- readPeakFile("Original_Data/FINAL_PEAK_FILES/Y_Endo_rep1-rep2.naive_overlap.filt.narrowPeak")
a <- annotatePeak(a, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
write.table(a@annoStat, file="Output_Data/Functional_Enrichment/Y_Endo_functional_enrichment.txt", sep="\t", quote=FALSE)


#####I manually combine the outputs of the above to generate a .csv file called "ALL_functional_enrichment_condensed.csv" and plot it here

table <- read.csv("Output_Data/Functional_Enrichment/ALL_functional_enrichment_condensed.csv", header=TRUE, check.names=FALSE)
table$Conditions <- reorder(table$Conditions, rev(table$Order))
table_reshape <- melt(table[,c(1:8)], id.var = "Conditions")

table.NSC <- table[c(5:8),]
table.NSC$Conditions <- reorder(table.NSC$Conditions, rev(table.NSC$Order))
table_reshape.NSC <- melt(table.NSC[,c(1:8)], id.var = "Conditions")


theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size=rel(4)),axis.text.y=element_text(colour="black", size=rel(2)),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
             legend.title = element_text(color = "black", size = 0),legend.text = element_text(color = "black", size = 26))

table_reshape.NSC$Conditions <- factor(table_reshape.NSC$Conditions,levels = c("Y_qNSC", "O_qNSC", "Y_aNSC", "O_aNSC"))

pdf("Output_Figs/Functional_Enrichment_Barplot_NSC.pdf", height=12, width=18)
p <- ggplot(table_reshape.NSC, aes(Conditions, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_brewer(palette="RdYlBu") + guides(fill = guide_legend(reverse=TRUE)) + theme
p
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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2                       reshape_0.8.8                            OneR_2.2                                
# [4] pheatmap_1.0.12                          ggplot2_3.3.2                            DESeq2_1.22.2                           
# [7] org.Mm.eg.db_3.7.0                       TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4 GenomicFeatures_1.34.8                  
# [10] AnnotationDbi_1.44.0                     ChIPseeker_1.18.0                        DiffBind_2.10.0                         
# [13] SummarizedExperiment_1.12.0              DelayedArray_0.8.0                       BiocParallel_1.16.6                     
# [16] matrixStats_0.56.0                       Biobase_2.42.0                           GenomicRanges_1.34.0                    
# [19] GenomeInfoDb_1.18.2                      IRanges_2.16.0                           S4Vectors_0.20.1                        
# [22] BiocGenerics_0.28.0                      data.table_1.12.8                       
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.1.8                         GOstats_2.48.0                          Hmisc_4.3-1                            
# [4] fastmatch_1.1-0                         plyr_1.8.6                              igraph_1.2.5                           
# [7] GSEABase_1.44.0                         splines_3.5.2                           BatchJobs_1.8                          
# [10] gridBase_0.4-7                          urltools_1.7.3                          amap_0.8-16                            
# [13] digest_0.6.25                           htmltools_0.5.0                         GOSemSim_2.8.0                         
# [16] viridis_0.5.1                           GO.db_3.7.0                             gdata_2.18.0                           
# [19] magrittr_1.5                            checkmate_2.0.0                         memoise_1.1.0                          
# [22] BBmisc_1.11                             cluster_2.1.0                           limma_3.38.3                           
# [25] Biostrings_2.50.2                       annotate_1.60.1                         graphlayouts_0.7.0                     
# [28] systemPipeR_1.16.1                      enrichplot_1.2.0                        prettyunits_1.1.1                      
# [31] colorspace_1.4-1                        blob_1.2.1                              ggrepel_0.8.2                          
# [34] xfun_0.15                               dplyr_1.0.0                             crayon_1.3.4                           
# [37] RCurl_1.98-1.2                          jsonlite_1.6.1                          TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
# [40] graph_1.60.0                            genefilter_1.64.0                       brew_1.0-6                             
# [43] survival_3.2-3                          sendmailR_1.2-1                         glue_1.4.1                             
# [46] polyclip_1.10-0                         gtable_0.3.0                            zlibbioc_1.28.0                        
# [49] XVector_0.22.0                          UpSetR_1.4.0                            Rgraphviz_2.26.0                       
# [52] scales_1.1.1                            DOSE_3.8.2                              DBI_1.1.0                              
# [55] edgeR_3.24.3                            Rcpp_1.0.4.6                            plotrix_3.7-8                          
# [58] htmlTable_2.0.0                         viridisLite_0.3.0                       xtable_1.8-4                           
# [61] progress_1.2.2                          gridGraphics_0.5-0                      foreign_0.8-76                         
# [64] bit_1.1-15.2                            europepmc_0.4                           Formula_1.2-3                          
# [67] AnnotationForge_1.24.0                  htmlwidgets_1.5.1                       httr_1.4.1                             
# [70] fgsea_1.8.0                             gplots_3.0.3                            acepack_1.4.1                          
# [73] ellipsis_0.3.1                          pkgconfig_2.0.3                         XML_3.99-0.3                           
# [76] farver_2.0.3                            nnet_7.3-14                             locfit_1.5-9.4                         
# [79] labeling_0.3                            ggplotify_0.0.5                         tidyselect_1.1.0                       
# [82] rlang_0.4.6                             reshape2_1.4.4                          munsell_0.5.0                          
# [85] tools_3.5.2                             generics_0.0.2                          RSQLite_2.2.0                          
# [88] ggridges_0.5.2                          stringr_1.4.0                           yaml_2.2.1                             
# [91] knitr_1.29                              bit64_0.9-7                             tidygraph_1.2.0                        
# [94] caTools_1.17.1.2                        purrr_0.3.4                             ggraph_2.0.3                           
# [97] RBGL_1.58.2                             DO.db_2.9                               xml2_1.3.2                             
# [100] biomaRt_2.38.0                          compiler_3.5.2                          rstudioapi_0.11                        
# [103] geneplotter_1.60.0                      tibble_3.0.1                            tweenr_1.0.1                           
# [106] stringi_1.4.6                           lattice_0.20-41                         Matrix_1.2-18                          
# [109] vctrs_0.3.1                             pillar_1.4.4                            lifecycle_0.2.0                        
# [112] BiocManager_1.30.10                     triebeard_0.3.0                         cowplot_1.0.0                          
# [115] bitops_1.0-6                            rtracklayer_1.42.2                      qvalue_2.14.1                          
# [118] R6_2.4.1                                latticeExtra_0.6-28                     hwriter_1.3.2                          
# [121] ShortRead_1.40.0                        KernSmooth_2.23-16                      gridExtra_2.3                          
# [124] boot_1.3-25                             MASS_7.3-51.6                           gtools_3.8.2                           
# [127] Category_2.48.1                         rjson_0.2.20                            withr_2.2.0                            
# [130] GenomicAlignments_1.18.1                Rsamtools_1.34.1                        GenomeInfoDbData_1.2.0                 
# [133] hms_0.5.3                               rpart_4.1-15                            grid_3.5.2                             
# [136] tidyr_1.1.0                             rvcheck_0.1.8                           ggforce_0.3.2                          
# [139] base64enc_0.1-3     


