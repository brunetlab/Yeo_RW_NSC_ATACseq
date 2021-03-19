#06_05_2018
#Running ngs.plot.r in order to get TSS enrichment for young/old ATAC-seq data generated on Rep1 of paired young/old


OUT_DIR="~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/Output_Figs"
IN_DIR="~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking/Original_Data"
NGS_PATH="/usr/local/lib/ngsplot/bin"

cd ${NGS_PATH}


ngs.plot.r -G mm10 -R tss -C "${IN_DIR}/Multiplot_Pooled_ngs_tss.txt" -O "${OUT_DIR}/tssrefseq_Pooled" -D refseq