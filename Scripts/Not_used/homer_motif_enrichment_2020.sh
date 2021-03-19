#!/bin/bash

#04/23/2020

#This script uses HOMER in order to find motif enrichment in narrowpeak files
#NOTE: I modified column 4 in the differential peak files ("*_FDR_0.05.narrowPeak") in order to give them unique IDs (since otherwise Homer thinks there are redundant peaks)



PEAKS_DIR="/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/Output_Data/Differential_Peaks/Homer_Peaks"
OUT_DIR="/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/Output_Data/Motifs_Output"
HOMER_PATH="/usr/local/lib/homer/bin"

MOTIF_ENRICHMENT="$HOMER_PATH/findMotifsGenome.pl"
GENOME_PATH="/usr/local/lib/homer/data/genomes/mm10"


cd ${HOMER_PATH}

for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*_FDR_0.05.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\_FDR_0.05.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size200_nobkg/${OFPREFIX}.motifs/" -size 200 -mask
done

for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*_FDR_0.05.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\_FDR_0.05.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size_gvn_nobkg/${OFPREFIX}.motifs/" -size given -mask
done

for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*_FDR_0.05.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\_FDR_0.05.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size_gvn_consensus/${OFPREFIX}.motifs/" -size given -mask -bg /Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/consensus_peakset_homer.bed
done



for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*_FDR_0.05.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\_FDR_0.05.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size200_YOQA/${OFPREFIX}.motifs/" -size 200 -mask -bg /Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/Y_O_qNSC_aNSC_peakset_homer.bed
done

for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*_FDR_0.05.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\_FDR_0.05.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size_gvn_YOQA/${OFPREFIX}.motifs/" -size given -mask -bg /Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/Output_Data/FINAL_PEAKSETS_COUNT_MATRICES/Y_O_qNSC_aNSC_peakset_homer.bed
done