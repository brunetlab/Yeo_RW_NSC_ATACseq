#!/bin/bash

#04/27/2020

#This script uses HOMER in order to find motif enrichment in narrowpeak files
#NOTE: I modified column 4 in the differential peak files ("*_FDR_0.05.narrowPeak") in order to give them unique IDs (since otherwise Homer thinks there are redundant peaks)



PEAKS_DIR="/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Output_Data_InVitro/Differential_Peaks/Homer_Peaks"
OUT_DIR="/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Output_Data_InVitro/Motifs_Output"
HOMER_PATH="/usr/local/lib/homer/bin"

MOTIF_ENRICHMENT="$HOMER_PATH/findMotifsGenome.pl"
GENOME_PATH="/usr/local/lib/homer/data/genomes/mm10"


cd ${HOMER_PATH}

for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size200_nobkg/${OFPREFIX}.motifs/" -size 200 -mask
done

for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size_gvn_nobkg/${OFPREFIX}.motifs/" -size given -mask
done


for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size200_YOQA/${OFPREFIX}.motifs/" -size 200 -mask -bg /Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Output_Data_InVitro/invitro_consensus_peakset_homer.txt
done

for PEAKS_FILE in $(find "$PEAKS_DIR" -name '*.narrowPeak')
do
	OFPREFIX=$(basename "${PEAKS_FILE}" | sed 's/\.narrowPeak//g')    

		perl $MOTIF_ENRICHMENT $PEAKS_FILE $GENOME_PATH "${OUT_DIR}/size_gvn_YOQA/${OFPREFIX}.motifs/" -size given -mask -bg /Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Output_Data_InVitro/invitro_consensus_peakset_homer.txt
done