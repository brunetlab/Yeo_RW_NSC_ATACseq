#!/bin/bash

#04/24/2020

#This script will convert the de-duplicated, Tn5-shifted tagAlign (bed) files (from the multi-rep output of the in vitro libraries) 
#to bam format then sort them for downstream analysis

BED_DIR="/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Original_Data_InVitro/FINAL_TagAlign_FILES"
BAM_DIR="/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Original_Data_InVitro/FINAL_BAM_FILES"


cd ${BED_DIR}

for BED_FILE in $(find "$BED_DIR" -name '*.tagAlign')
do
	OFPREFIX=$(basename "${BED_FILE}" | sed 's/\.tagAlign//g')    

	bedToBam -i $BED_FILE -g "/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Original_Data_InVitro/mm10.chrom.sizes" > "${BAM_DIR}/${OFPREFIX}.bam"
	

done

cd ${BAM_DIR}

for BAM_FILE in $(find "$BAM_DIR" -name '*.tn5.bam')
do
	OFPREFIX=$(basename "${BAM_FILE}" | sed 's/\.tn5.bam//g')    


	samtools sort -o "${BAM_DIR}/${OFPREFIX}.sorted.tn5.bam" $BAM_FILE


done