#!/bin/bash

#07/19/2018

#This script will convert the de-duplicated, Tn5-shifted tagAlign (bed) files to bam format

BAM_DIR="/Volumes/brunetseq/Users/ryeo/Documents_brunetseq/YOUNG_OLD_ATAC_ANALYSIS_2018/FINAL_BAM_FILES"


cd ${BAM_DIR}

for BAM_FILE in $(find "$BAM_DIR" -name '*.tn5.bam')
do
	OFPREFIX=$(basename "${BAM_FILE}" | sed 's/\.tn5.bam//g')    


	samtools sort -o "${BAM_DIR}/${OFPREFIX}.sorted.tn5.bam" $BAM_FILE


done
