#4/29/2020

#Here I merge the high-quality (filtered,de-duplicated, NON-TN5-SHIFTED) bam files, sub-sample them, sort, and Tn5-shift them
#This is the same procedure that was done to create the in vivo pooled BAM files

# BAM_DIR="/Volumes/Genomic_TARDIS/IN_VITRO_ATAC_2019/Data/BAM_FILES/"
# cd ${BAM_DIR}

# #This manually merges the high-quality de-duplicated bam files
# samtools merge Y_qNSC_invitro_pooled.bam RY-vitro-1-TAAGGCGA_S1_R1_001.trim.PE2SE.nodup.bam RY-vitro-5-GGACTCCT_S5_R1_001.trim.PE2SE.nodup.bam
# samtools merge O_qNSC_invitro_pooled.bam RY-vitro-6-TAGGCATG_S6_R1_001.trim.PE2SE.nodup.bam RY-vitro-10-CGAGGCTG_S10_R1_001.trim.PE2SE.nodup.bam
# samtools merge Y_aNSC_invitro_pooled.bam RY-vitro-7-CTCTCTAC_S7_R1_001.trim.PE2SE.nodup.bam RY-vitro-11-AAGAGGCA_S11_R1_001.trim.PE2SE.nodup.bam RY-vitro-15-TGGATCTG_S15_R1_001.trim.PE2SE.nodup.bam
# samtools merge O_aNSC_invitro_pooled.bam RY-vitro-4-TCCTGAGC_S4_R1_001.trim.PE2SE.nodup.bam RY-vitro-8-CAGAGAGG_S8_R1_001.trim.PE2SE.nodup.bam RY-vitro-16-CCGTTTGT_S16_R1_001.trim.PE2SE.nodup.bam


#Here I manually transfer the newly-pooled BAM files from my external hard drive to the Dropbox folder

BAM_DIR="/Users/ryeo/Dropbox/RWY_ATAC_Code_Checking/IN_VITRO/Original_Data_InVitro/POOLED_BAM_FILES"/
cd ${BAM_DIR}

# # #This reads out the number of reads per pooled sample
# for BAM_FILE in $(find "$BAM_DIR" -name '*.bam')
# do
# 	OFPREFIX=$(basename "${BAM_FILE}" | sed 's/\.bam//g')    

# 	echo "${OFPREFIX}"
# 	samtools idxstats ${BAM_FILE} | awk '{s+=$3+$4} END {print s}'

# done


# #This sub-samples the pooled bam down to 30M (these fractions were manually calculated)
# samtools view -s 0.5153139806 -b Y_Astrocyte_pooled.bam > Y_Astrocyte_pooled_30M.bam 	#58216934	
# samtools view -s 0.2950093736 -b Y_qNSC_pooled.bam > Y_qNSC_pooled_30M.bam 				#101691684	
# samtools view -s 0.8454266925 -b Y_aNSC_pooled.bam > Y_aNSC_pooled_30M.bam 				#35485040
# samtools view -s 0.2782876265 -b Y_NPC_pooled.bam > Y_NPC_pooled_30M.bam  				#107802134
# samtools view -s 0.8610071062 -b Y_Endo_pooled.bam > Y_Endo_pooled_30M.bam 				#34842918
# samtools view -s 0.3172237055 -b O_Astrocyte_pooled.bam > O_Astrocyte_pooled_30M.bam 	#94570486
# samtools view -s 0.3974230875 -b O_qNSC_pooled.bam > O_qNSC_pooled_30M.bam 				#75486304
# samtools view -s 0.340898802 -b O_aNSC_pooled.bam > O_aNSC_pooled_30M.bam 				#88002656
# samtools view -s 0.4633622983 -b O_NPC_pooled.bam > O_NPC_pooled_30M.bam 				#64744154
# samtools view -s 0.4357718487 -b O_Endo_pooled.bam > O_Endo_pooled_30M.bam 				#68843364


# samtools view -s 0.537480158 -b Y_qNSC_invitro_pooled.bam > Y_qNSC_invitro_pooled_30M.bam 	#55816014
# samtools view -s 0.500254246 -b O_qNSC_invitro_pooled.bam > O_qNSC_invitro_pooled_30M.bam  #59969506
# samtools view -s 0.510061264 -b Y_aNSC_invitro_pooled.bam > Y_aNSC_invitro_pooled_30M.bam  #58816464
# samtools view -s 0.340464553 -b O_aNSC_invitro_pooled.bam > O_aNSC_invitro_pooled_30M.bam 	#88114900			


# # #Sort bam files
# for BAM_FILE in $(find "$BAM_DIR" -name '*_30M.bam')
# do
# 	OFPREFIX=$(basename "${BAM_FILE}" | sed 's/\_30M.bam//g')    
# 	samtools sort -o "${OFPREFIX}.sorted_30M.bam" ${BAM_FILE}
# done

# # #This reads out the number of reads per pooled sample
# for BAM_FILE in $(find "$BAM_DIR" -name '*.sorted_30M.bam')
# do
# 	OFPREFIX=$(basename "${BAM_FILE}" | sed 's/\.sorted_30M.bam//g')    

# 	echo "${OFPREFIX}"
# 	samtools idxstats ${BAM_FILE} | awk '{s+=$3+$4} END {print s}'

# done


# #Using deepTools alignmentSieve.py to perform Tn5 shifting
# for SORTED_BAM_FILE in $(find "$BAM_DIR" -name '*.sorted_30M.bam')
# do
# 	OFPREFIX=$(basename "${SORTED_BAM_FILE}" | sed 's/\.sorted_30M.bam//g')  
# 	samtools index ${SORTED_BAM_FILE}
# 	alignmentSieve -b ${SORTED_BAM_FILE} -o "${OFPREFIX}.tn5.sorted_30M.bam" --ATACshift
# done


# #Need to re-sort the bam files after Tn5 shifting
# for TRUE_SORTED_BAM_FILE in $(find "$BAM_DIR" -name '*.tn5.sorted_30M.bam')
# do
# 	OFPREFIX=$(basename "${TRUE_SORTED_BAM_FILE}" | sed 's/\.tn5.sorted_30M.bam//g')  
# 	echo "${OFPREFIX}"
# 	samtools sort -o "${OFPREFIX}.tn5.truesorted_30M.bam" ${TRUE_SORTED_BAM_FILE}
# done

#Indexing the file and checking to ensure than downsampling worked by reading out number of reads
#Also generates a bigWig signal track (with bin-size 1bp) for track visualization
for FINAL_BAM_FILE in $(find "$BAM_DIR" -name '*.tn5.truesorted_30M.bam')
do
	OFPREFIX=$(basename "${FINAL_BAM_FILE}" | sed 's/\.tn5.truesorted_30M.bam//g')  

	samtools index ${FINAL_BAM_FILE}

	echo "${OFPREFIX}"
	samtools idxstats ${FINAL_BAM_FILE} | awk '{s+=$3+$4} END {print s}'

	bamCoverage --bam ${FINAL_BAM_FILE} -o "${BAM_DIR}/${OFPREFIX}.tn5.truesorted_30M.bw" --binSize 1 --numberOfProcessors 6
done

