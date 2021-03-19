#1/18/2019

#I had initally used the pooled tagAlign files which I subsampled and converted to bam files but this caused information loss and wouldn't allow me to run nucleoATAC
#So I decided to start from scratch and manually merge the high-quality (filtered,de-duplicated, NON-TN5-SHIFTED) bam files and then subsample them directly
#I then sorted and Tn5-shifted them

BAM_DIR="/Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/POOLED_BAM"
cd ${BAM_DIR}

#This manually merges the high-quality de-duplicated bam files
#samtools merge Y_Astrocyte_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-1-TAAGGCGA_S1_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-2_S2_R1_001.trim.PE2SE.nodup.bam
#samtools merge Y_qNSC_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-2-CGTACTAG_S2_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-2-CGTACTAG_S2_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-3_S3_R1_001.trim.PE2SE.nodup.bam
#samtools merge Y_aNSC_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-3-AGGCAGAA_S3_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-4_S4_R1_001.trim.PE2SE.nodup.bam
#samtools merge Y_NPC_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-4-TCCTGAGC_S4_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-4-TCCTGAGC_S4_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-5_S5_R1_001.trim.PE2SE.nodup.bam
#samtools merge Y_Endo_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-5-GGACTCCT_S5_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-5-GGACTCCT_S5_R1_001.trim.PE2SE.nodup.bam
#samtools merge O_Astrocyte_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-6-TAGGCATG_S6_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-6-TAGGCATG_S6_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-7_S7_R1_001.trim.PE2SE.nodup.bam
#samtools merge O_qNSC_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-7-CTCTCTAC_S7_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-7-CTCTCTAC_S7_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-8_S8_R1_001.trim.PE2SE.nodup.bam
#samtools merge O_aNSC_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-8-CAGAGAGG_S8_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-8-CAGAGAGG_S8_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-9_S9_R1_001.trim.PE2SE.nodup.bam
#samtools merge O_NPC_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-9-GCTACGCT_S9_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP3/RY-10_S10_R1_001.trim.PE2SE.nodup.bam
#samtools merge O_Endo_pooled.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP1/RY-YoungOld-10-CGAGGCTG_S10_R1_001.trim.PE2SE.nodup.bam /Volumes/Drive6TB/YOUNG_OLD_ATAC_ANALYSIS_2018/BAM_FILES_ALL_REPS/REP2/RY-YO-Rep-10-CGAGGCTG_S10_R1_001.trim.PE2SE.nodup.bam



# #This reads out the number of reads per pooled sample
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


# #Sort bam files
# for BAM_FILE in $(find "$BAM_DIR" -name '*.bam')
# do
# 	OFPREFIX=$(basename "${BAM_FILE}" | sed 's/\.bam//g')    
# 	samtools sort -o "${OFPREFIX}.sorted.bam" ${BAM_FILE}
# done


# #Using deepTools alignmentSieve.py to perform Tn5 shifting
# for SORTED_BAM_FILE in $(find "$BAM_DIR" -name '*.sorted.bam')
# do
# 	OFPREFIX=$(basename "${SORTED_BAM_FILE}" | sed 's/\.sorted.bam//g')  
# 	samtools index ${SORTED_BAM_FILE}
# 	alignmentSieve -b ${SORTED_BAM_FILE} -o "${OFPREFIX}.tn5.sorted.bam" --ATACshift
# done


# #Need to re-sort the bam files after Tn5 shifting
# for TRUE_SORTED_BAM_FILE in $(find "$BAM_DIR" -name '*.tn5.sorted.bam')
# do
# 	OFPREFIX=$(basename "${TRUE_SORTED_BAM_FILE}" | sed 's/\.tn5.sorted.bam//g')  
# 	echo "${OFPREFIX}"
# 	samtools sort -o "${OFPREFIX}.tn5.truesorted.bam" ${TRUE_SORTED_BAM_FILE}
# done

# #Indexing the file and checking to ensure than downsampling worked by reading out number of reads
# #Also generates a bigWig signal track (with bin-size 1bp) for track visualization
# for FINAL_BAM_FILE in $(find "$BAM_DIR" -name '*.tn5.truesorted.bam')
# do
# 	OFPREFIX=$(basename "${FINAL_BAM_FILE}" | sed 's/\.tn5.truesorted.bam//g')  

# 	samtools index ${FINAL_BAM_FILE}

# 	echo "${OFPREFIX}"
# 	samtools idxstats ${FINAL_BAM_FILE} | awk '{s+=$3+$4} END {print s}'

# 	bamCoverage --bam ${FINAL_BAM_FILE} -o "${BAM_DIR}/${OFPREFIX}.tn5.truesorted.bw" --binSize 1 --numberOfProcessors 6
# done

