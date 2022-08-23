#intersect_ChIP_ATAC.sh
#7/6/18
#The purpose of this script is to generate a subset of NSC ATAC-seq peaks that contain H3K27ac and p300 ChIP-seq peaks from Martynoga et al., 2013

ATAC_PATH="/Users/robinyeo//Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/Original_Data/FINAL_PEAK_FILES"
CHIP_PATH="/Users/robinyeo//Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/Original_Data/CHIP/Peak_Files"
OUT_PATH="/Users/robinyeo//Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/Output_Data/ChIP_Seq_Intersections/"

cd ~/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/Output_Data/ChIP_Seq_Intersections/

#Peak Numbers:
#Y_O_qNSC_ATAC 			64193
#Y_O_aNSC_ATAC 			58801
#qNSC_H3K27ac 			170993
#aNSC_H3K27ac 			94351
#qNSC_p300 				12953
#aNSC_p300 				9574
#qNSC_NFI 				38473



# # #the "-wa" option ensures that the peaks that get outputted are ATAC-seq peaks (that share an overlap with the ChIP-seq peaks)
bedtools intersect -wa -u -a "${ATAC_PATH}/Y_O_qNSC_peakset.txt" -b "${CHIP_PATH}/H3K27ac_quiNSC_rep1-pr.naive_overlap.filt.narrowPeak" > "${OUT_PATH}/qNSC_ATAC_peaks_with_H3K27ac.bed"
bedtools intersect -wa -u -a "${ATAC_PATH}/Y_O_aNSC_peakset.txt" -b "${CHIP_PATH}/H3K27ac_proNSC_rep1-pr.naive_overlap.filt.narrowPeak" > "${OUT_PATH}/aNSC_ATAC_peaks_with_H3K27ac.bed"

#qNSC_ATAC_peaks_with_H3K27ac 	23361/64193
#aNSC_ATAC_peaks_with_H3K27ac 	22392/58801


bedtools intersect -wa -u -a "${ATAC_PATH}/Y_O_qNSC_peakset.txt" -b "${CHIP_PATH}/p300_quiNSC_rep1-pr.naive_overlap.filt.regionPeak" > "${OUT_PATH}/qNSC_ATAC_peaks_with_p300.bed"
bedtools intersect -wa -u -a "${ATAC_PATH}/Y_O_aNSC_peakset.txt" -b "${CHIP_PATH}/p300_proNSC_rep1-pr.naive_overlap.filt.regionPeak" > "${OUT_PATH}/aNSC_ATAC_peaks_with_p300.bed"

#qNSC_ATAC_peaks_with_p300 	5401/64193
#aNSC_ATAC_peaks_with_p300 	3875/58801


bedtools intersect -wa -u -a "${OUT_PATH}/qNSC_ATAC_peaks_with_H3K27ac.bed" -b "${CHIP_PATH}/p300_quiNSC_rep1-pr.naive_overlap.filt.regionPeak" > "${OUT_PATH}/qNSC_ATAC_peaks_with_H3K27ac_p300.bed"
bedtools intersect -wa -u -a "${OUT_PATH}/aNSC_ATAC_peaks_with_H3K27ac.bed" -b "${CHIP_PATH}/p300_proNSC_rep1-pr.naive_overlap.filt.regionPeak" > "${OUT_PATH}/aNSC_ATAC_peaks_with_H3K27ac_p300.bed"

#qNSC_ATAC_peaks_with_H3K27ac_p300	5401/47424
#aNSC_ATAC_peaks_with_H3K27ac_p300 	3875/42829

