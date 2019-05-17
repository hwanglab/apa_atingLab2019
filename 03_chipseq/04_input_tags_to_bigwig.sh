#!/bin/bash -l

bamCvgBin="bamCoverage"
workD="fastq/kundaje_encode"
outD="bigwig_tracks"

bam=${workD}/01_0307_00CWCCF_HCT116_AHT_CTCF_hs_i72/align/ctl1/37_030J_00CWCCF_HCT116_AHT_Input_hs_i84.nodup.bam
outFileName=${outD}/37_030J_00CWCCF_HCT116_AHT_Input_hs_i84.nodup.bw
$bamCvgBin --bam $bam --binSize 5 --outFileFormat bigwig --outFileName $outFileName

bam=${workD}/02_0308_00CWCCF_DKO_AHT_CTCF_hs_i73/align/ctl1/38_030K_00CWCCF_DKO_AHT_Input_hs_i85.nodup.bam
outFileName=${outD}/38_030K_00CWCCF_DKO_AHT_Input_hs_i85.nodup.bw
$bamCvgBin --bam $bam --binSize 5 --outFileFormat bigwig --outFileName $outFileName

bam=${workD}/03_0309_00CWCCF_DU145_AHT_CTCF_hs_i75/align/ctl1/39_030L_00CWCCF_DU145_AHT_Input_hs_i86.nodup.bam
outFileName=${outD}/39_030L_00CWCCF_DU145_AHT_Input_hs_i86.nodup.bw
$bamCvgBin --bam $bam --binSize 5 --outFileFormat bigwig --outFileName $outFileName

bam=${workD}/04_030A_00CWCCF_LNCaP_AHT_CTCF_hs_i78/align/ctl1/40_030M_00CWCCF_LNCaP_AHT_Input_hs_i88.nodup.bam
outFileName=${outD}/40_030M_00CWCCF_LNCaP_AHT_Input_hs_i88.nodup.bw
$bamCvgBin --bam $bam --binSize 5 --outFileFormat bigwig --outFileName $outFileName

