#!/bin/bash -l

script="./chipseq_pipeline_arg.sh"
ncpu=8

bash "$script" 01_0307_00CWCCF_HCT116_AHT_CTCF_hs_i72 37_030J_00CWCCF_HCT116_AHT_Input_hs_i84 HCT116_CTCF $ncpu TF
bash "$script" 13_031N_00CWCCF_HCT116_AHT_RAD21_hs_i72 37_030J_00CWCCF_HCT116_AHT_Input_hs_i84 HCT116_RAD21 $ncpu TF
bash "$script" 25_031J_00CWCCF_HCT116_AHT_Pol2Ser5_hs_i68 37_030J_00CWCCF_HCT116_AHT_Input_hs_i84 HCT116_Pol2Ser5 $ncpu histone
bash "$script" 07_030P_00CWCCF_HCT116_AHT_SMC1_hs_i86 37_030J_00CWCCF_HCT116_AHT_Input_hs_i84 HCT116_SMC1 $ncpu TF
bash "$script" 19_030B_00CWCCF_HCT116_AHT_H3K27Ac_hs_i79 37_030J_00CWCCF_HCT116_AHT_Input_hs_i84 HCT116_H3K27Ac $ncpu histone
bash "$script" 31_031F_00CWCCF_HCT116_AHT_Pol2Ser2_hs_i93 37_030J_00CWCCF_HCT116_AHT_Input_hs_i84 HCT116_Pol2Ser2 $ncpu histone

bash "$script" 02_0308_00CWCCF_DKO_AHT_CTCF_hs_i73 38_030K_00CWCCF_DKO_AHT_Input_hs_i85 DKO_CTCF $ncpu TF
bash "$script" 14_031O_00CWCCF_DKO_AHT_RAD21_hs_i73 38_030K_00CWCCF_DKO_AHT_Input_hs_i85 DKO_RAD21 $ncpu TF
bash "$script" 26_031K_00CWCCF_DKO_AHT_Pol2Ser5_hs_i69 38_030K_00CWCCF_DKO_AHT_Input_hs_i85 DKO_Pol2Ser5 $ncpu histone
bash "$script" 08_030Q_00CWCCF_DKO_AHT_SMC1_hs_i87 38_030K_00CWCCF_DKO_AHT_Input_hs_i85 DKO_SMC1 $ncpu TF
bash "$script" 20_030C_00CWCCF_DKO_AHT_H3K27Ac_hs_i80 38_030K_00CWCCF_DKO_AHT_Input_hs_i85 DKO_H3K27Ac $ncpu histone
bash "$script" 32_031G_00CWCCF_DKO_AHT_Pol2Ser2_hs_i94 38_030K_00CWCCF_DKO_AHT_Input_hs_i85 DKO_Pol2Ser2 $ncpu histone

bam_linkd="./fastq/kundaje_encode/bam_link"
if [ ! -d $bam_linkd ]; then
  mkdir -p $bam_linkd
fi

ln -s ./fastq/kundaje_encode/*_HCT116_*/align/rep1/*.nodup.bam ${bam_linkd}/
ln -s ./fastq/kundaje_encode/*_HCT116_*/align/rep1/*.nodup.bam.bai ${bam_linkd}/
ln -s ./fastq/kundaje_encode/*_DKO_*/align/rep1/*.nodup.bam ${bam_linkd}/
ln -s ./fastq/kundaje_encode/*_DKO_*/align/rep1/*.nodup.bam.bai ${bam_linkd}/
