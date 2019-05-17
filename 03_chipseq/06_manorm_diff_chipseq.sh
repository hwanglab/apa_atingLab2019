#!/bin/bash -l

# script="/mnt/isilon/data/QHS/home/hongc2/projects/apa/activemotif_fastq/encode_analy_out/scripts/diff_tests/manorm_python_args.sh"
script="manorm_python_args.sh"

bash "$script" 01_0307_00CWCCF_HCT116_AHT_CTCF_hs_i72 02_0308_00CWCCF_DKO_AHT_CTCF_hs_i73 HCT116 DKO CTCF 500
bash "$script" 07_030P_00CWCCF_HCT116_AHT_SMC1_hs_i86 08_030Q_00CWCCF_DKO_AHT_SMC1_hs_i87 HCT116 DKO SMC1 500
bash "$script" 13_031N_00CWCCF_HCT116_AHT_RAD21_hs_i72 14_031O_00CWCCF_DKO_AHT_RAD21_hs_i73 HCT116 DKO RAD21 500
bash "$script" 19_030B_00CWCCF_HCT116_AHT_H3K27Ac_hs_i79 20_030C_00CWCCF_DKO_AHT_H3K27Ac_hs_i80 HCT116 DKO H3K27Ac 1000
bash "$script" 25_031J_00CWCCF_HCT116_AHT_Pol2Ser5_hs_i68 26_031K_00CWCCF_DKO_AHT_Pol2Ser5_hs_i69 HCT116 DKO Pol2Ser5 1000
bash "$script" 31_031F_00CWCCF_HCT116_AHT_Pol2Ser2_hs_i93 32_031G_00CWCCF_DKO_AHT_Pol2Ser2_hs_i94 HCT116 DKO Pol2Ser2 1000
