#!/bin/bash -l

echo "setting python configurations ..."

export PPD=$HOME/projects/apa_atingLab2019/01_polyAseq
export PYTHONPATH=$PPD:$PYTHONPATH
export PADD_GIT=$PPD/paddle-git

echo "running polyAseq pipelne ..."

echo "
python $PPD/polyAseq_lite.py \
	-C $PPD/configs/program.conf \
	-i ${PPD}/01_wkd/fastq \
	-e DKO \
	-c HCT \
	-d ${PPD}/01_wkd/comp_group.csv \
	-o ${PPD}/01_wkd/out \
	-b ${PPD}/01_wkd/barcode_list.csv \
	-p 6 \
	-r 0
"

python $PPD/polyAseq_lite.py \
	-C $PPD/configs/program.conf \
	-i ${PPD}/01_wkd/fastq \
	-e DKO \
	-c HCT \
	-d ${PPD}/01_wkd/comp_group.csv \
	-o ${PPD}/01_wkd/out \
	-b ${PPD}/01_wkd/barcode_list.csv \
	-p 6 \
	-r 0
