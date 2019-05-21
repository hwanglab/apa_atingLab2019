#!/usr/bin/env bash

wkd_dir=./546goi_all

echo "1. Check if RD directory exists ..."
if [ ! -d $wkd_dir ]; then
	mkdir $wkd_dir
	echo "created ${wkd_dir}"
fi

echo "2. Downloading TCGA pAs usage and methylation beta value tables."

dnfile=github_Supplementary_Table_S6.tsv.xz

data_url="https://s3.amazonaws.com/apa2019/${dnfile}"

echo "wget $data_url -P ${wkd_dir}"

wget $data_url -P ${wkd_dir}

cd ${wkd_dir}

xz -d ${dnfile}
