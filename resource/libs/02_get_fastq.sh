#!/usr/bin/env bash

## TODO
polyaseq_fastq_url="ftp://apa_atingLab2019/fastq/polyseq"

wkd=01_wkd
fastqd=${wkd}/fastq

echo "Check if fastq directory exists ..."
if [ ! -d $fastqd ]; then
	mkdir -p $fastqd
	echo "created ${fastqd}"
fi

rm -rf ${fastqd}/*.fastq.gz

echo "Start to download PolyA-seq FASTQ files ..."
wget -r --no-parent -A '*.tar.gz' $polyaseq_fastq_url

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
