#!/usr/bin/env bash

fastqd=fastq/FASTQ

echo "Check if fastq directory exists ..."
if [ ! -d $fastqd ]; then
	mkdir -p $fastqd
	echo "created ${fastqd}"
fi

rm -rf ${fastqd}/*.fastq.gz

echo "Start to download ChIP-seq/MBD-seq FASTQ files ..."

# https://www.ncbi.nlm.nih.gov/sra/SRP001414

srrs=( SRR030224 SRR030225 SRR030220 SRR_Chipseq1 SRR_Chipseq2 )

for srr in "${srrs[@]}"
do
  echo "fastq-dump -F -O ${fastqd} $srr"
  fastq-dump -F -O ${fastqd} $srr
done

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
