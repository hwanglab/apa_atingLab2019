#!/usr/bin/env bash

fastqd=fastq/FASTQ

echo "Check if fastq directory exists ..."
if [ ! -d $fastqd ]; then
	mkdir -p $fastqd
	echo "created ${fastqd}"
fi

rm -rf ${fastqd}/*.fastq.gz

echo "Start to download ChIP-seq/MBD-seq FASTQ files ..."

srrs=( SRR1234 SRR1235 SRR1236 SRR1237 ) #WARNING: modify this line!!!

for srr in "${srrs[@]}"
do
  echo "fastq-dump -F -O ${fastqd} $srr"
  fastq-dump -F -O ${fastqd} $srr
done

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
