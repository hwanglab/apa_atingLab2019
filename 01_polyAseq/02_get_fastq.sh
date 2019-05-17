#!/usr/bin/env bash

wkd=01_wkd
fastqd=${wkd}/fastq

echo "Check if fastq directory exists ..."
if [ ! -d $fastqd ]; then
	mkdir -p $fastqd
	echo "created ${fastqd}"
fi

rm -rf ${fastqd}/*.fastq.gz

echo "Start to download PolyA-seq FASTQ files ..."

srrs=( SRR1234 SRR1235 SRR1236 SRR1236 ) #WARNING: modify this line!!!

for srr in "${srrs[@]}"
do
  echo "fastq-dump -F -O ${fastqd} $srr"
  fastq-dump -F -O ${fastqd} $srr
done

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
