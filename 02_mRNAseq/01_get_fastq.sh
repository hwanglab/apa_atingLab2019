#!/usr/bin/env bash

fastqd=fastq

echo "Check if fastq directory exists ..."
if [ ! -d $fastqd ]; then
	mkdir -p $fastqd
	echo "created ${fastqd}"
fi

rm -rf ${fastqd}/*.fastq.gz

echo "Start to download mRNA-seq FASTQ files ..."

echo "Contact to authors to obtain SRR IDs!"

srrs=( SRR1234 SRR1235 )

for srr in "${srrs[@]}"
do
  echo "fastq-dump -F -O ${fastqd} $srr"
  fastq-dump -F --split-files -O ${fastqd} $srr
done

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
