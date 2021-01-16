#!/bin/bash -l
echo "Make sure that fastq-dump in $PATH and configure it appropriately"

wkd=01_wkd
fastqd=${wkd}/fastq

echo "Check if fastq directory exists ..."
if [ ! -d $fastqd ]; then
	mkdir -p $fastqd
	echo "created ${fastqd}"
fi

rm -rf ${fastqd}/*.fastq.gz

echo "Start to download PolyA-seq FASTQ files ..."

#	- SRP083252 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86178)
#	- SRP083254 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86180)

srrs=(SRR4091084 SRR4091085 SRR4091086 SRR4091087 SRR4091088 SRR4091089 SRR4091090 SRR4091091 SRR4091104 SRR4091105 SRR4091106 SRR4091107 SRR4091108 SRR4091109 SRR4091110 SRR4091111 SRR4091113 SRR4091115 SRR4091117 SRR4091119)

for srr in "${srrs[@]}"
do
  echo "fastq-dump -F -O ${fastqd} $srr"
  fastq-dump -F -O ${fastqd} $srr
done

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
