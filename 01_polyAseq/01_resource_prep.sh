#!/usr/bin/env bash

wkd=../resource
refd=${wkd}/ref

echo "1. Check if fastq directory exists ..."
if [ ! -d $refd ]; then
	mkdir $refd
	echo "created ${refd}"
fi

hg19_2bit_url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"

echo "2. download ucsc hg19 reference file ..."
twoBitToFa $hg19_2bit_url ${refd}/hg19.fa

echo "3. Merging hg19 + phix control ..."
cat ${refd}/hg19.fa ${refd}/PhiX.fa > ${refd}/hg19_PhiX.fa

echo "4. building ref index ..."
samtools faidx ${refd}/hg19_PhiX.fa

echo "5. generating chromosome size ..."
cut -f1,2 ${refd}/hg19_PhiX.fa.fai > ${refd}/ChromSizes.hg19.txt

echo "6. downloading ensembl hg19/b37 gtf file ..."
ensembl_gtf_url="http://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
rm -rf ${refd}/ensembl/*.gtf.gz*
wget -P ${refd}/ensembl $ensembl_gtf_url
gunzip ${refd}/ensembl/Homo_sapiens.GRCh37.87.gtf.gz
echo "Done."

echo "7. adding chr to the contig headers in b37 gtf ..."
grep -P '^[0-9]+' ${refd}/ensembl/Homo_sapiens.GRCh37.87.gtf | sed 's/^/chr/' > ${refd}/ensembl/Homo_sapiens.GRCh37.87.chr.gtf
sed 's/^MT/chrM' ${refd}/ensembl/Homo_sapiens.GRCh37.87.gtf >> ${refd}/ensembl/Homo_sapiens.GRCh37.87.chr.gtf

echo "8. building bowtie2 index ..."
bowtie2 ${refd}/hg19_PhiX.fa ${refd}/hg19_PhiX

echo "9. Download 01_polyAseq.data ..."
apa_atingLab2019_data_url="http://apa_atingLab2019_data_url"

polyAseq_data_url="https://s3.amazonaws.com/apa2019/01_polyAseq.data.tar.gz"

wget -P $polyAseq_data_url
tar zxvf ./01_polyAseq.data.tar.gz
