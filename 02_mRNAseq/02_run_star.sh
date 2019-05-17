#!/bin/bash -l

ncpu=8

ref_prefix="../resource/ref/star"

if [ ! -d $ref_prefix ]; then
	mkdir -p $ref_prefix
fi

echo "check if a STAR index file exists ..."

sample_idx_file="${ref_prefix}/SAindex"

if [ -f $sample_idx_file ]; then
	echo "star index files are already built."
else
	echo "star index files do not exist!"
	
	ref_fasta="../resource/ref/hg19_PhiX.fa"
	if [ ! -f $ref_fasta ]; then
		echo "reference FASTA file [$ref_fasta] does not exist!"
		exit(1)
	fi
	
	gtf_file="../resource/ensembl/Homo_sapiens.GRCh37.87.chr.gtf"
	if [ ! -f $gtf_file ]; then
		echo "GTF annotation [$gtf_file]  file does not exist!"
		exit(1)
	fi
	
	echo "building a reference genome index ..."
	echo "STAR --runThreadN $ncpu --runMode genomeGenerate --genomeDir $ref_prefix --genomeFastaFiles $ref_fasta --sjdbGTFfile $gtf_file"
	STAR --runThreadN $ncpu --runMode genomeGenerate --genomeDir $ref_prefix --genomeFastaFiles $ref_fasta --sjdbGTFfile $gtf_file
fi

fastqDir="fastq"
if [ ! -d $fastqDir ]; then
	echo "Run 01_get_fastq.sh first!"
	exit(1)
fi

bamDir="${fastqDir}/star"

if [ ! -d $bamDir ]; then
	mkdir -p $bamDir
fi

readPrefix=("HI.1159.005.Index_13.HCT116" "HI.1159.005.Index_16.DKO")

for read in "${readPrefix[@]}"
do
	out_bam_prefix="${bamDir}/${read}_"
	echo "STAR --genomeDir $ref_prefix \
					--runThreadN $ncpu \
					--readFilesIn ${fastqDir}/${read}_R1.fastq.gz ${fastqDir}/${read}_R2.fastq.gz \
					--readFilesCommand zcat \
					--outSAMtype BAM SortedByCoordinate \
					--outSAMattributes Standard \
					--outSAMunmapped Within \
					--outSAMattrRGline ID:4 SM:20 LB:${read}_mRNA PU:unit PL:illumina \
					--twopassMode Basic \
					--outFileNamePrefix $out_bam_prefix"
					
	out_bam="${out_bam_prefix}Aligned.sortedByCoord.out.bam"
		
	echo "samtools index $out_bam"
	
	samtools index $out_bam

	echo "Done ${read}."
done
echo "Done."
