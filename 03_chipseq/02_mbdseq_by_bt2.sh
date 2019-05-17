#!/bin/bash -l

ncpu=8
fastqd=fastq
bamd=${fastqd}/bt2

if [ ! -d $bamd ]; then
	mkdir -p $bamd
fi

bt2_idx="../resource/ref/hg19_PhiX"

read1=SRR030224
read2=SRR030225
read_out=HCT116

bowtie2 -p 8 --rg-id $read_out --rg $read1 --rg $read2 -x $bt2_idx -U $fastqd/${read1}.fastq.gz,$fastqd/${read2}.fastq.gz | samtools view -u - | samtools sort - -o ${bamd}/${read_out}_Mbd.bam
samtools index ${bamd}/${read_out}_Mbd.bam
samtools flagstat -@ 8 ${bamd}/${read_out}_Mbd.bam > ${bamd}/${read_out}_Mbd.bam.stats

read1=SRR030220
read_out=DKO

bowtie2 -p 8 --rg-id $read_out --rg $read1 -x $bt2_idx -U $fastqd/${read1}.fastq.gz | samtools view -u - | samtools sort - -o ${bamd}/${read_out}_Mbd.bam
samtools index ${bamd}/${read_out}_Mbd.bam
samtools flagstat -@ 8 ${bamd}/${read_out}_Mbd.bam > ${bamd}/${read_out}_Mbd.bam.stats

deseq2_outd=${bamd}/HCT116_DKO_deseq2

if [ ! -d $deseq2_outd ]; then
	mkdir -p $deseq2_outd
fi

Rscript deseq2_log2fc_diffcall.r -c ${bamd}/HCT116_Mbd.bam -e ${bamd}/DKO_Mbd.bam -t HCT116_DKO -o ${deseq2_outd}

# ----------------------------
# The following commandlines are for PrEC, DU145, and LNCaP which are not discussed in the manuscript submssion!

# read1=SRR402877
# read2=SRR402878
# read_out=PrEC
# 
# bowtie2 -p 8 --rg-id $read_out --rg $read1 --rg $read2 -x $bt2_idx -U $fastqd/${read1}.fastq.gz,$fastqd/${read2}.fastq.gz | samtools view -u - | samtools sort - -o ${bamd}/${read_out}.bam
# samtools index ${bamd}/${read_out}.bam
# samtools flagstat -@ 8 ${bamd}/${read_out}.bam > ${bamd}/${read_out}.bam.stats
# 
# read1=SRR402881
# read2=SRR402882
# read_out=DU145
# 
# bowtie2 -p 8 --rg-id $read_out --rg $read1 --rg $read2 -x $bt2_idx -U $fastqd/${read1}.fastq.gz,$fastqd/${read2}.fastq.gz | samtools view -u - | samtools sort - -o ${bamd}/${read_out}.bam
# samtools index ${bamd}/${read_out}.bam
# samtools flagstat -@ 8 ${bamd}/${read_out}.bam > ${bamd}/${read_out}.bam.stats
# 
# read1=SRR402879
# read2=SRR402880
# read_out=LNCaP
# 
# bowtie2 -p 8 --rg-id $read_out --rg $read1 --rg $read2 -x $bt2_idx -U $fastqd/${read1}.fastq.gz,$fastqd/${r