#!/bin/bash -l
readbase=$1
ctrlbase=$2
title=$3
ncpu=$4
type=$5

fastqd="./fastq"
outbase=${fastqd}/kundaje_encode
blacklist_fn="../resource/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

unset PATH
unset PYTHONPATH
unset R_LIBS

# modify the following settings match with your configuration
# -------------------------------------------------------------------->
progd="$HOME/apps/chipseqs/TF_chipseq_pipeline"
species_conf="/opt/bds_pipeline_genome_data/aquas_chipseq_species.conf"

export PATH=$HOME/.local/bin
export PATH=$HOME/.bds:$PATH
export PATH=$HOME/apps/miniconda3/bin:$PATH
export PATH=/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:$PATH
# <--------------------------------------------------------------------

export _JAVA_OPTIONS="-Xms256M -Xmx728M -XX:ParallelGCThreads=1"

outd=${outbase}/${readbase}
if [ ! -d $outd ]; then
	mkdir -p $outd
fi

echo "python \
${progd}/chipseq.py \
--se \
--species hg19 \
--species-file ${species_conf} \
--peak-caller macs2 \
--blacklist ${blacklist_fn} \
--out-dir ${outd} \
--title ${title} \
--nth ${ncpu} \
--fastq1 ${fastqd}/FASTQ/${readbase}.fastq.gz \
--ctl_fastq1 ${fastqd}/FASTQ/${ctrlbase}.fastq.gz \
--type $type
"

python ${progd}/chipseq.py \
--se \
--species hg19 \
--species-file ${species_conf} \
--peak-caller macs2 \
--blacklist ${blacklist_fn} \
--out-dir ${outd} \
--title ${title} \
--nth ${ncpu} \
--fastq1 ${fastqd}/FASTQ/${readbase}.fastq.gz \
--ctl_fastq1 ${fastqd}/FASTQ/${ctrlbase}.fastq.gz \
--type $type
