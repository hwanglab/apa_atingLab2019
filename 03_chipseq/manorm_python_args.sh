#!/bin/bash -l

ctrl=$1
expr=$2
ctrl_cell=$3
expr_cell=$4
protein=$5
width=$6

wkd="fastq/kundaje_encode"

comp_tag=${ctrl_cell}_${expr_cell}_${protein}

outd=${wkd}/manorm/${comp_tag}
if [ ! -d $outd ]; then
	mkdir -p $outd
fi

zcat ${wkd}/${ctrl}/peak/macs2/rep1/${ctrl}*.filt.narrowPeak.gz | \
	cut -f1,2,3,10 > ${wkd}/${ctrl}/peak/macs2/rep1/${ctrl}.nodup.tagAlign.pval0.01.500K.filt.narrowPeak.4.xls

zcat ${wkd}/${expr}/peak/macs2/rep1/${expr}*.filt.narrowPeak.gz | \
	cut -f1,2,3,10 > ${wkd}/${expr}/peak/macs2/rep1/${expr}.nodup.tagAlign.pval0.01.500K.filt.narrowPeak.4.xls

gunzip -fc ${wkd}/${ctrl}/align/rep1/${ctrl}.nodup.tagAlign.gz > ${wkd}/${ctrl}/align/rep1/${ctrl}.nodup.tagAlign.bed

gunzip -fc ${wkd}/${expr}/align/rep1/${expr}.nodup.tagAlign.gz > ${wkd}/${expr}/align/rep1/${expr}.nodup.tagAlign.bed

echo "manorm --p2 ${wkd}/${ctrl}/peak/macs2/rep1/${ctrl}.nodup.tagAlign.pval0.01.500K.filt.narrowPeak.4.xls \
	--p1 ${wkd}/${expr}/peak/macs2/rep1/${expr}.nodup.tagAlign.pval0.01.500K.filt.narrowPeak.4.xls \
	--r2 ${wkd}/${ctrl}/align/rep1/${ctrl}.nodup.tagAlign.bed \
	--r1 ${wkd}/${expr}/align/rep1/${expr}.nodup.tagAlign.bed \
	-w $width -m 1 -p 0.01 \
	--name2 ${ctrl_cell}_${protein} \
	--name1 ${expr_cell}_${protein} \
	-o $outd"

manorm --p2 ${wkd}/${ctrl}/peak/macs2/rep1/${ctrl}.nodup.tagAlign.pval0.01.500K.filt.narrowPeak.4.xls \
	--p1 ${wkd}/${expr}/peak/macs2/rep1/${expr}.nodup.tagAlign.pval0.01.500K.filt.narrowPeak.4.xls \
	--r2 ${wkd}/${ctrl}/align/rep1/${ctrl}.nodup.tagAlign.bed \
	--r1 ${wkd}/${expr}/align/rep1/${expr}.nodup.tagAlign.bed \
	-w $width -m 1 -p 0.01 \
	--name2 ${ctrl_cell}_${protein} \
	--name1 ${expr_cell}_${protein} \
	-o $outd
