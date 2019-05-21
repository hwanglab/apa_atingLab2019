#!/bin/bash -l

rbin=./02b_tcga_apa_paur_plot.r

input_xlsx=../resource/pau_testing_region/apa546_2018_11_16.xlsx
topk=8000

prepropa_rd=../01_polyAseq/01_wkd/out/03_CallApa/output/prepropa.rd
ann_rd=../01_polyAseq/01_wkd/out/04_AnnoApa/output/apa.ann.rd

in_rd_pref=./546goi_all
outd=./546goi_all

if [ ! -d $outd ]; then
  mkdir -p $outd;
fi

echo "
Rscript $rbin \
	-n 2 \
	-t $in_rd_pref \
	-T $topk \
	-b ${input_xlsx} \
	-rd1 ${prepropa_rd} \
	-rd2 ${ann_rd} \
	-o $outd \
	-s gtop${topk} \
	-r 4 \
	-d 0
"

Rscript $rbin \
	-n 2 \
	-t $in_rd_pref \
	-T $topk \
	-b ${input_xlsx} \
	-rd1 ${prepropa_rd} \
	-rd2 ${ann_rd} \
	-o $outd \
	-s gtop${topk} \
	-r 4 \
	-d 0
