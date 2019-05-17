#!/bin/bash -l

ncpu=$1
outd="star/output"

if [ ! -d $outd ]; then
	mkdir -p $outd
fi

Rscript ./04_apafactorexp.r \
  -i ../01_polyAseq/01_wkd/out/03_CallApa/output/prepropa.rd \
  -a ../01_polyAseq/data/apa_factors.csv \
  -c HCT116 \
  -e DKO \
  -cb star/HI.1159.005.Index_13.HCT116_Aligned.sortedByCoord.out.bam \
  -eb star/HI.1159.005.Index_16.DKO_Aligned.sortedByCoord.out.bam \
  -g ../resource/ensembl/Homo_sapiens.GRCh37.87.chr.gtf \
  -n $ncpu \
  -o $outd
