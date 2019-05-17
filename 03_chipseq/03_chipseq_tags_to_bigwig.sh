#!/bin/bash -l

chipseq_tags=( "CTCF" "SMC1" "RAD21" "H3K27Ac" "Pol2Ser5" "Pol2Ser2" )

for chipseq_tag in "${chipseq_tags[@]}"
do
	echo "Rscript ./chipseq_tags_to_bigwig.r -n 8 -s ./sample_sheet.csv -p $chipseq_tag -f ftp://account:passwd@hostname/apa_atingLab2019/chipseq_bigwigs -o bigwig_tracks"
	
	Rscript ./chipseq_tags_to_bigwig.r -n 8 -s ./sample_sheet.csv -p $chipseq_tag -f ftp://account:passwd@hostname/apa_atingLab2019/chipseq_bigwigs -o bigwig_tracks
done


chipseq_tag="Mbd"
echo "Rscript ./chipseq_tags_to_bigwig.r -n 8 -s ./sample_sheet.csv -p $chipseq_tag -f ftp://account:passwd@hostname/apa_atingLab2019/chipseq_bigwigs -o bigwig_tracks -u 0"

Rscript ./chipseq_tags_to_bigwig.r -n 8 -s ./sample_sheet.csv -p $chipseq_tag -f ftp://account:passwd@hostname/apa_atingLab2019/chipseq_bigwigs -o bigwig_tracks -u 0
