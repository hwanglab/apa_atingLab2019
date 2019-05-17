#!/bin/bash -l

#SBATCH -J tcgaplot
#SBATCH -o /mnt/isilon/data/QHS/home/hongc2/projects/apa/src/s11_tcga/scripts/logs/04_tcga_analy_apa_paur_plot.%N.%j.out
#SBATCH -e /mnt/isilon/data/QHS/home/hongc2/projects/apa/src/s11_tcga/scripts/logs/04_tcga_analy_apa_paur_plot.%N.%j.err
#SBATCH -p TingA
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hongc2@ccf.org

export LABD=/mnt/isilon/data/w_QHS/hwangt-share
export APP=$LABD/apps
export R_LIBS=$APP/R/lib

module load htslib/current
module load pcre/8.38
module load R/3.4.3

mountd=/mnt/isilon/data/QHS/home/hongc2
rbin=${mountd}/projects/apa/src/s11_tcga/04_tcga_analy_apa_paur_plot.r

which Rscript
which R
Rscript --version

input_xlsx=apa546_2018_11_16.xlsx
topk=8000

in_rd_pref=${mountd}/projects/apa/src/s11_tcga/output/03_tcga_analy_apa/546goi
outd=${mountd}/projects/apa/src/s11_tcga/output/03_tcga_analy_apa/546goi_all

if [ ! -d $outd ]; then
  mkdir -p $outd;
fi

echo "
Rscript $rbin \
	-n 2 \
	-t $in_rd_pref \
	-T $topk \
	-b ${mountd}/projects/apa/src/s11_tcga/input/${input_xlsx} \
	-o $outd \
	-s gtop${topk} \
	-r 4 \
	-d 0
"

Rscript $rbin \
	-n 2 \
	-t $in_rd_pref \
	-T $topk \
	-b ${mountd}/projects/apa/src/s11_tcga/input/${input_xlsx} \
	-o $outd \
	-s gtop${topk} \
	-r 4 \
	-d 0
