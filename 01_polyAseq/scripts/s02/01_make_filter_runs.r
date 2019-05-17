library(data.table)
library(stringr)
library(argparse)

parser <- ArgumentParser(description='make_filter_runs')

parser$add_argument("-i", "--bam_dir", type="character", required=TRUE,
										dest="bam_dir",
										help="bowtie BAM file directory")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output directory")

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 1,
										help="number of cpus to utilize [1]")

args <- parser$parse_args()

bams <- dir(args$bam_dir,pattern=".bam",full.names=T)
samples <- str_replace(basename(bams),"\\.bam","")

#cmd <- paste0(sprintf("samtools sort -@ %d -o ",args$ncpu),paste0(args$output_dir,samples,".01_sort.bam")," -T ",samples," ",bams)
prep_sorted_bams <- file.path(args$output_dir,sprintf("%s.01_sort.bam",samples))
cmd <- paste0(sprintf("samtools sort -@ %d -o ",args$ncpu),prep_sorted_bams," -T ",file.path('/tmp',samples)," ",bams)

#dir.create("output/run_scripts_sort",recursive=TRUE)
dir.create(sprintf("%s/run_scripts_sort",args$output_dir),recursive=TRUE)
for(i in 1:length(cmd))
{
	#writeLines(cmd[[i]],paste0("output/run_scripts_sort/",i,".sh"))
	writeLines(cmd[[i]],paste0(sprintf("%s/run_scripts_sort/",args$output_dir),i,".sh"))
}

cmd <- paste0(sprintf("samtools view -q 20 -b -o %s/",args$output_dir),paste0(samples,".02_mapq20.bam")," ",prep_sorted_bams)

#dir.create("output/run_scripts_filter",recursive=TRUE)
dir.create(sprintf("%s/run_scripts_filter",args$output_dir),recursive=TRUE)
for(i in 1:length(cmd))
{
	#writeLines(cmd[[i]],paste0("output/run_scripts_filter/",i,".sh"))
	writeLines(cmd[[i]],paste0(sprintf("%s/run_scripts_filter/",args$output_dir),i,".sh"))
}
