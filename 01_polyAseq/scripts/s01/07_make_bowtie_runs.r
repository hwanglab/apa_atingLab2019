library(data.table)
library(stringr)
library(argparse)

parser <- ArgumentParser(description='a_chop_plot')

parser$add_argument("-f", "--sample_sheet", type="character", required=TRUE,
										dest="sample_sheet",
										help="fastq_index.csv")

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 1,
										help="number of cpus to utilize [1]")

parser$add_argument("-r", "--ref_prefix", type="character", required=TRUE,
										dest="ref_prefix",
										help="target reference genome index prefix")

parser$add_argument("-i", "--qc_fastq_dir", type="character", required=TRUE,
										dest="qc_fastq_dir",
										help="polyA trimmed FASTQ file directory")

parser$add_argument("-o", "--script_dir", type="character", required=TRUE,
										dest="script_dir",
										help="script output directory")
args <- parser$parse_args()

script_dir1 = sprintf("%s/",args$script_dir)
bamdir = dirname(args$script_dir)
bamdir1 = sprintf("%s/",bamdir)
logdir = file.path(bamdir,'log')
dir.create(logdir, showWarnings = FALSE)
logdir1 = sprintf("%s/",logdir)

#fastq <- dir("output/chop_a_tails",pattern="_A.fastq.gz",full.names=T)
fastq <- dir(args$qc_fastq_dir,pattern="_A.fastq.gz",full.names=T)

# Want to start using sample names at this point, so don't want to bother aligning the OTHER and ones where we didn't expect to see that barcode
#index <- fread("./input/fastq_index.csv")
index <- fread(args$sample_sheet)
index$file <- str_replace(str_replace(index$file,"fastq/",""),"\\.fastq\\.gz","")

#index$key <- paste0("output/chop_a_tails/",index$file,"_",index$barcode,"_A.fastq.gz")
index$key <- paste0(args$qc_fastq_dir,"/",index$file,"_",index$barcode,"_A.fastq.gz")

stopifnot(index$key %in% fastq)
bt2cmd = sprintf("bowtie2 -p %d -x %s -U ",args$ncpu,args$ref_prefix)

cmds <- paste0(bt2cmd,index$key," 2> ",logdir1,index$sample,".log.txt",sprintf(" | samtools view -b - 1> %s",bamdir1),index$sample,".bam")

#cmds <- lapply(cmd,function(x) c("module add bowtie2/2.2.5",x))

for(i in 1:length(cmds))
{
	#writeLines(cmds[[i]],paste0("bowtie_alignments/run_scripts/",i,".sh"))
	writeLines(cmds[i],paste0(script_dir1,i,".sh"))
}
