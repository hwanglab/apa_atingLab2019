library(argparse)
library(data.table)

debug <- FALSE
if (!debug) {
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  source_root <- dirname(dirname(scriptPath))
  parser <- ArgumentParser(description='macs2_deseq2')
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                      dest="ncpu", default = 1,
                      help="number of cpus to utilize [1]")
  
  parser$add_argument("-q", "--min_mapq", type="integer", required=FALSE,
                      dest="min_mapq", default = 1,
                      help="min mapq to be applied to bams [1]")
  
  parser$add_argument("-g", "--genome", type="character", default='hg19', required=FALSE,
                      dest="genome",
                      help="genome[hg19]")
  
  parser$add_argument("-c", "--bam_co", type="character", required=TRUE,
                      dest="bam_co",
                      help="bam control")
  
  parser$add_argument("-t", "--bam_tr", type="character", required=TRUE,
                      dest="bam_tr",
                      help="bam treatment")
  
  parser$add_argument("-p", "--peak_type", type="character", default='sharp', required=FALSE,
                      dest="peak_type",
                      help="[sharp],broad")
  parser$add_argument("-m", "--peak_merge_method", type="character", default='concat', required=FALSE,
                      dest="peak_merge_method",
                      help="[concat],merge")
  
  parser$add_argument("-o", "--out_dir", type="character", required=TRUE,
                      dest="out_dir",
                      help="output directory")
  
  args <- parser$parse_args()

} else {
args <- data.table(bam_co="~/projects/apa_atingLab2019/03_chipseq/fastq/bt2/HCT116_Mbd.bam",
										 bam_tr="~/projects/apa_atingLab2019/03_chipseq/fastq/bt2/DKO_Mbd.bam",
										 peak_type="broad",
										 peak_merge_method="merged",
										 min_mapq=0,
										 out_dir="~/projects/apa_atingLab2019/03_chipseq/fastq/bt2/HCT116_DKO_deseq2")
}

library(parallel)
library(Rsubread)
source(file.path(Sys.getenv('R_UTIL_APA'),'lib_apps.R'))

dir.create(args$out_dir, showWarnings = FALSE)

message('bam --> MACS2 --> peaks')

samp <- data.table(bam=c(args$bam_co,args$bam_tr),
                   name=c('control','experiment'),
                   group=c('control','experiment'))

qc_bam <- function(bam,min_mapq=1){
  qc_bam <- sprintf('%s.q%d.bam',bam,min_mapq)
  cmd <- sprintf('samtools view -q %d -hb -o %s %s',min_mapq,qc_bam,bam)
  message(cmd)
  system(cmd)
  return(qc_bam)
}
  
if (args$min_mapq>0){
  upd_bams <- mapply(qc_bam,
                     samp$bam,
                     MoreArgs = list(min_mapq=args$min_mapq))
  samp$bam <- upd_bams
}

samp$peak_tsv <- mapply(macs2_peak,
                        samp$bam,
                        samp$name,
                        file.path(args$out_dir,samp$name),
                        MoreArgs = list(peak_type=args$peak_type,
                                        out_prefix=args$peak_type))
if (args$peak_merge_method=='concat'){
	message('narrow_peak1, narrow_peak2 --> concatenating --> unified peaks')
	peak12 <- macs2_narrowpeak_concat(samp$peak_tsv)

} else {
	message('narrow_peak1, narrow_peak2 --> merging --> unified peaks w/o overlap')
	peak12 <- macs2_narrowpeak_merge(samp$peak_tsv)
}
message('preparing deseq2 input')
deseq2_input <- prep_deseq2(samp,peak12,ncore=8)

message('feature matrix --> deseq2 --> log2fc table')
out_file <- run_deseq2(deseq2_input,args$out_dir)
