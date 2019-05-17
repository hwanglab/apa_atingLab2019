library(argparse)
library(data.table)

if (T) {
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  
  parser <- ArgumentParser(description='prepropa')
  
  parser$add_argument("-i", "--encode_out_base_dir", type="character", required=TRUE,
                      dest="encode_out_base_dir",
                      help="encode_out_base_dir")
  
  parser$add_argument("-o", "--out_dir", type="character", required=TRUE,
                      dest="out_dir",
                      help="output directory")
  
  args <- parser$parse_args()
} else {
	args <- data.table(encode_out_base_dir='fastq/kundaje_encode',
										 out_dir='fastq/kundaje_encode/alignment_stats')
}


source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
rep1i <- c(1,3,5)
# ctl1i <- c(2,4,6)

jsons <- Sys.glob(sprintf('%s/*00CWCCF_*/ENCODE_summary.json',args$encode_out_base_dir))

N <- length(strsplit(jsons[1],'/')[[1]]) - 1

out_prefix<-c('raw','filt','dup')
sample_base <- sapply(jsons,function(x){strsplit(x,'/')[[1]][N]})

if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir)
}

for (i in 1:3) {
  tabs <- lapply(jsons,parse_encode_summary_json,row2print=rep1i[i])
  # tabs <- lapply(jsons,parse_encode_summary_json,row2print=ctl1i[i])
  contents <- sapply(tabs,function(x){x$content})
  dt <- as.data.table(contents)
  rownames(dt)<-tabs[[1]]$head
  colnames(dt)<-sample_base
  # 
  # dt <- as.data.table(tabs)
  # J <- dim(dt)[2]
  # arownames <- t(dt[,1])
  # dt <- dt[,seq(2,J,by=2),with=F]
  # # if (i==3){
  # #   dt <- dt[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
  # # } else {
  # #   dt <- dt[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
  # # }
  # 
  # rownames(dt) <- arownames
  # colnames(dt) <- sample_base
  out_fn <- sprintf('%s/%s_input_flagstat.tsv',args$out_dir,out_prefix[i])
  fwrite(dt,file=out_fn,quote=F,sep="\t",row.names=T,col.names=T)
}
