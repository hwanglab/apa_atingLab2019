source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
library(data.table)

samples <- c('HCT116','DKO')
bamd <- "./fastq/bt2"
fpaths <- lapply(samples,function(sample) {
  list.files(path=bamd,pattern = sprintf("%s_Mbd.bam$",sample),full.names = T)
})

mdup_fpaths <- lapply(fpaths,function(fpath) {
  mkdup_dir <- file.path(dirname(fpath),'mkdup')
  if (!dir.exists(mkdup_dir)) {
    dir.create(mkdup_dir)
  }
  return(picard_dedup(fpath,mkdup_dir))
})

bam_dt <- data.table(sample=samples,
                     fpath=mdup_fpaths)

flagstat_out <- lapply(bam_dt$fpath,samtools_flagstat)

bam_stats <- rbindlist(lapply(1:length(flagstat_out),function(i) {t(flagstat_out[[i]]$val)}))
colnames(bam_stats) <- unlist(flagstat_out[[1]]$item)

bam_stats_dt <- cbind(bam_dt,bam_stats)

outd<- file.path(bamd,"alignment_stats")

if (!dir.exists(outd)) {
	dir.create(outd)
}

tsv_file <- file.path(outd,'Supplementary_Table_S4_Mbd.tsv')

fwrite(bam_stats_dt,file=tsv_file,sep='\t')

message(tsv_file)
