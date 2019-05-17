library(data.table)

source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))

bam_dt <- data.table(sample=c('HCT116','DKO'),
                     fpath=c('./star/HI.1159.005.Index_13.HCT116_Aligned.sortedByCoord.out.bam',
                             './star/HI.1159.005.Index_16.DKO_Aligned.sortedByCoord.out.bam'))

flagstat_out <- lapply(bam_dt$fpath,samtools_flagstat)

bam_stats <- rbindlist(lapply(1:length(flagstat_out),function(i) {t(flagstat_out[[i]]$val)}))
colnames(bam_stats) <- unlist(flagstat_out[[1]]$item)

bam_stats_dt <- cbind(bam_dt,bam_stats)
outd<-'./star/output'
if (!dir.exists(outd)) {
  dir.create(outd)
}

tsv_file <- './star/output/Supplementary_Table_S4_mRNA.tsv'

fwrite(bam_stats_dt,file=tsv_file,sep='\t')

message(tsv_file)
