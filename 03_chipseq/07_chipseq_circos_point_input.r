library(GenomicRanges)
library(data.table)
library(goldmine)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))

# -----------------------
message('preparing input setting ...')
wkd <- 'circos'
dir.create(wkd, showWarnings = FALSE)

encode_dir<-"fastq/kundaje_encode"
message(sprintf('encode directory [%s]',encode_dir))

# -----------------------
message("registering query genomic region ...")
if (F) {
  prepropa_rd <- '../01_polyAseq/01_wkd/out/03_CallApa/output/prepropa.rd'
  library(R.utils)
  
  red_ens <- loadToEnv(prepropa_rd)[["red_ens"]]
  goi_genes <- c('NFYA','DNAJB6','HEATR2')
  
  stopifnot(all(goi_genes %in% red_ens$tu$name))
  
  query_gr <- red_ens$tu[match(goi_genes,red_ens$tu$name),]
} else {
  #the following genomic ranges are provided by Vishal Nanvaty
  #https://bdcta.slack.com/archives/D8SAKEAP8/p1554390702003700
  query_gr <- GRanges(
    seqnames = c("chr7", "chr6","chr7"),
    ranges = IRanges(start=c(704741,41036213,157197743),
                     end=c(824968,41158907,157262434)),
    strand = c("+", "+","+"),
    gene = c("HEATR2","NFYA","DNAJB6"))
}
message('converting query region to a bed file format ...')
query_bed <- file.path(wkd,'query_3genes.bed')
grange_to_bed_file0(query_gr,query_bed)

# -----------------------
message('locating chipseq BAM files ...')
source(file.path(Sys.getenv('R_UTIL'),'lib_chipseq.r'))

bam_dt <- get_chipseq_bam(get_mount_dir())
bam_dt <- subset(bam_dt,sample %in% c('HC','DC'))

# -------------------
message('getting meand depth by sliding a window ...')

bam_dt$md_bed_fpath <- lapply(1:nrow(bam_dt),function(i) {
  bam <- bam_dt[i,fpath]
  extracted_bam <- sprintf('%s.query.bam',bam)
  
  message(sprintf('extracting bam [%s]',bam))
  samtools_extract_by_bed(bam,query_bed,extracted_bam,reuse=F,samtools_bin='samtools')
  
  out_prefix <- file.path(wkd,sprintf("%s_query",bam_dt[i,sample]))
  
  message(sprintf('computing a mean depth on bam [%s]',extracted_bam))
  
  md_bed_fpath <- mosdepth_md_by_window(extracted_bam,out_prefix=out_prefix,window=30,ncpu=1,minMapQuality=1)
  
  message(sprintf('Done[%s]',bam))
  
  if (file.exists(extracted_bam)) {
    #unlink(extracted_bam)
  }
  
  extracted_bai <- sprintf("%s.bai",extracted_bam)
  if (file.exists(extracted_bai)) {
    #unlink(extracted_bai)
  }
  
  return(md_bed_fpath)
})


lapply(1:nrow(bam_dt),function(i) {
  md_bed <- bam_dt[i,md_bed_fpath]
  sample <- bam_dt[i,sample]
  
  md_dt <- fread(md_bed[[1]],col.names=c('chr','start','end','value1'))
  md_gr <- makeGRanges(md_dt)
  
  ij <- findOverlaps(query_gr,md_gr)
  md_gr$query_gene <- NA
  md_gr$query_gene[subjectHits(ij)] <- query_gr$gene[queryHits(ij)]

  md_dt <- as.data.table(md_gr)[,c('seqnames','start','end','value1','query_gene')]

  colnames(md_dt) <- c('chr','start','end','value1','gene')
  md_dt <- md_dt[!is.na(gene),]
  md_dt$start <-  md_dt$start + 1
  dt_by_genes <- split(md_dt,by='gene')
  genes <- names(dt_by_genes)
  for (j in 1:length(dt_by_genes)) {
    csv_fn <- file.path(wkd,sprintf('%s_%s.csv',sample,genes[[j]]))
    
    fwrite(dt_by_genes[[j]][,c('chr','start','end','value1')],file=csv_fn,sep=',',col.names = FALSE)
    message(csv_fn)
  }
  return(NA)
})
