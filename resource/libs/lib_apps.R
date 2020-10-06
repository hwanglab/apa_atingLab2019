#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018

library(data.table)
library(stringr)
library(DESeq2)
library(Rsubread)
library(Rsamtools)
library(reshape)
library(GenomicRanges)
library(jsonlite)
# library(gamlss)

add_chr_contig_header <- function(bam){
  out_bam <- sprintf('%s_chr.bam',bam)
  
  cmd <- sprintf("samtools view -H %s | sed -e 's/SN:\\([0-9XY]\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - %s > %s",bam,bam,out_bam)
  if (!file.exists(out_bam)){
    system(cmd)
  }

  cmd <- sprintf("samtools index %s",out_bam)
  bai <- sprintf("%s.bai",out_bam)
  if (!file.exists(bai)){
    system(cmd)
  }
  
  return(out_bam)
}

get_gunzip_fn <- function(file) {
	S = str_length(file)
	L = length(grep('\\.gz$',file,value=TRUE))
	if (L>0){
		file2 = substr(file,1,S-3)
		system(sprintf('gunzip -fc %s > %s',file,file2))
	} else {
		file2 = file
	}
	return(file2)
}

run_dreme_motif_search <- function(fasta_gz,control_fasta_gz,outdir) {
	message(sprintf('run dreme motif analysis on [%s;%s]',fasta_gz,control_fasta_gz))
	if (!dir.exists(outdir)){
		dir.create(outdir, showWarnings = FALSE)
	}

	fasta = get_gunzip_fn(fasta_gz)
	control_fasta = get_gunzip_fn(control_fasta_gz)

	cmd = sprintf("dreme -oc %s -p %s -n %s -png -dna -norc",outdir,fasta,control_fasta)
	message(cmd)
	system(cmd)
	message('Done.')
	if (fasta_gz!=fasta){unlink(fasta)}
	if (control_fasta_gz!=control_fasta){unlink(control_fasta)}
}

macs2_peak <- function(bam,sample_name,out_dir,peak_type='sharp',out_prefix='narrow',in_bam=NA,macs2_path=NA) {
  
  
  
  if (is.na(macs2_path)){
    macs2 <- "macs2"
  } else {
    macs2 <- macs2_path
  }
  
  cmd <- sprintf("%s callpeak -t %s -f BAM -g hs -n %s -B --outdir %s",macs2,bam,sample_name,out_dir)
  
  #cmd <- sprintf("%s callpeak -t %s -f BAM -g hs -n %s -B --keep-dup all -m 3 100 --outdir %s",macs2,bam,sample_name,out_dir)
  
  if (!is.na(in_bam)){cmd <- sprintf("%s -c %s",cmd,in_bam)}
  
  if (peak_type=='broad'){
    cmd <- sprintf("%s --broad",cmd)
  }
  
  message(cmd)
  system(cmd) #debug
  
  if (out_prefix=='narrow'){
    peak_tsv <- file.path(out_dir,sprintf('%s_peaks.narrowPeak',sample_name))
  } else {
    peak_tsv <- file.path(out_dir,sprintf('%s_peaks.xls',sample_name))
  }
  
  return(peak_tsv)
}

macs2_merge <- function(macs2_peaks) {
  
  dts <- list()
  if (startsWith(macs2_peaks[1], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  dts[[1]] <- fread(sprintf("%s -v ^# %s",grep_cmd,macs2_peaks[1]))
  
  if (startsWith(macs2_peaks[2], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  dts[[2]] <- fread(sprintf("%s -v ^# %s",grep_cmd,macs2_peaks[2]))
  
  peaks <- lapply(dts,function(x){
    GRanges(seqnames = x$chr,
            ranges = IRanges(start = x$start,
                             end = x$end,
                             names = x$name))
  })
  merged <- union(peaks[[1]],peaks[[2]])
  merged_dt <- data.table(GeneID=paste0('macs2_',1:length(merged)),
                          Chr=as.character(seqnames(merged)),
                          Start=start(merged),
                          End=end(merged),
                          Strand='*')
  return(merged_dt)
}

macs2_narrowpeak_concat <- function(macs2_peaks) {
  
  dts <- list()
  
  tmp_out1 <- sprintf('%s.tmp',macs2_peaks[1])
  
  if (endsWith(macs2_peaks[1], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n > %s",grep_cmd,macs2_peaks[1],tmp_out1)

  message(cmd)
  system(cmd)
  
  tmp_out2 <- sprintf('%s.tmp',macs2_peaks[2])
  if (endsWith(macs2_peaks[2], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n > %s",grep_cmd,macs2_peaks[2],tmp_out2)
  
  message(cmd)
  system(cmd)

  cmd <- sprintf('cat %s %s | sort -k1,1V -k2,2n -k3,3n',tmp_out1,tmp_out2)
  
  merged <- fread(cmd=cmd,col.names=c('Chr','Start','End'))
  merged$GeneID <- paste0('macs2_',1:dim(merged)[1])
  merged$Strand <- '*'

  unlink(tmp_out1)
  unlink(tmp_out2)
  
  return(merged)
}

macs2_narrowpeak_merge <- function(macs2_peaks) {
  
  dts <- list()
  if (endsWith(macs2_peaks[1], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n",grep_cmd,macs2_peaks[1])
  message(cmd)
  dts[[1]] <- fread(cmd=cmd)
  
  if (endsWith(macs2_peaks[2], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n",grep_cmd,macs2_peaks[2])
  message(cmd)
  dts[[2]] <- fread(cmd=cmd)
  
  peaks <- lapply(dts,function(x){
    GRanges(seqnames = x$chr,
            ranges = IRanges(start = x$start,
                             end = x$end))
  })
  merged <- union(peaks[[1]],peaks[[2]])
  merged_dt <- data.table(GeneID=paste0('macs2_',1:length(merged)),
                          Chr=as.character(seqnames(merged)),
                          Start=start(merged),
                          End=end(merged),
                          Strand='*')
  return(merged_dt)
}

get_total_reads_from_bam <- function(bam){
  cmd <- sprintf("samtools flagstat %s | head -n1 | cut -d\' \' -f1",bam)
  message(cmd)
  out <- system(cmd,intern = TRUE)
  num_reads <- as.integer(out)
  return(num_reads)
}

samtools_flagstat <- function(bam,samtools_bin='samtools',reuse=F){
  
  flag_stat_file <- paste0(bam,'.flagstat')
  cmd <- sprintf("samtools flagstat %s > %s",bam,flag_stat_file)
  message(cmd)
  
  if (!reuse | !file.exists(flag_stat_file)) {
    out <- system(cmd,intern = TRUE)
  }
  
  con  <- file(flag_stat_file, open = "r")
  
  stat.item <- list()
  stat.val <- list()
  
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    val <- as.numeric(str_extract(oneLine, pattern = "[0-9]+(?=\\s\\+)"))
    
    if (grepl('total', oneLine)) {
      stat.item <- c(stat.item,'total')
      stat.val <- c(stat.val,val)
    } else if (grepl('secondary',oneLine)) {
      stat.item <- c(stat.item,'secondary')
      stat.val <- c(stat.val,val)
    } else if (grepl('supplementary',oneLine)) {
      stat.item <- c(stat.item,'supplementary')
      stat.val <- c(stat.val,val)
    } else if (grepl('duplicates',oneLine)) {
      stat.item <- c(stat.item,'duplicates')
      stat.val <- c(stat.val,val)
    } else if (grepl('with itself and mate mapped',oneLine)) {
      stat.item <- c(stat.item,'both_mapped')
      stat.val <- c(stat.val,val)
    } else if (grepl('paired in sequencing',oneLine)) {
      stat.item <- c(stat.item,'paired_reads')
      stat.val <- c(stat.val,val)
    } else if (grepl('read1',oneLine)) {
      stat.item <- c(stat.item,'read1')
      stat.val <- c(stat.val,val)
    } else if (grepl('read2',oneLine)) {
      stat.item <- c(stat.item,'read2')
      stat.val <- c(stat.val,val)
    } else if (grepl('properly paired',oneLine)) {
      stat.item <- c(stat.item,'properly_paired')
      stat.val <- c(stat.val,val)
    } else if (grepl('singletons',oneLine)) {
      stat.item <- c(stat.item,'singletons')
      stat.val <- c(stat.val,val)
    } else if (grepl('with mate mapped to a different chr \\(mapQ>=5\\)',oneLine)) {
      stat.item <- c(stat.item,'mate_in_diffchr_mapq5')
      stat.val <- c(stat.val,val)
    } else if (grepl('with mate mapped to a different chr',oneLine)) {
      stat.item <- c(stat.item,'mate_in_diffchr')
      stat.val <- c(stat.val,val)
    }  else if (grepl('mapped',oneLine)) {
      stat.item <- c(stat.item,'mapped')
      stat.val <- c(stat.val,val)
    } else {
      message(sprintf('WARNING[%s]',oneLine))
    }
  } 
  
  close(con)
  
  stat_dt <- data.table(item=stat.item,val=stat.val)
  if (file.exists(flag_stat_file)) {
    if (!reuse) {
      unlink(flag_stat_file)
    }
  }
  return(stat_dt)
}


macs2_bdgdiff <- function(samp,
                          oprefix,
                          out_dir,
                          loglr_cutoff=1,
                          mode='TF',
                          macs2_bin='macs2'){
  # samp <- data.table(bam=c(bam1,bam2),
  # bdg=c(bdg1,bdg2),
  # bdgi=c(bdgi1,bdgi2),
  # nread=c(nread1,nread2),
  # nreadi=c(nreadi1,nreadi2),
  # group = c("ctrl","experiment"),
  # name=strsplit(args$name,",")[[1]])

  # ref: https://github.com/taoliu/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands
  
  if (mode == 'TF'){
    minlen <- 120
    maxgap <- 60
  } else if (mode =='histone'){
    minlen <- 120
    maxgap <- 60
  }
  
  cmd <- sprintf('%s bdgdiff --t1 %s --t2 %s --c1 %s --c2 %s --d1 %d --d2 %d -l %d -g %d --o-prefix %s --outdir %s',
                 macs2_bin,samp$bdg[1],samp$bdg[2],
                 samp$bdgi[1],samp$bdgi[2],
                 as.integer(round(samp$nread[1]/1e6)),
                 as.integer(round(samp$nread[2]/1e6)),
                 minlen,maxgap,oprefix,out_dir)
  message(cmd)
  system(cmd)
  
}

encode_narrowPeak <- function(encode_narrow_peak_file){
  #contig_name, start, end, peak_name, -10log10(qval), ?, fold-change, -log10(pval), -log10(qval), relative_summit_position
  if (file.exists(encode_narrow_peak_file)) {
    encode_npeak <- fread(encode_narrow_peak_file,header = F,col.names = c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))
    #browser()
    encode_npeak <- encode_npeak[order(chr,start,end),]
    
    return(encode_npeak)
  } else {
    stop(sprintf('the input file[%s] does not exist!',encode_narrow_peak_file))
  }
}

prep_deseq2 <- function(bam_dt,region2eval,countMultiMappingReads=F,allowMultiOverlap=F,fraction=F,ncore=8) {
  
  R <- dim(region2eval)[1]
  C <- dim(bam_dt)[1]
  
  message(sprintf('R=%d',R)) #debug
  message(sprintf('C=%d',C))
  
  fc <- mclapply(bam_dt$bam,
                 function(bam){
                   Rsubread::featureCounts(bam,
                                 annot.ext = region2eval,
				 countMultiMappingReads=countMultiMappingReads,
				 allowMultiOverlap=allowMultiOverlap,
				 fraction=fraction)
                   },
                 mc.cores=ncore)
  
  cts <- matrix(NA, nrow = R, ncol = C)
  
  rownames(cts) <- paste0(fc[[1]]$annotation$Chr,':',
                          fc[[1]]$annotation$Start,'-',
                          fc[[1]]$annotation$End)
  colnames(cts) <- bam_dt$name
  
  for (i in 1:C){
    cts[,i] <- fc[[i]]$counts
  }

  cts <- round(cts)
  
  coldata <- DataFrame(group=bam_dt$group)
  rownames(coldata) <- bam_dt$name
  coldata$group <- factor(coldata$group)
  return(list(cts=cts,coldata=coldata))
}

read_count <- function(bam,region2eval,
                       countMultiMappingReads=F,allowMultiOverlap=F,fraction=F,ignoreDup=F) {
  
  R <- dim(region2eval)[1]

  message(sprintf('R=%d',R)) #debug

  fc <- Rsubread::featureCounts(bam,
                          annot.ext = region2eval,
                          countMultiMappingReads=countMultiMappingReads,
                          allowMultiOverlap=allowMultiOverlap,
                          fraction=fraction,
                          ignoreDup=ignoreDup)
  
  cts <- matrix(NA, nrow = R, ncol = 1)
  
  rownames(cts) <- paste0(fc$annotation$Chr,':',
                          fc$annotation$Start,'-',
                          fc$annotation$End)
  
  colnames(cts) <- 'read_count'
  cts[,1] <- round(fc$counts)

  return(cts)
}

run_deseq2 <- function(deseq2_input,out_dir,shrink_lfc=FALSE,alpha=0.1) {
  
  dds <- DESeqDataSetFromMatrix(countData = deseq2_input$cts,
                                colData = deseq2_input$coldata,
                                design = ~group)
  dds <- DESeq(dds)

  if (shrink_lfc){ #apply fc shrinkage with apeglm
    res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
  } else {
    res <- results(dds,alpha=alpha)
  }
  
  # plotMA(res)

  res <- res[res$baseMean>0.,]
  df <- as.data.frame(res[order(res$pvalue),])

  read_cnt <- counts(dds)
  colnames(read_cnt) <- paste0('raw_',colnames(read_cnt))
  raw_read_cnt_cols <- colnames(read_cnt)
  df <- merge(df,read_cnt,by=0,all.x=TRUE)

  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  read_cnt <- counts(dds,normalized=TRUE)
  colnames(read_cnt) <- paste0('norm_',colnames(read_cnt))
  df <- merge(df,read_cnt,by=0,all.x=TRUE)
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  
  names(df)[1] <- sprintf("region\t%s",names(df)[1])
  out_file <- file.path(out_dir,'deseq2.tsv')
  df0 <- df 
  write.table(df,
              quote=FALSE,sep="\t",
              na="NA",
              row.names=TRUE,
              col.names=TRUE,
              file=out_file)
  
  return(df0)
}

picard_dedup <- function(bam,bam_out_dir,min_mapq=0,picard_jar_path='/media/sammy/projects/apps/picard/latest/picard.jar', reuse=TRUE){
  
  
  bam_fname <- basename(bam)
  mdup_bam <- file.path(bam_out_dir,sprintf('%s.mdup.bam',bam_fname))
  out_metrics <- file.path(bam_out_dir,sprintf('%s.mdup.metrics',bam_fname))
  
  if (is.na(picard_jar_path)){
    cmd_prefix <- "java -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT"
  } else {
    cmd_prefix <- sprintf("java -jar %s MarkDuplicates VALIDATION_STRINGENCY=SILENT",picard_jar_path)
  }
  if(file.exists(mdup_bam)&reuse) {
    message('reuse the previous result ...')
  } else {
    cmd <- sprintf('%s I=%s O=%s M=%s',cmd_prefix,bam,mdup_bam,out_metrics)
    message(cmd)
    system(cmd)
  }
  
  if (min_mapq>0) {
    mdup_mapq_bam <- file.path(bam_out_dir,sprintf('%s.mapq%d.bam',bam_fname,min_mapq))
    ncpu <- 4
    #cmd <- sprintf('samtools view -@ %d -b -h -F0x404 %s > %s',ncpu,mdup_bam,dedup_bam)
    if(file.exists(mdup_mapq_bam)&reuse){
      message('reuse the previous result ...')
    } else {
      cmd <- sprintf('samtools view -@ %d -b -h -q %d %s > %s',ncpu,min_mapq,mdup_bam,mdup_mapq_bam)
      message(cmd)
      system(cmd)
    }
    out_bam <- mdup_mapq_bam
  } else {
    out_bam <- mdup_bam
  }
  
  cmd <- sprintf('samtools index %s',out_bam)
  bai <- sprintf('%s.bai',out_bam)
  if (file.exists(bai)&reuse) {
    message('reuse the previous result ...')
  } else {
    system(cmd)
  }
  return(out_bam)
}

log2fc <- function(ctrl_cnt,expl_cnt){
  if (ctrl_cnt>0){
    if (expl_cnt>0){
      fcL2 <- log2(expl_cnt/ctrl_cnt)
    } else {
      fcL2 <- -Inf
    }
  } else {
    fcL2 <- Inf
  }
  return(fcL2)
}

bam_to_bdg <- function(bam,out_bdg,bin_size=5,bamCoverage_bin='bamCoverage'){
  
  cmd <- sprintf('%s --bam %s --binSize %d --outFileFormat bedgraph --outFileName %s',bamCoverage_bin,bam,bin_size,out_bdg)
  
  message(cmd)
  
  if (!file.exists(out_bdg)) {
    system(cmd)
  }
}

bedtools_bamtobed_qc <- function(in_bam,tag_bed,bedtools_bin="bedtools"){
  
  cmd <- sprintf("samtools view -hb -f 0x404 %s | %s bamtobed -i stdin > %s",
                 in_bam,
                 bedtools_bin,
                 tag_bed)
  message(cmd)
  system(cmd)
}

bedtools_bamtobed <- function(in_bam,out_bed,bedtools_bin="bedtools",optStr=NA){
  cmd <- sprinf("%s bamtobed",bedtools_bin)
  if (is.na(optStr)){
    cmd <- sprintf("%s %s",cmd,optStr)
  }
  cmd <- sprintf("%s -i %s > %s",cmd,in_bam,out_bed)
  message(cmd)
  system(cmd)
}


parse_encode_summary_json <- function(json_fn,row2print){
  
  js1 <- fromJSON(json_fn)

  head <- strsplit(js1$qc_files$head[row2print],'\t')[[1]]

  content <- strsplit(js1$qc_files$contents[row2print],'\t')[[1]]
  
  
  return (data.table(head=head,content=content))
}

diffReps <- function(bed_co,bed_tr,outfile,
                     genome="hg19",
                     test_model='nb',
                     window=200,
                     frag=150,
                     diffReps_bin='diffReps.pl'){
  cmd <- diffReps_bin
  
  cmd <- sprintf("%s -tr %s -co %s -gn hg19 -re %s -me %s --frag %d --nproc 8 --window %d",diffReps_bin,bed_tr,bed_co[1],outfile,test_model,frag,window)
  
  message(cmd)
  system(message)
}


macs2_read_narrow_peak <- function(narrow_peak_fn){
  macs_narrow_peak_cols <- c('chr','start','end','name','ucsc_score','strand','fc','pval_nlog10','qval_nlog10','summit')
  return(fread(narrow_peak_fn,col.names=macs_narrow_peak_cols))
}

sort_matrix_by_rowsum <- function(m2) {
  L <- dim(m2)[2]
  m2.row.sum <- cbind(m2,rowSums(m2))
  o2 <- rev(order(m2.row.sum[,L+1]))
  m2.row.sum <- m2.row.sum[o2,]
  return(m2.row.sum[,1:L])
}

getMeanDepthFromBam <- function(bam,get,minMapQuality=30,debug=F)
{
  if(debug){browser()}
  message(bam)
  param <- ApplyPileupsParam(which=get, 
                             what=c("seq", "qual"),
                             yieldBy="position",
                             yieldAll=TRUE,
                             minMapQuality=minMapQuality,
                             maxDepth=1000)
  
  fls <- PileupFiles(bam, param=param)
  calcInfo <- function(x)
  {
    if(debug){browser()}
    info <- apply(x[["seq"]], 2, function(y) {
      y <- y[c("A", "C", "G", "T"), , drop=FALSE]
      cvg <- colSums(y)
    })
    info
  }
  res <- applyPileups(fls, calcInfo, param=param)
  genecov <- t(do.call(cbind,res))
  return(list(mean_dp=mean(genecov),cov_vec=genecov))
}

# --------------------------------------------------------------------
grange_to_bed_file <- function(gr, file, name=NULL, offset=0)
{
  if(file.exists(file)){file.remove(file)}
  
  df <- data.frame(chr=as.character(seqnames(gr)),
                   start=start(gr)-1+offset,
                   end=end(gr)-offset,
                   anot=sprintf("%d_%s_%s",gr$test_id,gr$type,gr$gene_name))
  if(!is.null(name))
  {
    df$name <- values(gr)[,name]
  }
  #df <- df[df$chr %in% c(sapply(seq(1,22),function(x) paste("chr",x,sep="")),"chrX","chrY"),]
  write.table(df, file=file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
}

parse_apa_bed_file <- function(mosdepth_region_fpath) {
  
  dt <- fread(cmd=sprintf("zcat %s",mosdepth_region_fpath),
              col.names=c('chr','start','end','annot','mean_depth'))
  
  tmp <- tstrsplit(dt$annot,'_',fixed=TRUE)
  dt$test_id <- as.numeric(tmp[[1]])
  dt$type <- tmp[[2]]
  dt$gene_name <- tmp[[3]]
  dt$annot <- NULL
  return(dt)
}

compPauRatioFromBam_mosdepth <- function(bam,apa_test_gr,minMapQuality=30, min_depth=0.,ncpu=4,mosdepth_bin="mosdepth", comp_ratio=T,debug=F)
{
  #apa_test_reg can be multiple regions
  # test_id chrom    start      end   type gene_name
  # 1:       1  chr4  8160406  8160455 shared    ABLIM2
  # 2:       1  chr4  8147865  8148132 unique    ABLIM2
  # 3:       2 chr16 89371613 89371752 shared   ANKRD11
  # 4:       2 chr16 89365283 89367370 unique   ANKRD11
  # 5:       3  chr2  9347148  9347349 shared     ASAP2
  # 6:       3  chr2  9347350  9349629 unique     ASAP2
  # 7:       4 chr12 58022496 58022686 shared  B4GALNT1
  # 8:       4 chr12 58022210 58022495 unique  B4GALNT1
  # 9:       5 chr12 58022496 58022686 shared  B4GALNT1
  # 10:       5 chr12 58019917 58020744 unique  B4GALNT1
  # 11:       6  chr8 41525785 41526082 shared      ANK1
  # 12:       6  chr8 41529872 41530430 shared      ANK1
  # 13:       6  chr8 41510744 41513272 unique      ANK1
  
  if(debug){browser()}
  message(bam)
  
  apa_test_bed_fpath <- paste0(bam,'.goi.bed')
  
  out_prefix <- paste0(bam,'.md')
  region_bed_fpath <- sprintf("%s.regions.bed.gz",out_prefix)
  if (!file.exists(region_bed_fpath)) {

    if (file.exists(apa_test_bed_fpath)) {
      unlink(apa_test_bed_fpath)
    }
    
    grange_to_bed_file(apa_test_gr,apa_test_bed_fpath)
  
    if (!file.exists(apa_test_bed_fpath)){
      stop('%s is not available and probably goi_gr to bed conversion is failed',apa_test_bed_fpath)
    }
    
    cmd <- sprintf('%s -t %d -b %s -n -Q %d %s %s',mosdepth_bin,ncpu,apa_test_bed_fpath,minMapQuality,out_prefix,bam)
    message(cmd)
    
    system(cmd)
  }
  
  bnd <- parse_apa_bed_file(region_bed_fpath)
  bnds <- split(bnd,by=c('test_id'),sorted=FALSE,drop=FALSE)
  
  eps2 <- 1e-7
  rm(bnd)
  
  pau_ratios <- sapply(1:length(bnds), function(i){
    reg <- bnds[[i]]
    
    width <- reg[type=='shared',end] - reg[type=='shared',start]+ 1
    shared_mdp0 <- sum(reg[type=='shared',mean_depth] * width) / sum(width)
    
    if (!is.finite(shared_mdp0)) {
      shared_mdp0 <- 0
    }
    shared_mdp <- shared_mdp0 + eps2 #shared
    
    width <- reg[type=='unique',end] - reg[type=='unique',start]+ 1
    unique_mdp <- sum(reg[type=='unique',mean_depth] * width) / sum(width)
    
    if (!is.finite(unique_mdp)) {
      unique_mdp <- 0
    }
    
    if (comp_ratio) {
      if (shared_mdp > 0) {
        paur <- c(unique_mdp/shared_mdp,shared_mdp,unique_mdp)
      } else {
        paur <- c(NA,shared_mdp,unique_mdp)
      }
    } else {
      paur <- c(-1.,shared_mdp,unique_mdp)
    }
    
    return(paur)
  })
  colnames(pau_ratios)<-names(bnds)
  rownames(pau_ratios) <-c('ratio','common','unique')
  message(sprintf("Done mean depth calculation on [%s]",bam))

  unlink(apa_test_bed_fpath)

  return(pau_ratios)
}

compPauRatioFromBam <- function(bam,apa_test_gr,minMapQuality=30, min_depth=0.,comp_ratio=T,debug=F)
{
  #apa_test_reg can be multiple regions
  # test_id chrom    start      end   type gene_name
  # 1:       1  chr4  8160406  8160455 shared    ABLIM2
  # 2:       1  chr4  8147865  8148132 unique    ABLIM2
  # 3:       2 chr16 89371613 89371752 shared   ANKRD11
  # 4:       2 chr16 89365283 89367370 unique   ANKRD11
  # 5:       3  chr2  9347148  9347349 shared     ASAP2
  # 6:       3  chr2  9347350  9349629 unique     ASAP2
  # 7:       4 chr12 58022496 58022686 shared  B4GALNT1
  # 8:       4 chr12 58022210 58022495 unique  B4GALNT1
  # 9:       5 chr12 58022496 58022686 shared  B4GALNT1
  # 10:       5 chr12 58019917 58020744 unique  B4GALNT1
  # 11:       6  chr8 41525785 41526082 shared      ANK1
  # 12:       6  chr8 41529872 41530430 shared      ANK1
  # 13:       6  chr8 41510744 41513272 unique      ANK1
  
  if(debug){browser()}
  message(bam)
  param <- ApplyPileupsParam(which=apa_test_gr,
                             what=c("seq", "qual"),
                             yieldBy="position",
                             yieldAll=TRUE,
                             minMapQuality=minMapQuality,
                             maxDepth=5000)
  
  fls <- PileupFiles(bam, param=param)
  
  calcInfo <- function(x) {
    info <- apply(x[["seq"]], 2, function(y) {
      y <- y[c("A", "C", "G", "T"), , drop=FALSE] #TODO: what about indels?
      cvg <- colSums(y)
    })
    info
  }
  res <- applyPileups(fls, calcInfo, param=param)
  
  if(debug){browser()}

  cvg_end <- cumsum(width(apa_test_gr))
  L <- length(apa_test_gr)

  apa_test_gr<-as.data.table(apa_test_gr)
  apa_test_gr$cvg_start <- c(1,cvg_end[1:(L-1)]+1)
  apa_test_gr$cvg_end <- cvg_end
  
  bnd <- apa_test_gr[,.(start=min(cvg_start),end=max(cvg_end)),by=c('test_id','type')]
  
  bnds <- split(bnd,by=c('test_id'),sorted=FALSE,drop=FALSE)
  
  genecov <- t(do.call(cbind,res))
  
  eps2 <- 1e-7
  rm(bnd)
  rm(res)
  
  pau_ratios <- sapply(1:length(bnds), function(i){
    reg <- bnds[[i]]
    shared_mdp0 <- mean(genecov[reg[type=='shared',start]:reg[type=='shared',end]])
    if (!is.finite(shared_mdp0)) {
      shared_mdp0 <- 0
    }
    shared_mdp <- shared_mdp0 + eps2 #shared
    
    unique_mdp <- mean(genecov[reg[type=='unique',start]:reg[type=='unique',end]]) #unique
    if (!is.finite(unique_mdp)) {
      unique_mdp <- 0
    }
    
    if (comp_ratio) {
      if (unique_mdp > shared_mdp) {
        if (unique_mdp >= min_depth){
          paur <- c(1.,shared_mdp,unique_mdp)
        } else {
          paur <- c(NA,shared_mdp,unique_mdp)
        }
      } else {
        if (shared_mdp >= min_depth){
          paur <- c(unique_mdp/shared_mdp,shared_mdp,unique_mdp)
        } else {
          paur <- c(NA,shared_mdp,unique_mdp)
        }
      }
    } else {
      paur <- c(-1.,shared_mdp,unique_mdp)
    }
    
    return(paur)
  })
  colnames(pau_ratios)<-names(bnds)
  rownames(pau_ratios) <-c('ratio','common','unique')
  message(sprintf("Done meand depth calculation on [%s]",bam))
  return(pau_ratios)
}

gamlss_general_with_tryCatch <- function(dat2) {
  fitspm <- tryCatch(
    {
      fitspm <- gamlss(x ~ y, data=dat2, family=NO)
    },
    error=function(cond){
      message('gamlss failed in handling input dat2!')
      return(NULL)
    },
    finally={
      message('gamlss processes input dat2!')
    }
  )
  return(fitspm)
}

gamlss_with_tryCatch <- function(dat2,mode=0) {
  fitspm <- tryCatch(
    {
      if (mode==1){
        fitspm <- gamlss(pau ~ 1, data=dat2)
      } else {
        fitspm <- gamlss(nbeta ~ pau, data=dat2, family=NO)
      }
    },
    error=function(cond){
      message(sprintf('gamlss failed in handling input dat2![%s]',cond))
      return(NA)
    },
    finally={
      message('gamlss processes input dat2!')
    }
  )
  return(fitspm)
}

regression_fit <- function(dat2,method='spline_monotonic',debug=False) {
  
  if (debug){browser()}
  
  if (method=='spline_monotonic'){
    fitspm_up <- gamlss(pau ~ pbm(nbeta,mono='up'), data=dat2)
    fitspm_dn <- gamlss(pau ~ pbm(nbeta,mono='down'), data=dat2)
    
    if (fitspm_up$sbc < fitspm_dn$sbc){
      fitspm <- fitspm_up
    } else {
      fitspm <- fitspm_dn
    }
  } else { #use a simple linear regression
    fitspm <- gamlss_with_tryCatch(dat2)
  }
  
  if (data.class(fitspm)!="gamlss") {
    regfit_info <- list(rsq=0.,pval=1.)
  } else {
    if (F){
      plot(pau ~ nbeta, data=dat2, col='lightblue')
      lines(fitted(fitspm)[order(dat2$nbeta)] ~ dat2$nbeta[order(dat2$nbeta)])
      summary(fitspm)
    }
  
    fitspmi <- list()
    fitspmi$R2 <- Rsq(fitspm)
    fitspmi$rss <- sum(residuals(fitspm, what='mu',type='simple')^2)
    
    fit0 <- gamlss_with_tryCatch(dat2,mode=1)
    if (data.class(fit0)!="gamlss") {
      regfit_info <- list(rsq=0.,pval=1.)
    } else {
      if (F){
        plot(pau ~ nbeta, data=dat2, col='lightblue')
        lines(fitted(fit0)~dat2$nbeta)
        summary(fit0)
      }
  
      fit0i <- list()
      fit0i$R2 <- Rsq(fit0)
      fit0i$rss <- sum(residuals(fit0, what='mu',type='simple')^2)
      
      n <- dim(dat2)[1]
      
      if (method=='spline_monotonic') {
        d <- round(edf(fitspm)[[1]][1])
        if (d==1) {d <- 2}
      } else {
        d <- 2
      }
      
      Fvalue <- (fit0i$rss-fitspmi$rss)/(d-1)/fitspmi$rss*(n-d)
      regfit_info <- list(rsq=fitspmi$R2,pval=pf(Fvalue,2,7,lower.tail=FALSE)[1])
    }
  }
  return(regfit_info)
}

regression_fit_general <- function(dat2,method='spline_monotonic',debug=False) {
  
  if (debug){browser()}
  
  if (method=='spline_monotonic'){
    fitspm_up <- gamlss(y ~ pbm(x,mono='up'), data=dat2)
    fitspm_dn <- gamlss(y ~ pbm(x,mono='down'), data=dat2)
    
    if (fitspm_up$sbc < fitspm_dn$sbc){
      fitspm <- fitspm_up
    } else {
      fitspm <- fitspm_dn
    }
  } else { #use a simple linear regression
    fitspm <- gamlss_general_with_tryCatch(dat2)
  }
  
  if (is.null(fitspm)) {
    regfit_info <- list(rsq=0.,pval=1.)
  } else {
    if (F){
      plot(y ~ x, data=dat2, col='lightblue')
      lines(fitted(fitspm)[order(dat2$x)] ~ dat2$x[order(dat2$x)])
      summary(fitspm)
    }
    
    fitspmi <- list()
    fitspmi$R2 <- Rsq(fitspm)
    fitspmi$rss <- sum(residuals(fitspm, what='mu',type='simple')^2)
    
    fit0 <- gamlss(y ~ 1, data=dat2)
    if (F){
      plot(y ~ x, data=dat2, col='lightblue')
      lines(fitted(fit0)~dat2$x)
      summary(fit0)
    }
    
    fit0i <- list()
    fit0i$R2 <- Rsq(fit0)
    fit0i$rss <- sum(residuals(fit0, what='mu',type='simple')^2)
    
    n <- dim(dat2)[1]
    
    if (method=='spline_monotonic') {
      d <- round(edf(fitspm)[[1]][1])
      if (d==1) {d <- 2}
    } else {
      d <- 2
    }
    
    Fvalue <- (fit0i$rss-fitspmi$rss)/(d-1)/fitspmi$rss*(n-d)
    regfit_info <- list(rsq=fitspmi$R2,pval=pf(Fvalue,2,7,lower.tail=FALSE)[1])
  }
  return(regfit_info)
}

pathoscope2 <- function(ps2_info,reads,outdir,readL=0,max_num_ti_report=100,debug=F) {
  if (debug) {browser()}
  
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  
  if (!file.exists(ps2_info$conf_path)) {
    stop(sprintf('check PS2 conf file [%s] exists!',ps2_info$conf_fpath))
  }
  
  cmd <- sprintf("python %s -c %s -r Illumina -a bt2 -t %d -o %s -p %d --adjreflen -b very-sensitive-local",ps2_info$prog_path,ps2_info$conf_path,max_num_ti_report,outdir,ps2_info$ncpu)
  
  if (readL > 0) {
    cmd <- sprintf("%s -L %d",cmd,readL)
  }
  
  if (length(reads)>1) {
    if (file.exists(reads[[1]]) & file.exists(reads[[2]])) {
      cmd <- sprintf('%s -1 %s -2 %s',cmd,reads[[1]],reads[[2]])
    } else {
      if (!debug) {
	      stop('check if both reads [%s] and [%s] exist!',reads[[1]],reads[[2]])
      }
    }
  } else {
    if (file.exists(reads[[1]])) {
      cmd <- sprintf('%s -u %s',cmd,reads[[1]])
    } else {
      if (!debug) {
      	stop('check if read [%s] exists!',reads[[1]])
      }
    }
  }
  message(cmd)
  if (!debug) {
  	system(cmd)
  }
}

pathoscope2_get_rankings <- function(ps2_report_file,sample=NA) {
  # browser()
  ps2_dir <- Sys.getenv("PATHOSCOPE2")
  pybin2litereport <- file.path(ps2_dir,'pathoreport','postprocess_ps2_report_lite.py')
  ps2_lite_tmp_fn <- paste0(ps2_report_file,'.tmp')
  cmd <- sprintf('python %s -i %s -o %s',
                 pybin2litereport, 
                 ps2_report_file, 
                 ps2_lite_tmp_fn)
  
  if (!is.na(sample)) {
    cmd <- sprintf('%s -s %s',sample)
  }
  
  message(cmd)
  system(cmd)
  ps2_ranking_dt <- fread(ps2_lite_tmp_fn)
  unlink(ps2_lite_tmp_fn)
  return(ps2_ranking_dt)
}

pathoqc <- function(ps2_info,raw_reads,read_base,outdir,ncpu=4,debug=F) {

	qc_read1 <- file.path(outdir,paste0(read_base,'_1_qc.fq'))
	
	if (!file.exists(qc_read1)){
		
		if (length(raw_reads)>1) {
			qc_reads <- file.path(outdir,paste0(read_base,c('_1_qc.fq','_2_qc.fq')))
			cmd <- sprintf("python %s -1 %s -2 %s -t 33 -m 25 -e 50 -g 1 -a Y -a2 Y -d 0 -q 1 -p %d -o %s",
				ps2_info$pathoqc_bin,
				raw_reads[[1]],
				raw_reads[[2]],
				ncpu,
				outdir)
			message(cmd)
			if (!debug){
			  system(cmd)
			}
		} else {
				qc_reads <- file.path(outdir,paste0(read_base,'_1_qc.fq'))
				cmd <- sprintf("python %s -1 %s -t 33 -m 25 -e 50 -g 1 -a Y -d 0 -q 1 -p %d -o %s",
				ps2_info$pathoqc_bin,
				raw_reads[[1]],
				ncpu,
				outdir)
				
				message(cmd)
				
				if (!debug) {
				  system(cmd)
				}
				
		}
	} else {
		if (length(raw_reads)>1) {
			qc_reads <- file.path(outdir,paste0(read_base,c('_1_qc.fq','_2_qc.fq')))
		} else {
			qc_reads <- file.path(outdir,paste0(read_base,'_1_qc.fq'))
		}
	}
	return(qc_reads)
}

samtools_extract_by_bed <- function(bam,bed,out_bam,reuse=F,samtools_bin='samtools'){
  if (reuse & file.exists(out_bam)) {
    message(sprintf('reuse prev out_bam [%s] ...',out_bam))
  } else {
    cmd <- sprintf("%s view -O BAM -q 1 -F 4 -@ 4 -L %s -o %s %s",samtools_bin,bed,out_bam,bam)
    message(cmd)
    system(cmd)
  }
}

bedtools_genomecov <- function(bam,out_file,reuse=F,bedtools_bin='bedtools') {
  if (reuse & file.exists(out_file)) {
    message(sprintf('reuse prev out_file [%s] ...',out_file))
  } else {
    cmd <- sprintf("%s genomecov -bga -split -ibam %s | gzip -fc > %s",bedtools_bin,bam,out_file)
    message(cmd)
    system(cmd)
  }
}

mosdepth_md_in_bed <- function(bam,roi_bed,mosdepth_bin,ncpu=1,minMapQuality=1) {
  
  out_prefix <- paste0(bam,'.md')
  
  if (!file.exists(roi_bed)){
    stop('%s is not available and probably goi_gr to bed conversion is failed',roi_bed)
  }
  
  cmd <- sprintf('%s -t %d -b %s -n -Q %d %s %s',mosdepth_bin,ncpu,roi_bed,minMapQuality,out_prefix,bam)
  message(cmd)
  system(cmd)
  
  region_bed_fpath <- sprintf("%s.regions.bed.gz",out_prefix)
  if (!file.exists(region_bed_fpath)) {
    stop('%s is not available and probably mosdepth run on [%s] failed',region_bed_fpath,bam)
  }
  
  bnd <- parse_apa_bed_file(region_bed_fpath)
  
  if (file.exists(roi_bed)) {
    unlink(roi_bed)
  }
  
}

mosdepth_md_by_window <- function(bam,out_prefix=NA,mosdepth_bin="mosdepth",window=30,ncpu=1,minMapQuality=1) {
	
	if (is.na(out_prefix)) {
		out_prefix <- paste0(bam,'.md')
	}

  cmd <- sprintf('%s -t %d -b %d -n -Q %d %s %s',mosdepth_bin,ncpu,window,minMapQuality,out_prefix,bam)
  message(cmd)
  system(cmd)
  
  region_bed_fpath <- sprintf("%s.regions.bed.gz",out_prefix)
  
  if (!file.exists(region_bed_fpath)) {
    stop('%s is not available and probably mosdepth run on [%s] failed',region_bed_fpath,bam)
  }
	return(region_bed_fpath)
}
