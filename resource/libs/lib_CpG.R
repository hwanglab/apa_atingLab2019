#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018

library(data.table)
library(stringr)
library(GenomicRanges)

lcpg_load_ucsc_CpG_bed <- function(cpg_fn,decoy=FALSE) {
	#cpg_fn = '../../../apa/data/cpgs/cpgIslandExt_hg37_20171027.txt.gz'
	message(sprintf('loading [%s] into data.table ...',cpg_fn))
  if (decoy) {
	  cpgDt = fread(input=sprintf('zcat %s',cpg_fn),header = TRUE)
  } else {
    cpgDt = fread(input=sprintf("zcat %s | awk '$2 !~ /_/'",cpg_fn),header = TRUE)
  }
	cpgDt$"#bin"=NULL
	cpgDt$idx = 1:dim(cpgDt)[1]
	cpgDt$chromStart = as.integer(cpgDt$chromStart)
	cpgDt$chromEnd = as.integer(cpgDt$chromEnd)
	message('Done.')
	cpgDt
}

add_flanking <- function(cpgDt,flanking_bp=10){
  
  cpgDt$chromStart = cpgDt$chromStart-flanking_bp
  cpgDt$chromEnd = cpgDt$chromStart+flanking_bp
  cpgDt
}


lcpg_2grange <- function(cpgDt,flanking_bp=0){
  
  GRanges(seqnames=Rle(cpgDt$chrom),
          ranges=IRanges(cpgDt$chromStart-flanking_bp,cpgDt$chromEnd+flanking_bp))
}

load_cpgislandExt <- function(cpgisland_fn,decoy=FALSE) {
  message(sprintf('Reading %s',cpgisland_fn))
  stopifnot(file.exists(cpgisland_fn))
  if (decoy) {
    cpgDt = fread(input=sprintf('zcat %s',cpgisland_fn),header = TRUE)
  } else {
    cpgDt = fread(input=sprintf("zcat %s | awk '$2 !~ /_/'",cpgisland_fn),header = TRUE)
  }
  cpgDt$`#bin` <- NULL
  return(cpgDt)
}

get_cpg_movavg <- function(groi,edge_offset=0,wlen=50){ #granges of interest
# -----------------
  # <testing codes>
  # cg_gr <- GRanges(seqnames = 'chr1',IRanges(start=c(1,3,5,9),end=c(2,4,6,10)),strand='*',score=c(1,1,1,1))
  # cg_cvg <- coverage(cg_gr,weight="score")
  # groi_pos <- GRanges(seqnames='chr1',IRanges(start=1,end=15),strand='+')
  # binned_windows <- slidingWindows(groi_pos, width=4, step=2)
  # cg_movavg <- binnedAverage(binned_windows[[1]], cg_cvg, "binned_avg")
  # GRanges object with 7 ranges and 1 metadata column:
  #   seqnames    ranges strand | binned_avg
  # <Rle> <IRanges>  <Rle> |  <numeric>
  #   [1]     chr1  [ 1,  4]      + |          1
  # [2]     chr1  [ 3,  6]      + |          1
  # [3]     chr1  [ 5,  8]      + |        0.5
  # [4]     chr1  [ 7, 10]      + |        0.5
  # [5]     chr1  [ 9, 12]      + |        0.5
  # [6]     chr1  [11, 14]      + |          0
  # [7]     chr1  [13, 15]      + |          0
# -----------------
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(goldmine)
  
  message('1. collapsing overlapped genomic regions ...')
  if (edge_offset>0){
    start(groi) <- start(groi)-edge_offset
    end(groi) <- end(groi)+edge_offset
  }
  groi <- reduce(groi)
  
  message('2. loading ref genome sequence ...')
  groi$name <- sprintf('region:%d',1:length(groi))
  groi_pos <- groi
  strand(groi_pos) <- '+' #to extract only watson strand
  seqs <- getSeq(Hsapiens,groi_pos)
  names(seqs) <- groi$name
  
  message('3. search every CG location in the sequences derived from expanded/collapsed PRs')
  matched_grs <- vmatchPattern('CG',seqs) #list of data.frame
  # browser() #debug
  message('4. getting the actual genomic coordinates ...')
  cgloc_dt <- lapply(1:length(matched_grs),function(i){
    pri <- groi_pos[i]
    gri <- matched_grs[[i]]
    if (isEmpty(gri)){
      return(NA)
    } else {
      gr <- GRanges(seqnames = seqnames(pri),
                    strand = '*',
                    ranges = IRanges(start=start(pri)+start(gri)-1,
                                     end=start(pri)+end(gri)-1)
      )
      return(as.data.table(gr))
    }
    }
  )
  
  cgloc_dt <- cgloc_dt[!is.na(cgloc_dt)]
  
  message('5. combining CpG positions across all PRs and generates a read-like pileup ...')
  
  cgloc_dt <- rbindlist(cgloc_dt)
  cgloc_dt$width <- NULL
  names(cgloc_dt) <- c('chr','start','end','strand')
  cg_gr <- makeGRanges(cgloc_dt,strand=F)
  cg_gr$score <- 1. #set a default magnitude
  
  message('6. generating sliding binned windows ...')
  binned_windows <- slidingWindows(groi_pos, width=wlen, step=wlen/2)
  bwdts <- lapply(binned_windows,function(binned_window){
    as.data.table(binned_window)
  })
  
  message('7. merging binned wondiws ...')
  bwdt <- rbindlist(bwdts)
  bwdt$width <- NULL
  names(bwdt) <- c('chr','start','end','strand')
  prgr <- makeGRanges(bwdt,strand=F)
  
  message('8. performing a moving average on CpG sites on unstranded prs templates ...')
  
  cg_cvg <- coverage(cg_gr,weight="score")
  cg_movavg <- binnedAverage(prgr, cg_cvg, "binned_avg")
  return(cg_movavg)
}

get_ucsc_cg_island_avg <- function(ucsc_cgi_fn){
  
  cgdt <- load_cpgislandExt(ucsc_cgi_fn,decoy=T)
  cg_movavg <- GRanges(seqnames = cgdt$chrom,
                       IRanges(start=cgdt$chromStart,
                               end=cgdt$chromEnd),
                       strand="*",
                       binned_avg=as.double(cgdt$perCpg))
return(cg_movavg)
}