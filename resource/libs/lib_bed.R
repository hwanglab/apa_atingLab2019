#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018

library(data.table)
library(stringr)
library(GenomicRanges)

overlap_genomic_region <- function(qchrom,qstart,qend,RegionT,hit='best') {
	#case 1.
	j1 = which(RegionT$chrom == qchrom &
						 	RegionT$chromStart < qstart &
						 	qstart < RegionT$chromEnd)
	overlap1 = RegionT[j1,chromEnd] - qstart
	
	#case 2.
	j2 = which(RegionT$chrom == qchrom &
						 	qstart < RegionT$chromStart&
						 	RegionT$chromEnd <= qend)
	overlap2 = RegionT[j2,chromEnd] - RegionT[j2,chromStart]
	
	#case 3.
	j3 = which(RegionT$chrom == qchrom &
						 	RegionT$chromStart < qend &
						 	qend < RegionT$chromEnd)
	overlap3 = qend - RegionT[j3,chromStart]
	
	#case 4.
	j4 = which(RegionT$chrom == qchrom &
						 	RegionT$chromStart < qstart &
						 	qend <= RegionT$chromEnd)
	overlap4 = rep(qend - qstart, length(j4))
	
	if ('index' %in% colnames(RegionT)) {
		js = c(RegionT[j1]$index,RegionT[j2]$index,RegionT[j3]$index,RegionT[j4]$index)
	} else {
		js = c(j1,j2,j3,j4)
	}
	
	overlaps = c(overlap1,overlap2,overlap3,overlap4)
	if (any(js)){
		ui = which(!duplicated(js))
		dt = data.table(idx=js[ui],overlap=overlaps[ui])
		if (hit == 'best'){
			setorder(dt, -overlap)
			return(dt[1]$idx)
		}
		else if (hit == 'all'){
			return(paste(dt$idx,collapse=","))
		}
	} else {
		return(NA)
	}
}

overlap_genomic_region2 <- function(qchrom,qstart,qend,qmid,qmid_orient,RegionT,hit='all') {
	#case 1.
	j1 <- which(RegionT$chrom == qchrom &
						 	RegionT$chromStart < qstart &
						 	qstart < RegionT$chromEnd)
	dist1 <- vector(mode="numeric", length=0)
	if (any(j1)) {
		dist1 <- RegionT[j1,chromMid] - qmid
		if (qmid_orient== 1) {dist1 <- -dist1} #qmid_orient[0:expected to locate left,1:expect to locate right]
	}

	#case 2.
	j2 = which(RegionT$chrom == qchrom &
						 	qstart < RegionT$chromStart&
						 	RegionT$chromEnd <= qend)
	dist2 <- vector(mode="numeric", length=0)
	if (any(j2)) {
		dist2 <- RegionT[j2,chromMid] - qmid
		if (qmid_orient== 1) {dist2 <- -dist2}
	}
	
	#case 3.
	j3 = which(RegionT$chrom == qchrom &
						 	RegionT$chromStart < qend &
						 	qend < RegionT$chromEnd)
	dist3 <- vector(mode="numeric", length=0)
	if (any(j3)) {
		dist3 <- RegionT[j3,chromMid] - qmid
		if (qmid_orient== 1) {dist3 <- -dist3}
	}
	
	#case 4.
	j4 = which(RegionT$chrom == qchrom &
						 	RegionT$chromStart < qstart &
						 	qend <= RegionT$chromEnd)
	dist4 <- vector(mode="numeric", length=0)
	if (any(j4)) {
		dist4 <- RegionT[j4,chromMid] - qmid
		if (qmid_orient== 1) {dist4 <- -dist4}
	}
	
	if ('index' %in% colnames(RegionT)) {
		js = c(RegionT[j1]$index,RegionT[j2]$index,RegionT[j3]$index,RegionT[j4]$index)
	} else {
		js = c(j1,j2,j3,j4)
	}
	
	dists = c(dist1,dist2,dist3,dist4)
	if (any(js)){
		ui = which(!duplicated(js))
		dt = data.table(idx=js[ui],
										dist=dists[ui],
										adist=abs(dists[ui]))
		if (hit == 'best'){
			setorder(dt, -adist)
			return(dt[1]$idx)
		}
		else if (hit == 'all'){
			return(paste(dt$idx,collapse=","))
		}
	} else {
		return(NA)
	}
}

annot_distance_btwn_gr <-function(mygr){
  mygr <- sort(mygr)
  strands <- c('+','-')
  out_gr<-list()
  ext_bp <- 1e3
  for (i in 1:2){
    mygr2 <- mygr[as.character(strand(mygr))==strands[[i]],]
    L <- length(mygr2)
    
    #distance to preceding
    dist_btwn_region <- start(mygr2)[2:L] - end(mygr2)[1:L-1]
    
    mygr2$precede_dist <- append(0,dist_btwn_region)
    mygr2$follow_dist <- append(mygr2$precede_dist[2:L],ext_bp)

    out_gr[[i]] <- mygr2
  }
  
  return(unlist(GRangesList(out_gr)))
}

# --------------------------------------------------------------------

const_exp_delim <- function(my_gr,extbp){
  message('mid point, expand, same length of region ...')
  mid <- round(0.5*(start(my_gr)+end(my_gr)))
  start(my_gr) <- mid-extbp
  end(my_gr) <- mid+extbp
  return(my_gr)
}

make_gr_grid <- function(gr,window_bp=50,slide_bp=25,ncpu=1){
  library(parallel)
  if (T){
    split_grs <- mclapply(gr, function(r){
      GRanges(seqnames=seqnames(r),
              IRanges(start=seq(start(r),end(r),by=slide_bp),
                      width=window_bp),
              strand=strand(r))
    },mc.cores = 16)
  } else {
    split_grs <- lapply(gr, function(r){
      GRanges(seqnames=seqnames(r),
              IRanges(start=seq(start(r),end(r),by=slide_bp),
                      width=window_bp),
              strand=strand(r))
    })
  }
  return(GRangesList(split_grs))
}

get_phast_scores_surrounding_roi <- function(window_gr,phast,window_bp=50,slide_bp=25,ncpu=1){
  
  #browser()
  phast_grids <- make_gr_grid(window_gr,window_bp=window_bp,slide_bp=slide_bp,ncpu=ncpu)
  
  if (T){
    phast_scores <- mclapply(phast_grids,function(phast_grid){
      phastCons2 <- scores(phast,phast_grid)
      pcon_scores <- phastCons2$scores
      if (as.character(strand(phast_grid[1,])) == "-"){
        #browser()
        pcon_scores <- pcon_scores[length(phastCons2):1]
      }
      pcon_scores[is.na(pcon_scores)] <- 0.
      return(pcon_scores)
    },mc.cores = ncpu)
  } else {
    phast_scores <- lapply(phast_grids,function(phast_grid){
      phastCons2 <- scores(phast,phast_grid)
      pcon_scores <- phastCons2$scores
      if (as.character(strand(phast_grid[1,])) == "-"){
        #browser()
        pcon_scores <- pcon_scores[length(phastCons2):1]
      }
      pcon_scores[is.na(pcon_scores)] <- 0.
      return(pcon_scores)
    })
  }
  return(do.call(rbind, phast_scores))
}

get_phast_scores_roi <- function(window_gr,phast,window_bp=50,slide_bp=25,ncpu=1){
  
  phast_grids <- make_gr_grid(window_gr,window_bp,slide_bp)
  
  if (T){
    phast_scores <- mclapply(phast_grids,function(phast_grid){
      phastCons2 <- scores(phast,phast_grid)
      pcon_scores <- phastCons2$scores
      pcon_scores[is.na(pcon_scores)] <- 0.
      return(pcon_scores)
    },mc.cores = ncpu)
  } else {
    phast_scores <- lapply(phast_grids,function(phast_grid){
      phastCons2 <- scores(phast,phast_grid)
      pcon_scores <- phastCons2$scores
      pcon_scores[is.na(pcon_scores)] <- 0.
      return(pcon_scores)
    })
  }
  return(list(phastcon_score=phast_scores,
              window_gr=window_gr,
              window_bp=window_bp,
              slide_bp=slide_bp))
  # return(list(phastcon_score=do.call(rbind, phast_scores),
  #             window_gr=window_gr,
  #             window_bp=window_bp,
  #             slide_bp=slide_bp))
}

comp_movavg <- function(pas_gr, groi, edge_offset=0, wlen=50){

  message('generating sliding windows ...')
  groi_sw <- to_slide_windows_gr(groi,edge_offset,wlen=wlen)
  
  message('performing a moving average on CpG sites on unstranded prs templates ...')
  pas_cvg_gr <- coverage(pas_gr,weight="score")
  
  message('moving average ...')
  movavg <- binnedAverage(groi_sw, pas_cvg_gr, "binned_avg")
  return(movavg)
}

to_slide_windows_gr <- function(groi,edge_offset,wlen=50) {
  
  groi_pos <- groi
  
  strand(groi_pos) <- '+' #to extract only watson strand
  
  if (edge_offset>0){
    start(groi_pos) <- start(groi_pos)-edge_offset
    end(groi_pos) <- end(groi_pos)+edge_offset
  }
  
  message('generating binned windows ...')
  binned_windows <- slidingWindows(groi_pos, width=wlen, step=wlen/2)
  bwdts <- lapply(binned_windows,function(binned_window){
    as.data.table(binned_window)
  })
  
  message('merging binned wondiws ...')
  bwdt <- rbindlist(bwdts)
  bwdt$width <- NULL
  names(bwdt) <- c('chr','start','end','strand')
  prgr <- makeGRanges(bwdt,strand=F)
  
  return(prgr)
}

expand_gr <- function(my_gr,extbp5,extbp3){
  #we could use flank and union but more complicated. thus decided to write a simple method
  message('expanding goi in both 5prime and 3prime by strand sepcific manner ...')
  
  j <- which(as.character(strand(my_gr))=='+')
  start(my_gr)[j] <- start(my_gr)[j] - extbp5
  end(my_gr)[j] <- end(my_gr)[j] + extbp3
  
  j <- which(as.character(strand(my_gr))=='-')
  start(my_gr)[j] <- start(my_gr)[j] - extbp3
  end(my_gr)[j] <- end(my_gr)[j] + extbp5
  
  return(my_gr)
}

expand_reduced_gr <- function(goi_gr,extbp5=1e2,extbp3=1e2,do_reduce=F){

  prsE <- expand_gr(goi_gr,extbp5,extbp3)
  if (do_reduce){
    message('[WARNING] reducing gr and so metainfo will be removed by collapsing ...')
    prsE <- reduce(prsE,ignore.strand=FALSE)
  }
  prsE_dt <- as.data.table(prsE)
  setnames(prsE_dt,'seqnames','chr')
  prsE_dt$width <- NULL
  return(prsE_dt)
}


bed_genomecov <- function(gr,hg19_chrom_len_fn='/home/hongc2/projects/apa/activemotif_fastq/encode_analy_out/scripts/chipseq_integration_v1/data/genome.txt',wkd='/tmp') {
  
  #grs <- sortSeqlevels(c(gr1,gr2,gr3,gr4))
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)
  dt <- as.data.table(gr)
  dt$width <- NULL
  wkd <- '~/temp'
  
  temp_bed_fn <- file.path(wkd,'test.bed')
  fwrite(dt,file=temp_bed_fn,sep = '\t',col.names = FALSE)
  
  bgr_out_fn <- file.path(wkd,'test.bgr')
  
  cmd <- sprintf("bedtools genomecov -bga -split -i %s -g %s | awk '$4 != 0 { print $0 }' > %s",temp_bed_fn,hg19_chrom_len_fn,bgr_out_fn)
  
  system(cmd)
  
  dt <- fread(bgr_out_fn)
  colnames(dt) <- c('chr','start','end','ovr')
  dt$end <- dt$end - 1
  unlink(temp_bed_fn)
  unlink(bgr_out_fn)
  
  return(makeGRanges(dt))
}
