library(data.table)
library(GenomicRanges)
library(IRanges)
library(goldmine)

get_region_btwn_apa <- function(apa.calls) {
  # browser()
  message("obtain genomic region between two PRs ...")
  # Shorts
  proximalp <- apa.calls[mcall=="ShorterInB" & strand=='+', list(chr=chr,start=up_end,end=down_start,strand=strand)]
  proximalm <- apa.calls[mcall=="ShorterInB" & strand=='-', list(chr=chr,start=down_end,end=up_start,strand=strand)]
  proximal_btwn <- makeGRanges(rbindlist(list(proximalp,proximalm)),strand = T)
  
  # Longers
  distalp <- apa.calls[mcall=="LongerInB" & strand=='+', list(chr=chr,start=down_end,end=up_start,strand=strand)]
  distalm <- apa.calls[mcall=="LongerInB" & strand=='-', list(chr=chr,start=up_end,end=down_start,strand=strand)]
  distal_btwn <- makeGRanges(rbindlist(list(distalp,distalm)),strand = T)
  return(list(proximal=proximal_btwn,distal=distal_btwn))
}

get_region_btwn_apa_mid <- function(apa.calls) {
  # browser()
  message("obtain genomic region between two PRs ...")
  # Shorts
  proximalp <- apa.calls[mcall=="ShorterInB" & strand=='+', list(chr=chr,start=round((up_start+up_end)/2),end=round((down_start+down_end)/2),strand=strand)]
  proximalm <- apa.calls[mcall=="ShorterInB" & strand=='-', list(chr=chr,start=round((down_start+down_end)/2),end=round((up_start+up_end)/2),strand=strand)]
  proximal_btwn <- makeGRanges(rbindlist(list(proximalp,proximalm)),strand = T)
  
  # Longers
  distalp <- apa.calls[mcall=="LongerInB" & strand=='+', list(chr=chr,start=round((down_start+down_end)/2),end=round((up_start+up_end)/2),strand=strand)]
  distalm <- apa.calls[mcall=="LongerInB" & strand=='-', list(chr=chr,start=round((up_start+up_end)/2),end=round((down_start+down_end)/2),strand=strand)]
  distal_btwn <- makeGRanges(rbindlist(list(distalp,distalm)),strand = T)
  return(list(proximal=proximal_btwn,distal=distal_btwn))
}

callMUD <- function(apa)
{
  
  # wt mean method
  calls <- apa[,list(gene_name=gene_name[1],chr=chr[1],strand=strand[1],wa=weighted.mean(x=end,w=mean_a_frac),wb=weighted.mean(x=end,w=mean_b_frac)),by="tu"]
  calls$wb_minus_wa <- calls$wb-calls$wa
  
  calls <- calls[!is.nan(wa)&!is.nan(wb),]
  # Call direction based on strand and value
  cp <- calls[strand=="+",]
  cp$wcall <- ""
  cp[wb>=wa,]$wcall <- "LongerInB"
  cp[wb<wa,]$wcall <- "ShorterInB"
  # browser()
  
  stopifnot(cp$wcall!="")
  cm <- calls[strand=="-",]
  cm$wcall <- ""
  cm[wb>=wa,]$wcall <- "ShorterInB"
  cm[wb<wa,]$wcall <- "LongerInB"
  
  stopifnot(cm$wcall!="")
  calls <- rbind(cp,cm)
  #calls[wcall=="",]$wcall <- "Error"
  
  # mud method
  ud <- apa[,list(MostUp=pr[b_minus_a_frac==max(b_minus_a_frac)],MostDown=pr[b_minus_a_frac==min(b_minus_a_frac)]),by="tu"]
  stopifnot(nrow(ud)==length(unique(apa$tu)))
  
  ud$chr <- apa[match(ud$MostUp,apa$pr),]$chr
  ud$up_start <- apa[match(ud$MostUp,apa$pr),]$start
  ud$up_end <- apa[match(ud$MostUp,apa$pr),]$end
  ud$down_start <- apa[match(ud$MostDown,apa$pr),]$start
  ud$down_end <- apa[match(ud$MostDown,apa$pr),]$end
  ud$strand <- apa[match(ud$MostDown,apa$pr),]$strand
  
  upp <- ud[strand=="+",]
  upp$mcall <- ""
  upp[up_start<=down_start,]$mcall <- "ShorterInB"
  upp[up_start>down_start,]$mcall <- "LongerInB"
  
  stopifnot(upp$mcall!="")
  
  upm <- ud[strand=="-",]
  upm$mcall <- ""
  upm[up_start<=down_start,]$mcall <- "LongerInB"
  upm[up_start>down_start,]$mcall <- "ShorterInB"
  
  stopifnot(upm$mcall!="")
  
  udmeth <- rbind(upp,upm)
  
  calls <- cbind(calls[match(udmeth$tu,calls$tu),],udmeth[,c("MostUp","MostDown","up_start","up_end","down_start","down_end","mcall"),with=F])
  calls$agree <- calls$wcall==calls$mcall
  calls
}

# --------------------------------------------------------------------

split_apa_prs <- function(mycalls,extbp=1e4){
  prpu <- mycalls[strand=="+",
                  list(tu=tu,strand=strand,
                       chr=chr,
                       start=round((up_start+up_end)/2)-extbp,
                       end=round((up_start+up_end)/2)+extbp,
                       usage='up',
                       pr=MostUp)]
  
  prpd <- mycalls[strand=="+",
                  list(tu=tu,strand=strand,
                       chr=chr,
                       start=round((down_start+down_end)/2)-extbp,
                       end=round((down_start+down_end)/2)+extbp,
                       usage='down',
                       pr=MostDown)]
  
  prmu <- mycalls[strand=="-",
                  list(tu=tu,strand=strand,
                       chr=chr,
                       start=round((up_start+up_end)/2)-extbp,
                       end=round((up_start+up_end)/2)+extbp,
                       usage='up',
                       pr=MostUp)]
  
  prmd <- mycalls[strand=="-",
                  list(tu=tu,strand=strand,
                       chr=chr,
                       start=round((down_start+down_end)/2)-extbp,
                       end=round((down_start+down_end)/2)+extbp,
                       usage='down',
                       pr=MostDown)]
  
  pr <- rbind(prpu,prpd,prmu,prmd)
}

# --------------------------------------------------------------------
load_region_per_type_surrounding_by_pr <- function(apa.calls,extbp=1e4){
  message('loading PR locations from apa.calls ...')
  
  # Annotate to class out the APA as intronic or UTR
  # For the UTR type, want them to both be after the inner-most cdsEnd
  
  # Split PR usages
  # Shorts
  mycalls <- apa.calls[mcall=="ShorterInB",]
  pr <- split_apa_prs(mycalls,extbp)
  pr.short.gr <- makeGRanges(pr,strand=TRUE)
  
  # Longs
  mycalls <- apa.calls[mcall=="LongerInB",]
  pr <- split_apa_prs(mycalls,extbp)
  pr.long.gr <- makeGRanges(pr,strand=TRUE)
  
  return(list(shorter=pr.short.gr,longer=pr.long.gr))
}

get_nonoverlap_btwn_gr <-function(mygr,max_extbp=1.5e3){
  
  #browser()
  
  mygr <- sort(mygr)
  strands <- c('+','-')
  
  novregions <- list()
  
  for (i in 1:2){
    mygr2 <- mygr[as.character(strand(mygr))==strands[[i]],]
    mygr2_by_chrom <- splitAsList(mygr2, seqnames(mygr2))
    browser()
    novregions[[i]] <- lapply(mygr2_by_chrom,function(mygr2c){
      L <- length(mygr2c)
      st1 <- start(mygr2c)[2:L]
      dist_precede <- start(mygr2c)[2:L] - end(mygr2c)[1:L-1]
      mid2 <- start(mygr2c)[2:L]+round(dist_precede/2.)
      dist1 <- 2
      if (start(mygr2c)[1] >= max_extbp) {
        dist1 <- start(mygr2c)[1] - max_extbp
      }
      dist_precede <- append(dist1,dist_precede)
      pos_dist_idx <- which(dist_precede>0)
      mygr2 <- mygr2c[pos_dist_idx>0,]
      
      mid2 <- append(dist1,mid2)
      novregion <- GRanges(seqnames = seqnames(mygr2c),
                                IRanges(start=mid2[1:L-1]-1,
                                        end=mid2[2:L])
      )
    })
    
    #browser()
  }
  
  x <- 1
  #browser()
  return(unlist(GRangesList(novregion)))
}

non_overlap_exp_prs <-function(mygr,max_extbp=1.5e3) {
  
  # browser()
  
  mygr <- sort(mygr)
  strands <- c('+','-')
  exp_grs <- list()
  
  for (i in 1:2){ #stran
    mygr2 <- mygr[as.character(strand(mygr))==strands[[i]],]
    mygr2_by_chrom <- splitAsList(mygr2, seqnames(mygr2))
    
    exp_gr_chrs <- list()
    j1 <- 1
    for (j in 1:length(mygr2_by_chrom)){
      
      mygr2c <- mygr2_by_chrom[[j]]

      exit_flag <- F
      while(!exit_flag){
        L <- length(mygr2c)
        dist_precede <- start(mygr2c)[2:L] - end(mygr2c)[1:L-1]
        mygr2c <- mygr2c[dist_precede>0,]
        if (!any(dist_precede<0)){exit_flag<-T}
      }
      
      if (L>1){
        # add the first delimiter
        # browser()
        mid2 <- end(mygr2c)[1:L-1]+round(dist_precede/2.)
        st_pos <- 2
        
        # message(sprintf('i=%d,j=%d,count=%d',i,j,length(mygr2c)))
        
        if (start(mygr2c)[1] >= max_extbp) {
          st_pos <- start(mygr2c)[1] - max_extbp + 1
        }
        
        # -------------------
        # determine start delimiter
        mid2st <- append(st_pos,mid2)
        st1 <- start(mygr2c) - max_extbp
        
        adj_starts <- sapply(1:L,function(m){
          if (st1[m] > mid2st[m]){
            new_st <- st1[m]
          } else {
            new_st <- mid2st[m]
          }
          return(new_st)
        })
        
        # -------------------
        # determine end delimiter
        ed_pos <- end(mygr2c)[L] + 3e3
        mid2ed <- append(mid2,ed_pos)
        ed1 <- end(mygr2c) + max_extbp
        
        adj_ends <- sapply(1:L,function(m){
          if (ed1[m] < mid2ed[m]){
            new_ed <- ed1[m]
          } else {
            new_ed <- mid2ed[m]
          }
          return(new_ed)
        })
        exp_gr_chrs[[j1]] <- GRanges(seqnames = seqnames(mygr2c),
                                     IRanges(start=adj_starts,
                                             end=adj_ends),
                                     strand=strands[[i]],
                                     tu=mygr2c$tu,
                                     pr=mygr2c$pr)
        j1 = j1 + 1
      }
    }
    # browser()
    exp_grs[[i]] <- unlist(GRangesList(exp_gr_chrs))
  }
  
  # x <- 1
  # browser()
  return(unlist(GRangesList(exp_grs)))
}

to_slide_windows_gr2 <- function(groi,edge_offset) {
  
  groi_pos <- groi
  
  strand(groi_pos) <- '+' #to extract only watson strand
  
  if (edge_offset>0){
    start(groi_pos) <- start(groi_pos)-edge_offset
    end(groi_pos) <- end(groi_pos)+edge_offset
  }
  
  message('generating binned windows ...')
  binned_windows <- slidingWindows(groi_pos, width=50, step=25)
  M <- length(binned_windows)
  bwdts <- lapply(1:M,function(m){
    dt <- as.data.table(binned_windows[[m]])
    R <- dim(dt)[1]
    #to make odd number
    if ((R>1) & (R%%2 ==1)) {
      dt <- dt[1:(R-1),]
    }
    dt$groi_idx <- sprintf('g%d',m)
    return(dt)
  })
  
  message('merging binned wondiws ...')
  bwdt <- rbindlist(bwdts)
  
  bwdt$width <- NULL
  names(bwdt) <- c('chr','start','end','strand','groi_idx')
  prgr <- makeGRanges(bwdt,strand=F)
  
  return(prgr)
}

comp_movavg2 <- function(signal_gr, groi, edge_offset=0){
  
  message('generating sliding windows ...')
  groi_sw <- to_slide_windows_gr2(groi,edge_offset)
  
  message('performing a moving average of the signal on goi_gr on unstranded prs templates ...')
  signal_cvg_gr <- coverage(signal_gr,weight="score")
  
  message('moving average ...')
  movavg <- binnedAverage(groi_sw, signal_cvg_gr, "binned_avg")
  return(movavg)
}

map_refgene_to_ensembl <- function(prepropa_rd,apa_ann_rd,debug=F) {
  
  #prepropa_rd <- file.path(mount_prefix,'apa_atingLab2019','01_polyAseq','01_wkd','out','03_CallApa','output','prepropa.rd')
  #apa_ann_rd <- file.path(mount_prefix,'apa_atingLab2019','01_polyAseq','01_wkd','out','04_AnnoApa','output','apa.ann.rd')
  
  if (debug){browser()}
  message('filtering out genes not of our interest!')
  
  if (file.exists(prepropa_rd)) {
    load(prepropa_rd)
  } else {
    stop(sprintf('[%s] not exist!',prepropa_rd))
  }
  
  if (file.exists(apa_ann_rd)) {
    load(apa_ann_rd)
  } else {
    stop(sprintf('[%s] not exist!',apa_ann_rd))
  }

  rm(apa.int,jtu,pr,reads,red_ens)
  
  apa.calls$ensembl_gene <- genes[match(apa.calls$gene_name,genes$name),gene.id]
  goi_map <- apa.calls[,list(tu,gene_name,strand,ensembl_gene)]
  
  updated <- data.table(gene_name=c('ZNF595','OTUD7B','SMEK2','DDX52','PPOX','RP11-103J17.2','RP11-30J20.1',"RAB26","MT1A","AC011526.1"),
                        ensembl_gene=c('ENSG00000272602','ENSG00000264522','ENSG00000275052','ENSG00000278053','ENSG00000143224','ENSG00000261761','ENSG00000254101','ENSG00000167964','ENSG00000205362','ENSG00000267107'))
  
  upd_entries <- updated[match(goi_map$gene_name,updated$gene_name),ensembl_gene]
  j <- which(!is.na(upd_entries))
  goi_map[j,ensembl_gene:=upd_entries[j]]
  return(goi_map)
}

get_tagMatrixList <- function(asignal,window_grs,bgnd1,bgnd2_sampled,extbp,rd_fn,reuse=T){
  
  if (reuse & file.exists(rd_fn)){
    message('use prev rd file ...')
    load(rd_fn)
  } else {
    # ==================
    load(file=sprintf('./rds/03_stratify_apa_prs_e%d.rd',extbp)) #window_grs
    load(file=sprintf('./rds/04_background_prs_e%d.rd',extbp)) # extbp,bgnd1,bgnd2
    
    message('2. get tag matrix for background ...')
    backgrnd1Matrix <- getTagMatrix(asignal,
                                    weightCol = 'binned_avg',
                                    windows = bgnd1)
    
    backgrnd2Matrix <- getTagMatrix(asignal,
                                    weightCol = 'binned_avg',
                                    windows = bgnd2_sampled)
    
    message('3. calling getTagMatrix() ...')

    L <- length(window_grs)
    tagMatrixList <- list()
    for (i in 1:L) {
      
      tagMatrix <- getTagMatrix(asignal,
                                weightCol = 'binned_avg',
                                windows = window_grs[[i]])
      tagMatrixList[[i]] <- sort_matrix_by_rowsum(tagMatrix)
    }
    
    tagMatrixList[[L+1]] <- sort_matrix_by_rowsum(backgrnd1Matrix)
    tagMatrixList[[L+2]] <- sort_matrix_by_rowsum(backgrnd2Matrix)
    save(tagMatrixList,file=rd_fn,compression=T)
  }
  return(tagMatrixList)
}

get_heatmaps <- function(tagMatrixList,axis_args2,xlab2,L){

  source(file.path(Sys.getenv('R_UTIL_APA'),'lib_util.R'))
  maxVal <- get_max_val_in_list(tagMatrixList)
  if (maxVal>0){
    minVal <- get_min_val_in_list(tagMatrixList)/maxVal
  } else {
    minVal <- 0.
  }
  
  goi_names <- names(tagMatrixList)
  extbp <- -axis_args2$at[1]
  tagHeatmap2(tagMatrixList[1:L],
              xlim=c(-extbp, extbp),
              zlim=c(minVal,1.),
              axis_args=axis_args2,
              color="red",
              xlab=xlab2,
              title=goi_names[1:L],
              denom=maxVal)
  
  tagHeatmap2(tagMatrixList[(L+1):(L+2)],
              xlim=c(-extbp, extbp),
              zlim=c(minVal,1.),
              axis_args=axis_args2,
              color="red",
              xlab=xlab2,
              title=goi_names[(L+1):(L+2)],
              denom=maxVal)
}

tagHeatmap2 <- function(tagMatrix, xlim, zlim, axis_args, xlab="", ylab="", title=NULL, color="red", denom=1.) {
  listFlag <- FALSE
  if (is(tagMatrix, "list")) {
    listFlag <- TRUE
  }
  peakHeatmap.internal3(tagMatrix, xlim, axis_args, listFlag, color, xlab, ylab, title, denom, zlim)
}

getCols2 <- function(n) {
  col <- c("#8dd3c7", "#ffffb3", "#bebada",
           "#fb8072", "#80b1d3", "#fdb462",
           "#b3de69", "#fccde5", "#d9d9d9",
           "#bc80bd", "#ccebc5", "#ffed6f")
  
  col2 <- c("#1f78b4", "#ffff33", "#c2a5cf",
            "#ff7f00", "#810f7c", "#a6cee3",
            "#006d2c", "#4d4d4d", "#8c510a",
            "#d73027", "#78c679", "#7f0000",
            "#41b6c4", "#e7298a", "#54278f")
  
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
            "#33a02c", "#fb9a99", "#e31a1c",
            "#fdbf6f", "#ff7f00", "#cab2d6",
            "#6a3d9a", "#ffff99", "#b15928")
  
  ## colorRampPalette(brewer.pal(12, "Set3"))(n)
  colorRampPalette(col3)(n)
}

peakHeatmap.internal3 <- function(tagMatrix, xlim, axis_args, listFlag, color, xlab, ylab, title, denom, zlim) {
  if ( is.null(xlab) || is.na(xlab))
    xlab <- ""
  if ( is.null(ylab) || is.na(ylab))
    ylab <- ""
  
  if (listFlag) {
    nc <- length(tagMatrix)
    if ( is.null(color) || is.na(color) ) {
      cols <- getCols2(nc)
    } else if (length(color) != nc) {
      cols <- rep(color[1], nc)
    } else {
      cols <- color
    }
    
    if (is.null(title) || is.na(title))
      title <- names(tagMatrix)
    if (length(xlab) != nc) {
      xlab <- rep(xlab[1], nc)
    }
    if (length(ylab) != nc) {
      ylab <- rep(ylab[1], nc)
    }
    if (length(title) != nc) {
      title <- rep(title[1], nc)
    }
    par(mfrow=c(1, nc))
    for (i in 1:nc) {
      peakHeatmap.internal4(tagMatrix[[i]], xlim, axis_args, cols[i], xlab[i], ylab[i], title[i], denom, zlim)
    }
  } else {
    if (is.null(color) || is.na(color))
      color <- "red"
    if (is.null(title) || is.na(title))
      title <- ""
    
    title <- sprintf("%s(total:%d)",dim(tagMatrix)[2])
    peakHeatmap.internal4(tagMatrix, xlim, axis_args, color, xlab, ylab, title, denom, zlim)
  }
}

peakHeatmap.internal4 <- function(tagMatrix, xlim=NULL, axis_args=NA, color="red", xlab="", ylab="", title="", denom=1., zlim=NULL) {
  
  tagMatrix <- t(apply(tagMatrix, 1, function(x) x/denom))
  ii <- order(rowSums(tagMatrix))
  tagMatrix <- tagMatrix[ii,]
  cols <- colorRampPalette(c("white",color))(200)
  if (is.null(xlim)) {
    xlim <- 1:ncol(tagMatrix)
  } else if (length(xlim) == 2) {
    xlim <- seq(xlim[1], xlim[2])
  }
  #axis_args <- c(-5e3,-2.5e3, 0, 2.5e3, 5e3)
  if (is.null(zlim)){
    image.plot(x=xlim, y=1:nrow(tagMatrix),z=t(tagMatrix),useRaster=TRUE, col=cols, yaxt="n", ylab="", xaxt="n", xlab="", main=title, horizontal=FALSE)
  } else {
    image.plot(x=xlim, y=1:nrow(tagMatrix),z=t(tagMatrix),zlim=zlim,useRaster=TRUE, col=cols, yaxt="n", ylab="", xaxt="n", xlab="", main=title, horizontal=FALSE)
  }
  
  axis(side = 1,at=axis_args$at,labels=axis_args$labels)
  
}


read_st3_sheet_to_gr <- function(most_updn_pr_xlsx,sheet,btn_flank=100){
  mud_dt <- as.data.table(read.xlsx(most_updn_pr_xlsx,sheet=sheet))
  
  # ----------------------------------------
  jmessage('defining extend apa region ...')
  mud_dt[,pr_up_mid:=round((pr_up_end+pr_up_start)/2)]
  mud_dt[,pr_down_mid:=round((pr_down_end+pr_down_start)/2)]
  mud_dt[,btn_width:=abs(pr_up_mid-pr_down_mid)]
  
  mud_dt$pr_start <- as.numeric()
  mud_dt$pr_end <- as.numeric()
  
  mud_dt[pr_up_mid>=pr_down_mid,pr_end:=(pr_up_mid+btn_flank)]
  mud_dt[pr_up_mid<pr_down_mid,pr_start:=(pr_up_mid-btn_flank)]
  
  mud_dt[pr_down_mid>pr_up_mid,pr_end:=(pr_down_mid+btn_flank)]
  mud_dt[pr_down_mid<=pr_up_mid,pr_start:=(pr_down_mid-btn_flank)]
  
  return(mud_dt)
}

read_supp_tab3 <- function(st3_xlsx,sheet_i=2) {
  
  mud_dt <- read_st3_sheet_to_gr(st3_xlsx,sheet=sheet_i)
  mud_dt$gidx <- 1:dim(mud_dt)[1]
  mud <- mud_dt[,list(tu_chr,pr_start,pr_end,tu_strand,gidx,gene_name,call)]
  colnames(mud) <- c('chr','start','end','strand','gidx','gene_name','call')
  mud$width <- mud$end - mud$start + 1
  
  return(mud)
}
