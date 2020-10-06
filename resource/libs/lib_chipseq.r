source(sprintf('%s/lib_util.R',Sys.getenv('R_UTIL_APA')))
library(data.table)
library(openxlsx)
library(goldmine)
library(stringr)
library(gsubfn)

mount_prefix <- get_mount_dir()

source(file.path(Sys.getenv('R_UTIL_APA'),'lib_PRs.r'))

get_chipseq_bam <- function(mount_prefix){
  
  bamd <- file.path(mount_prefix,'03_chipseq/fastq/kundaje_encode/bamlink')
  
  fbase <- c("01_0307_00CWCCF_HCT116_AHT_CTCF_hs_i72.nodup.bam",
             "02_0308_00CWCCF_DKO_AHT_CTCF_hs_i73.nodup.bam",
             "13_031N_00CWCCF_HCT116_AHT_RAD21_hs_i72.nodup.bam",
             "14_031O_00CWCCF_DKO_AHT_RAD21_hs_i73.nodup.bam",
             "19_030B_00CWCCF_HCT116_AHT_H3K27Ac_hs_i79.nodup.bam",
             "20_030C_00CWCCF_DKO_AHT_H3K27Ac_hs_i80.nodup.bam",
             "25_031J_00CWCCF_HCT116_AHT_Pol2Ser5_hs_i68.nodup.bam",
             "26_031K_00CWCCF_DKO_AHT_Pol2Ser5_hs_i69.nodup.bam")
  
  sample <-c('HC','DC','HR', 'DR', 'HH', 'DH', 'HPS5', 'DPS5')
  
  color2 <- c('gray','black',
              'gray','black',
              'gray','black',
              'gray','black')
  
  bam_dt <- data.table(fpath=file.path(bamd,fbase),
                       sample=sample,
                       color2=color2)
  return(bam_dt)
}


get_chipseq_bam_full <- function(mount_prefix){
  
  
  bamd <- file.path(mount_prefix,'03_chipseq/fastq/kundaje_encode/bamlink')
  
  fbase <- c("01_0307_00CWCCF_HCT116_AHT_CTCF_hs_i72.nodup.bam",
             "02_0308_00CWCCF_DKO_AHT_CTCF_hs_i73.nodup.bam",
             "07_030P_00CWCCF_HCT116_AHT_SMC1_hs_i86.nodup.bam",
             "08_030Q_00CWCCF_DKO_AHT_SMC1_hs_i87.nodup.bam",
             "13_031N_00CWCCF_HCT116_AHT_RAD21_hs_i72.nodup.bam",
             "14_031O_00CWCCF_DKO_AHT_RAD21_hs_i73.nodup.bam",
             "19_030B_00CWCCF_HCT116_AHT_H3K27Ac_hs_i79.nodup.bam",
             "20_030C_00CWCCF_DKO_AHT_H3K27Ac_hs_i80.nodup.bam",
             "25_031J_00CWCCF_HCT116_AHT_Pol2Ser5_hs_i68.nodup.bam",
             "26_031K_00CWCCF_DKO_AHT_Pol2Ser5_hs_i69.nodup.bam",
             "31_031F_00CWCCF_HCT116_AHT_Pol2Ser2_hs_i93.nodup.bam",
             "32_031G_00CWCCF_DKO_AHT_Pol2Ser2_hs_i94.nodup.bam")
  
  samples <-c('HC','DC','HS','DS','HR', 'DR', 'HH', 'DH', 'HPS5', 'DPS5','HPS2', 'DPS2')
  groups <- rep(c('HCT116','DKO'),6)
  proteins <- c(rep('CTCF',2),rep('SMC1',2),rep('RAD21',2),rep('H3K27Ac',2),rep('Pol2Ser5',2),rep('Pol2Ser2',2))
  
  color2 <- c('gray','black',
              'gray','black',
              'gray','black',
              'gray','black',
              'gray','black',
              'gray','black')
  
  bam_dt <- data.table(fpath=file.path(bamd,fbase),
                       group=groups,
                       protein=proteins,
                       sample=samples,
                       color2=color2)
  return(bam_dt)
}

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

to_slide_windows_gr3 <- function(groi,window=50,step=25,edge_offset=0) {
  # browser()
  groi2 <- groi
  
  if (edge_offset>0){
    start(groi2) <- start(groi2)-edge_offset
    end(groi2) <- end(groi2)+edge_offset
  }
  
  message('generating binned windows ...')
  binned_windows <- slidingWindows(groi2, width=window, step=step)
  M <- length(binned_windows)
  bwdts <- lapply(1:M,function(m){
    dt <- as.data.table(binned_windows[[m]])
    R <- dim(dt)[1]
    #to make odd number
    if ((R>1) & (R%%2 ==1)) {
      dt <- dt[1:(R-1),]
    }
    dt$gidx <- groi$gidx[m]
    return(dt)
  })
  
  message('merging binned wondiws ...')
  bwdt <- rbindlist(bwdts)
  
  bwdt$width <- NULL
  names(bwdt) <- c('chr','start','end','strand','gidx')
  slide_groi <- makeGRanges(bwdt,strand=T)
  
  slide_groi <- sortSeqlevels(slide_groi)
  slide_groi <- sort(slide_groi)
  
  return(slide_groi)
}

varv_map2_unifiedv <- function(binned_avg,umax_bp){
  
  # browser()
  
  unifiedv <- rep(NA,(2*umax_bp+1))
  umidx <- umax_bp + 1
  
  R <- length(binned_avg)
  if (R %% 2 == 0) {
    if (R < (2*umax_bp+1)) {
      binned_avg[R+1] <- binned_avg[R]
      R <- R + 1
    } else {
      binned_avg <- binned_avg[1:(R-1)]
      R <- R - 1
    }
  }
  
  mmidx <- round((R+1)/2)
  mmax_bp <- mmidx - 1
  
  if (umidx>=mmidx) {
    offset <- umax_bp - mmax_bp + 1
    unifiedv[offset:(offset+R-1)] <- binned_avg
    
    # if ((offset-1)>0) {
    #   unifiedv[offset-1] <- 1.
    # } else {
    #   unifiedv[1] <- 1.
    # }
    
    # if ((offset+R)<=(2*umax_bp+1)) {
    #   unifiedv[(offset+R)] <- 1.
    # } else {
    #   unifiedv[(2*umax_bp+1)] <- 1.
    # }
    
  } else {
    message('umidx < mmidx')
    stopifnot(F)
  }
  return(unifiedv)
}

map_to_unifmat <- function (mysignal,apa_btwn,max_width){
  
  # browser()
  mysignal_by_gene_list <- splitAsList(mysignal, mysignal$gidx)
  unif_densv_list <- lapply(mysignal_by_gene_list,function(sw_gr){
    unif_densv <- varv_map2_unifiedv(sw_gr$binned_avg,max_width)
  })
  
  if (F) {
    #debug
    for (i in 1:length(unif_densv_list)) {
      message(sprintf("i=%d,dim[1]=%d",i,length(unif_densv_list[[i]])))
    }
  }
  
  unifmt <- do.call(rbind, unif_densv_list)
  
  signal_gidx <- as.integer(names(mysignal_by_gene_list))
  unifmt <- unifmt[match(apa_btwn$gidx,signal_gidx),]
  
  # rownames(unifmat) <- sprintf('%s:%d-%d',
  #                              seqnames(apa_btwn)[j],
  #                              start(apa_btwn)[j],
  #                              end(apa_btwn)[j])
  
  rownames(unifmt) <- apa_btwn$gene_name
  colnames(unifmt) <- -max_width:max_width
  return(unifmt)
}


get_encode_narrow_peak_fpath <- function() {
  mount_prefix <- get_mount_dir()
  
  hct_np_fpat <- sprintf("%s/03_chipseq/fastq/kundaje_encode/*HCT116*/peak/macs2/rep1/*_HCT116*.nodup*.filt.narrowPeak.gz",mount_prefix)

  dko_np_fpat <- sprintf("%s/03_chipseq/fastq/kundaje_encode/*DKO*/peak/macs2/rep1/*_DKO*.nodup*.filt.narrowPeak.gz",mount_prefix)
  
  fpats <- data.table(npeak_fn=c(hct_np_fpat,dko_np_fpat),
                      sample=c('HCT116','DKO'))
  
  comp_exps <- lapply(fpats$npeak_fn, function(npeak_fn){
    cmd <- sprintf('ls %s',npeak_fn)
    fpaths <- system(cmd,intern = TRUE)
    
    proteins <- list()
    for (i in 1:length(fpaths)) {
      proteins[[i]] <- str_extract(fpaths[[i]],"(?<=AHT_)\\w+(?=_hs)") #r regexp extract
    }
    return(data.table(fpath=fpaths,name=proteins))
  })
  
  comp_exps[[1]]$group <- fpats$sample[1]
  comp_exps[[2]]$group <- fpats$sample[2]
  comp_exp <- rbindlist(comp_exps)
  rm(comp_exps)
  
  return(comp_exp)
}

readSampleInfo <- function(file,colors=NULL)
{
  if(is.null(file)){stop("Must provide a path to a file to open.")}
  if(!file.exists(file)){stop(paste("Could not open: ",file,sep=""))}
  samplesheet <- read.csv(file, header=TRUE, stringsAsFactors=FALSE)
  if(sum(c("sample","group","bam") %in% colnames(samplesheet))!=3){stop("CSV must contain columns: sample, group, bam")}
  nsamp <- nrow(samplesheet)
  ngroup <- length(unique(samplesheet$group))
  message(paste("Found ", nsamp, " samples in ", ngroup, " groups",sep=""))
  grouporder <- unique(samplesheet$group)
  samplesheet$group <- factor(samplesheet$group,levels=grouporder,ordered=TRUE)
  message(paste("Group order detected as: ", paste(grouporder,collapse=", ", sep=""), sep=""))
  message("Number of samples in each group:")
  print(summary(samplesheet$group))
  samplesheet <- samplesheet[order(samplesheet$group,samplesheet$sample),]
  if(is.null(colors) & sum(colnames(samplesheet)=="color")==0)
  {
    message("Auto-picking group colors from RColorBrewer")
    if(ngroup<=9)
    {
      colors <- RColorBrewer::brewer.pal(ngroup,"Set1")
    } else
    {
      colors <- colorRampPalette(brewer.pal(9,"Set1"))(ngroup)
    }
    samplesheet$color <- colors[as.numeric(samplesheet$group)]
  } else if(!is.null(colors))
  {
    if(length(colors)!=ngroup){stop("colors vector must equal number of groups in file")}
    samplesheet$color <- colors[as.numeric(samplesheet$group)]
  }
  samplesheet
}
