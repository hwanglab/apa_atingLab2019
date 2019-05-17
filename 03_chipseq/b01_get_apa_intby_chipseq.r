closeAllConnections()
rm(list=ls())

library(argparse)
library(data.table)
library(goldmine)
library(GenomicRanges)
library(openxlsx)

source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_chipseq.r'))

intersect_by_apa <- function(manorms,roi.gr) {
  # browser()

  apa_int_manorms <- list()
  
  chipseq_names <- names(manorms)
  manorm_total_bp <- 0
  apa_manorm_total_bp <- 0
  
  for (i in 1:length(manorms)){
    jmessage(chipseq_names[i])
    npeak <- manorms[[i]]
    
    #NOTE that roi.gr is query
    apa_int_manorm <- intersect(roi.gr,npeak,ignore.strand=T)
    apa_int_manorm$M_value <- 0
    apa_int_manorm$P_value <- 0
    apa_int_manorm$gidx <- 0
    
    #recover chipseq annot
    ij <- findOverlaps(apa_int_manorm,npeak,ignore.strand=T)
    apa_int_manorm[queryHits(ij)]$M_value <- npeak[subjectHits(ij)]$M_value
    apa_int_manorm[queryHits(ij)]$P_value <- npeak[subjectHits(ij)]$P_value
    
    #recover roi annot
    ij <- findOverlaps(apa_int_manorm,roi.gr,ignore.strand=T)
    apa_int_manorm[queryHits(ij)]$gidx <- roi.gr[subjectHits(ij)]$gidx
    apa_int_manorms[[chipseq_names[i]]] <- apa_int_manorm
    manorm_total_bp <- manorm_total_bp + sum(width(npeak))
    apa_manorm_total_bp <- apa_manorm_total_bp + sum(width(apa_int_manorm))
  }
  
  jmessage(sprintf('<rec>diff_chipseq_bp:%d',manorm_total_bp))
  jmessage(sprintf('<rec>apa_by_diff_chipseq_total_bp:%d',apa_manorm_total_bp))
  return(apa_int_manorms)
}

# =====================================
jmessage('input argument ........')
# goi_mode <- 'proximal'
reuse <- FALSE
manorm_rd <- "fastq/kundaje_encode/goi_apa_manorm.rd"
stopifnot(file.exists(manorm_rd))

exp_label <- 'chipseq'
apa_colon_xlsx <- "../01_polyAseq/01_wkd/out/04_AnnoApa/output/apa.ann_HCT116_vs_DKO.xlsx"
stopifnot(file.exists(apa_colon_xlsx))

macs2dir <- "fastq/bt2/HCT116_DKO_deseq2"
mbd_deseq2_fpath <- file.path(macs2dir,'deseq2.tsv.rc3_fc1.tsv')
stopifnot(file.exists(mbd_deseq2_fpath))

#define working directory
wkd <- 'fastq/cluster_analy'
if (!file.exists(wkd)) {dir.create(wkd)}

# =====================================

jmessage('loading diff chipseq binding sites analyzed by manorm')
load(manorm_rd) #manorms
manorms <- lapply(manorms,function(manorm){return(manorm$gr)})

# ------------------------------
jmessage('loading MBD raw signal sites ...')

mbd2_dt <- fread(mbd_deseq2_fpath)

mbd2_dt_ctrl <- mbd2_dt[,list(chrom,start,end,raw_control,pvalue)]
setnames(mbd2_dt_ctrl,c('chrom','raw_control','pvalue'),c('chr','M_value','P_value'))
manorms[['HCT116_Mbd']] <- makeGRanges(mbd2_dt_ctrl)

mbd2_dt_expr <- mbd2_dt[,list(chrom,start,end,raw_experiment,pvalue)]
setnames(mbd2_dt_expr,c('chrom','raw_experiment','pvalue'),c('chr','M_value','P_value'))
manorms[['DKO_Mbd']] <- makeGRanges(mbd2_dt_expr)

# --------------------------------------------------------------------
jmessage('loading APA genomic regions and choosing region of interest ...')
apa <- read_supp_tab3(apa_colon_xlsx,sheet_i=1) #keep apa to bring annotation back

jmessage('goi_mode=proximal')
roi <- apa[call=='ShorterInB',list(chr,start,end,gidx)]

roi.gr <- makeGRanges(roi)
roi.gr <- sortSeqlevels(roi.gr)
roi.gr <- sort(roi.gr)
jmessage(sprintf('<rec>total number of regions in roi.gr[%d]',length(roi.gr)))
jmessage(sprintf('<rec>total number of basepairs in roi.gr[%d]',sum(width(roi.gr))))
# -----------------------------------------
jmessage('intersect apa by chipseq np + (manorm)')

by_apa_rd_fn <- file.path(wkd,sprintf('%s.rd',exp_label))

if (reuse & file.exists(by_apa_rd_fn)) {
  load(by_apa_rd_fn)
} else {
  apa_int_manorms <- intersect_by_apa(manorms,roi.gr)
  save(apa_int_manorms,roi,file=by_apa_rd_fn,compress = T)
}
jmessage(sprintf('continue to work on [%s]',by_apa_rd_fn))
