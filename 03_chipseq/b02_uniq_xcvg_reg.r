library(argparse)
library(data.table)
library(goldmine)
library(GenomicRanges)
library(parallel)

source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_bed.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_chipseq.r'))

if (F) {
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  
  parser <- ArgumentParser(description='prepropa')
  parser$add_argument("-b", "--bpco", type="integer", required=FALSE,
                      dest="bpco",default=10,
                      help="minimum genomic region length in bp[10]")
  parser$add_argument("-c", "--chrom_len_fpath", type="character", required=FALSE,
                      dest="chrom_len_fpath",default="../resource/ref/hg19_chrom_length.txt",
                      help="chromosome length file path")
  parser$add_argument("-r", "--reuse", type="integer", required=FALSE,
                      dest="reuse",default=0,
                      help="Yes:1, No:[0]")
  
  args <- parser$parse_args()
} else {
  args <- data.table(reuse = 0,
                     bpco = 10,
                     chrom_len_fpath="../resource/ref/hg19_chrom_length.txt")
}

if (args$reuse==0) {
  reuse <- FALSE
} else {
  reuse <- TRUE
}

exp_label <- 'chipseq'
wkd <- 'fastq/cluster_analy'
by_apa_rd_fn <- file.path(wkd,sprintf('%s.rd',exp_label))
message(by_apa_rd_fn)
load(by_apa_rd_fn) #roi,apa_int_manorms

exp_label <- sprintf('%s_xing',exp_label)
ovr_apa_rd_fn <- file.path(wkd,sprintf('%s.rd',exp_label))
if (reuse & file.exists(ovr_apa_rd_fn)) {
  load(ovr_apa_rd_fn) #cvg_dt
} else {
  xing_diff_gr <- apa_int_manorms[[1]]
  message(sprintf('|xdg|=%d',length(xing_diff_gr)))
  proteins <- names(apa_int_manorms)
  for (i in 2:length(apa_int_manorms)) {
    xing_diff_gr <- c(xing_diff_gr,apa_int_manorms[[i]])
  }
  
  xing_diff_gr <- sortSeqlevels(xing_diff_gr)
  xing_diff_gr <- sort(xing_diff_gr)
  
  start(xing_diff_gr) <- start(xing_diff_gr) - 1 #start bp in BED is 0-based
  
  message('running bed_genomecov ...')
  xing_diff_gr <- bed_genomecov(xing_diff_gr,
                                hg19_chrom_len_fn=args$chrom_len_fpath,
                                wkd='/tmp')

  save(xing_diff_gr,file=ovr_apa_rd_fn,compress=T)
}

# -----------------------------------
jmessage('generate unique protein signal blocks ...')
exp_label <- sprintf('%s_ublocks',exp_label)
ublock_rd_fn <- file.path(wkd,sprintf('%s.rd',exp_label))

if (reuse & file.exists(ublock_rd_fn)) {
  load(ublock_rd_fn) #cvg_dt
} else {
  cvg_dts <- list()
  xing_diff.mval <- as.data.table(xing_diff_gr)
  xing_diff.pval <- as.data.table(xing_diff_gr)
  setnames(xing_diff.mval,'seqnames','chr')
  setnames(xing_diff.pval,'seqnames','chr')
  
  for (j in 1:length(apa_int_manorms)) {
    prot <- names(apa_int_manorms)[[j]]
    message(prot)
    ij <- findOverlaps(xing_diff_gr,apa_int_manorms[[j]])
    xing_diff.mval[,(prot):=0.]
    xing_diff.mval[queryHits(ij),(prot):=apa_int_manorms[[j]]$M_value[subjectHits(ij)]]
    
    xing_diff.pval[,(prot):=0.]
    xing_diff.pval[queryHits(ij),(prot):=apa_int_manorms[[j]]$P_value[subjectHits(ij)]]
  }
  
  xing_diff.mval$pos <- paste0(xing_diff.mval$chr,':',
                               xing_diff.mval$start,'-',
                               xing_diff.mval$end)
  
  xing_diff.pval$pos <- paste0(xing_diff.pval$chr,':',
                               xing_diff.pval$start,'-',
                               xing_diff.pval$end)
  
  save(xing_diff.mval,xing_diff.pval,file=ublock_rd_fn,compress=T)
}

jmessage(sprintf('<rec>total_uniq_cvg_region:%d',dim(xing_diff.mval)[1]))
jmessage(sprintf('<rec>total_uniq_cvg_region_bp:%d',sum(xing_diff.mval$width)))

pdf_file <- file.path(wkd,sprintf('%s_ublocks_width_hist.pdf',exp_label))
pdf(pdf_file)
hist(xing_diff.mval$width,breaks=100,
     main=sprintf('histogram of unique binding site width across \n ChIPSeq MANorm(M-value) + MBD binding scores [total:%d]',dim(xing_diff.mval)[1]))

dev.off()

summary(xing_diff.mval$width)

message(sprintf("retaining a minimum [%d]-bp genomic regions ...",args$bpco))
xing_diff.mval_filt <- xing_diff.mval[width>args$bpco,]

ublock_tsv_fn2 <- file.path(wkd,sprintf('%s_ublocks_filt%da.tsv',exp_label,args$bpco))
xing_diff.mval_filt2 <- xing_diff.mval_filt
xing_diff.mval_filt2$chr <- xing_diff.mval_filt2$pos
xing_diff.mval_filt2[,c('start','end','width','strand','ovr','pos'):=NULL]
setnames(xing_diff.mval_filt2,'chr','pos')

fwrite(xing_diff.mval_filt2,file=ublock_tsv_fn2,sep = '\t')

# -----------------------------
ublock_tsv_fn <- file.path(wkd,sprintf('%s_ublocks_filt%d.tsv',exp_label,args$bpco))
message(sprintf('annotating 3 regions of interest that we confirm that the ChIPSeq binding sites show a good evidence of cohesin complexes ...'))

goi3_gr <- GRanges(seqnames = c("chr6","chr7","chr7"),
                   strand = c("+", "+", "+"),
                   ranges = IRanges(start = c(41068475,807199,157204331), 
                                    end=c(41069346,808852,157204597)))
goi3_gr$gene <- c('NFYA','HEATR2','DNAJB6')
goi3_gr$idx <- 1:3

xing_diff.mval_filt.gr <- makeGRanges(xing_diff.mval_filt)
ij <- findOverlaps(xing_diff.mval_filt.gr, goi3_gr)
xing_diff.mval_filt.gr$gene <- ""
xing_diff.mval_filt.gr[queryHits(ij)]$gene <- goi3_gr[subjectHits(ij)]$gene
xing_diff.mval_filt <- as.data.table(xing_diff.mval_filt.gr)
# -----------------------------

xing_diff.mval_filt$strand <- NULL
xing_diff.mval_filt$ovr <- NULL
setnames(xing_diff.mval_filt,'seqnames','chr')
fwrite(xing_diff.mval_filt,file=ublock_tsv_fn,sep = '\t')

message(sprintf("use [%s] as an input to consensus clustering w/ NMF",ublock_tsv_fn))