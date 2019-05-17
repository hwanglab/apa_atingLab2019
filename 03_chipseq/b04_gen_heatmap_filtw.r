library(argparse)
library(data.table)

if (F) {
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  
  parser <- ArgumentParser(description='prepropa')
  
  parser$add_argument("-t", "--tag_file", type="character", required=TRUE,
                      dest="tag_file",
                      help="original NMF input TSV file")
  
  parser$add_argument("-c", "--cluster_file", type="character", required=TRUE,
                      dest="cluster_file",
                      help="NMF output TSV file")
  
  parser$add_argument("-c", "--col2analy", type="character", required=TRUE,
                      dest="col2analy",
                      help="column to analyze in the input_tsv file")
  
  parser$add_argument("-t", "--title", type="character", required=FALSE,
                      dest="title", default="nmf9_wo_filtering",
                      help="title")
  
  parser$add_argument("-r", "--reuse", type="integer", required=FALSE,
                      dest="reuse",default=1,
                      help="Yes:[1], No:0")
  
  args <- parser$parse_args()
} else {
  args <- data.table(tag_file='fastq/cluster_analy/chipseq_xing_ublocks_filt10.tsv',
                     cluster_file='fastq/cluster_analy/figures/chipseq_xing_ublocks_filt10a_NMF_groupInfo.txt',
                     col2analy='(K=9)',
                     title="NMF.k=9, log2(rawMBD+1), region_length>=10bp",
                     r=0)
}

library(openxlsx)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
library(goldmine)
source(file.path(Sys.getenv('R_UTIL'),'lib_chipseq.r'))

dt0 <- fread(args$cluster_file)
cluster_cidx <- which(colnames(dt0)==args$col2analy)

class_dt <- dt0[,list(ID,get(args$col2analy))]
setnames(class_dt,c('ID','V2'),c('pos','cluster_id'))

tag_dt<-fread(args$tag_file)

message('applying log2() to MBD mean depth ...')
tag_dt$HCT116_Mbd <- log2(tag_dt$HCT116_Mbd + 1.)
tag_dt$DKO_Mbd <- log2(tag_dt$DKO_Mbd + 1.)

tag_dt$cluster_id <- class_dt[match(tag_dt$pos,class_dt$pos),'cluster_id']
tag_dt <- tag_dt[order(chr,start,end)]

clustered.gr <- makeGRanges(tag_dt)

apa_call <- 'ShorterInB'

message("annotated genes of interest ...")
apa_colon_xlsx <- "../01_polyAseq/01_wkd/out/04_AnnoApa/output/apa.ann_HCT116_vs_DKO.xlsx"
stopifnot(file.exists(apa_colon_xlsx))
apa <- read_supp_tab3(apa_colon_xlsx,sheet_i=1) #keep apa to bring annotation back
apa.gr <- makeGRanges(apa[call==apa_call,])

ij <- findOverlaps(clustered.gr,apa.gr)
clustered.gr$gene[queryHits(ij)] <- apa.gr$gene_name[subjectHits(ij)]

tag_dt <- as.data.table(clustered.gr)
tag_dt$strand <- NULL
setnames(tag_dt,'seqnames','chr')

library(ggplot2)
library(gplots)
library(cluster)    # clustering algorithms
library(dendextend) # for comparing two dendrograms
library(factoextra)

cluster_dist_file <- './fastq/cluster_analy/chipseq_xing_ublocks_filt10a_NMF_DistanceMat_K9.txt'
exp_label <- 'chipseq'
wkd <- 'fastq/cluster_analy'

cmd <- sprintf("sed 's/(//g' %s | sed 's/)//g' | cut -f1 --complement",cluster_dist_file)

x <- fread(cmd = cmd)

hclust_dendro <- hclust(as.dist(x))

n <- length(unique(tag_dt$cluster_id))
ccls <- data.table(cluster_id=rev(hclust_dendro$order),cls_order=1:n)

ccls <- ccls[order(cluster_id)]
k9_pal <- c('#191970','#006400','#ff4500','#ffd700','#00ff00','#00bfff','#0000ff','#ff69b4','#ffe4c4')
ccls$col <- k9_pal

tag_dt$cls_order <- ccls[match(tag_dt$cluster_id,ccls$cluster_id),cls_order]

md_dt_cls <- tag_dt[order(cls_order,-width,gene)]

mycol <- colorpanel(20,"orange", "white", "blue")

regChip <- as.matrix(md_dt_cls[,5:12,with=F]) #chipseq manorm M-value
message('scaling ChIP-seq read depth ...')
chIPSeqs <- apply(regChip[,1:6],2,function(regChipj){(regChipj - mean(regChipj))/sd(regChipj)})
regChip[,1:6] <- chIPSeqs

message('scaling MBD read depth across two samples ...')
mbdSeqs <- (regChip[,7:8] - mean(regChip[,7:8]))/sd(regChip[,7:8])
regChip[,7:8] <- mbdSeqs

message(sprintf('dim of matrix = %d x %d',dim(regChip)[1],dim(regChip)[2]))

rd_fn <- file.path(wkd,sprintf('%s_cluster.rd',exp_label))
save(md_dt_cls,regChip,hclust_dendro,file=rd_fn,compress=T)

wb <- createWorkbook()

sheet_name <- 'input_chipseq'
addWorksheet(wb,sheet_name)
writeData(wb, sheet=sheet_name,x=md_dt_cls)

sheet_name <- 'scaled_chipseq'
addWorksheet(wb,sheet_name)
writeData(wb, sheet=sheet_name,x=regChip)

saveWorkbook(wb, file.path(wkd,sprintf('%s_cluster.xlsx',exp_label)), overwrite = TRUE)

pdf_fn <- file.path(wkd,sprintf('%s_heatmap.pdf',exp_label))
pdf(pdf_fn)

plot(hclust_dendro)

heatmap.2(regChip,
          Rowv=F,
          Colv=T,
          col=mycol,
          RowSideColors=ccls[match(md_dt_cls$cluster_id,ccls$cluster_id),col],
          scale="none",
          density.info="none",
          dendrogram = "column",
          trace="none",
          srtCol=25,
          main=args$title)

dev.off()
x<-1