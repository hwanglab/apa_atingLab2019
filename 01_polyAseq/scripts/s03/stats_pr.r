library(ggplot2)
library(gridExtra)
library(argparse)
library(data.table)
library(stringr)
library(parallel)
library(handy)
library(openxlsx)

# Test for APA with DEXSeq
testApa2 <- function(samp,pr,jtu,a,out_prefix)
{
  
  message("processing PR ...")
  # stopifnot(a %in% samp$group)
  # stopifnot(b %in% samp$group)
  
  message("Starting with ",length(unique(jtu$join$pr))," pA sites accross ",length(unique(jtu$join$tu))," TUs")

  # Require mean >= 10 in at least a or b before moving to testing
  cov <- 1
  matmean <- pr$means$raw
  pr.expr <- rownames(matmean)[(matmean[,a]>=cov)]
  pr.a <- rownames(matmean)[(matmean[,a]>=cov)]
  
  
  # Grab unique TUs only (all PR of TU are unique)
  
  join.a <- jtu$join
  join.a <- join.a[(pr %in% pr.a)&(unique_tu==TRUE),]
  length(unique(join.a$tu))
  
  #browser()
  
  join.want <- jtu$join[((pr %in% pr.expr) & unique_tu==TRUE),]
  stopifnot(!str_detect(join.want$over_tus,","))
  stopifnot(!str_detect(join.want$flank_tus,","))
  stopifnot(all(join.want$unique_pr))
  stopifnot(length(unique(join.want$pr))==nrow(join.want))
  stopifnot(class(join.want$pr)=="character")
  stopifnot(join.want$pr %in% rownames(pr$counts$raw))
  join.want <- join.want[,list(pr=pr,type=type,coding=coding,unique_pr=unique_pr,unique_tu=unique_tu,over_tus=over_tus,flank_tus=flank_tus,num_pr=length(pr)),by="tu"]
  
  message("Subsetting to TUs with >1 pA site")
  message("After subsetting to TUs with all PRs unique: ",length(unique(join.want$pr))," sites accross ",length(unique(join.want$tu))," TUs")
  #join.want <- join.want[num_pr>1,] # ignore any tu where at most one pr appears only
  message("After subsetting to >1 pA site TUs: ",length(unique(join.want$pr))," sites accross ",length(unique(join.want$tu))," TUs")
  
  j <- match(join.want$tu,red_ens$tu$tu)
  join.want$gene_symbol <- red_ens$tu$name[j]

  num_pr_gene <- data.table(unique(join.want[,list(gene_symbol,num_pr)]))
  ggplot(num_pr_gene,aes(num_pr)) + 
    geom_histogram()

  fn <- sprintf('%s_num_pr_per_gene.tsv',out_prefix)
  fwrite(num_pr_gene[,c('gene_symbol','num_pr')],file=fn)
  # -------------------------
  
  join.want$mean_rcount_pr <- pr$means$raw[join.want$pr,4]
  
  j <- match(join.want$pr,pr$pr$pr)
  join.want <- cbind(join.want,as.data.table(pr$pr[j,])[,c('seqnames','start','end')])
  # --> mean_rcount_pr available
  ggplot(join.want,aes(mean_rcount_pr)) + 
    geom_histogram() +
    scale_x_log10()
  fn <- sprintf('%s_mean_rcount_pr.tsv',out_prefix)
  fwrite(join.want[,c('gene_symbol','pr','seqnames','start','end','mean_rcount_pr')],file=fn)
  # -------------------------
  
  join.want$pr_width <- join.want$end - join.want$start + 1
  # --> pr width is available
  
  ggplot(join.want,aes(pr_width)) + 
    geom_histogram() +
    scale_x_log10()
  fn <- sprintf('%s_pr_width.tsv',out_prefix)
  fwrite(join.want[,c('gene_symbol','pr','seqnames','start','end','pr_width')],file=fn)
  # -------------------------
  
  join.want.pr <- join.want[,c('pr','seqnames','start','end')]
  chrs <- handy::chrs()
  
  space_btn_prs <- NA
  for (chr in chrs){
    pr_chr <- join.want.pr[seqnames==chr,][order(start,end)]
    L <- dim(pr_chr)[1]
    space_btn_pr <- pr_chr$start[2:L] - pr_chr$end[1:L-1]
    space_btn_prs <- c(space_btn_prs,space_btn_pr)
  }
  
  btn_pr <- data.table(space=space_btn_prs)
  ggplot(btn_pr,aes(space)) + 
    geom_histogram() +
    scale_x_log10()
  
  x <- 1
}


pr_rd_file="/media/tommy/cache/exp_out/ating/polyASeqs/03_CallApa/output/prepropa.rd"
load(pr_rd_file)
a <- "HCT116"
testApa2(samp,pr,jtu,a,out_prefix="/media/tommy/cache/exp_out/ating/polyASeqs/03_CallApa/output/hct116_pr")
