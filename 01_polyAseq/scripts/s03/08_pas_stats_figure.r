library(argparse)
library(data.table)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

if (T) {
  parser <- ArgumentParser(description='prepropa')
  
  parser$add_argument("-i", "--prepropa", type="character", required=TRUE,
                      dest="prepropa",
                      help="prepropa rd file path")
  
  parser$add_argument("-d", "--polyadb_fpath", type="character", required=TRUE,
                      dest="polyadb_fpath",
                      help="polya database file path")
  
  parser$add_argument("-c", "--ctrl", type="character", required=FALSE,
                      dest="ctrl", default="HCT116",
                      help="control group tag [HCT116]")
  
  parser$add_argument("-e", "--expr", type="character", required=FALSE,
                      dest="expr", default="DKO",
                      help="experimental group tag [DKO]")
  
  parser$add_argument("-t", "--tmp_dir", type="character", required=FALSE,
                      dest="tmp_dir", default="/tmp",
                      help="temporary directory [/tmp]")
  
  parser$add_argument("-o", "--out_dir", type="character", required=TRUE,
                      dest="out_dir",
                      help="output directory")

  args <- parser$parse_args()
} else {
  args <- data.table(prepropa='../../01_wkd/out/03_CallApa/output/prepropa.rd',
                     ctrl='HCT116',
                     expr='DKO',
                     tmp_dir="/tmp",
                     polyadb_fpath="../../data/polya_db2/PolyADB_2.bed",
                     out_dir="../../01_wkd/out/03_CallApa/output")
}

load(args$prepropa)

library(goldmine)

# Want to show HCT/DKO univ vs DB vs Null
hexlist <- list()
# HCT/DKO pA sites

pr_count_cln <- colnames(pr$counts$raw)

for (tag in c(args$ctrl,args$expr)) {
  pr_clns <- pr_count_cln[startsWith(pr_count_cln,tag)]
  for (pr_cln in pr_clns) {
    pr_dt <- pr$counts$raw[,eval(pr_cln)]
    pr_dt <- pr_dt[pr_dt>0]
    pr_of_interest <- names(pr_dt)
    hexlist[[pr_cln]] <- countHexes(pr$pr[(pr$pr$pr %in% pr_of_interest),])
  }
}


#allpr <- pr$pr
#allpr <- countHexes(allpr)
#hexlist$allpr <- allpr

# PolyA DB 2
pad <- fread(args$polyadb_fpath)
pa <- with(pad,GRanges(V1,IRanges(V2,V3),strand=V6))
pa <- pa[seqnames(pa) %in% handy::chrs()]
pa <- unique(pa)
padp <- countHexes(gr=pa)
hexlist$pAdb2 <- padp

# Null
tmp_dir <- file.path(args$tmp,"10_pas_stats_figure")
if (!dir.exists(tmp_dir)) {dir.create(tmp_dir)}

cnt <- pr$means$raw
hctdko <- pr$pr[(cnt[,args$ctrl]>=10)|(cnt[,args$expr]>=10)]
hctdko <- countHexes(hctdko)

dgp <- drawGenomePool(query=hctdko,
                      n=10,
                      genome="hg19",
                      cachedir=tmp_dir,
                      sync=FALSE)

null <- countHexes(gr=dgp)
hexlist$background <- null 

# Make tables
hextabs <- lapply(hexlist,function(x) as.data.frame(table(x$hex)))
for(h in 1:length(hextabs))
{
  hextabs[[h]]$run <- names(hexlist)[h]   
}
tab <- rbindlist(hextabs)

# Do plots
cols <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#969696","#666666")

hexlist_names <- names(hexlist)
L <- length(hexlist_names)
sample_names <- hexlist_names[1:(L-2)]

tab.per <- tab[,list(Var1=Var1,Freq=Freq/sum(Freq)),by="run"]
tab.per[run=="background",run:="null"]
tab.per$run <- factor(tab.per$run,levels=c(sample_names,"pAdb2","null"))
tab.per$motif <- tab.per$Var1

library(ggplot2)

if (!dir.exists(args$out_dir)) {dir.create(args$out_dir)}
out_fpath <- file.path(args$out_dir,'hexamer.pdf')

pdf(file=out_fpath,width=4,height=4)
#ggplot(tab,aes(x=run,y=Freq,fill=Var1)) + geom_bar(stat="identity") + handy::ggnice() + scale_fill_manual(values=cols)
p <- ggplot(tab.per,aes(x=run,y=Freq,fill=motif)) + 
  geom_bar(stat="identity") + 
  handy::ggnice() + 
  scale_fill_manual(values=cols) + 
  labs(y="Fraction of Regions",x="") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

print(p)

dev.off()

if (dir.exists(tmp_dir)) {unlink(tmp_dir,recursive = T)}
