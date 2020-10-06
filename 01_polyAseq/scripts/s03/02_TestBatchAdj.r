library(handy)
library(argparse)
library(data.table)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

if (TRUE) {

  parser <- ArgumentParser(description='testbatchadj')
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
  										dest="ncpu", default = 2,
  										help="number of cpus to utilize [2]")
  
  
  parser$add_argument("-d", "--prs_gt", type="integer", required=FALSE,
  										dest="prs_gt", default = 1,
  										help="greater than PRs to consider [1]")
  
  
  parser$add_argument("-i", "--input_file", type="character", required=TRUE,
  										dest="input_file",
  										help="rscript workspace rd file")
  
  parser$add_argument("-c", "--comp_file", type="character", required=TRUE,
  										dest="comp_file",
  										help="a file containing a pair of group samples to compare with")
  
  parser$add_argument("-o", "--output_file", type="character", required=TRUE,
  										dest="output_file",
  										help="output file")
  args <- parser$parse_args()
}

load(args$input_file)

# Compare with this adjuster and without it
adjust.var <- "batch"

# Do the tests and call the sigs w/ and w/o batch adj.
# take comp_group_fn as a user input

compdt <- fread(args$comp_file,header=TRUE)

#comps <- data.table(a=c("HCT116","HCT116","PrEC","PrEC","LNCaP"),b=c("DKO","DICER","LNCaP","DU145","DU145"))

comps <- data.table(a=compdt$group1,b=compdt$group2)

dotests <- function(a,b,ncpu=1,prs_gt=1)
{
	message("Testing: ",a," vs ",b)
	mycomp <- paste0(a,"_vs_",b)
	ta <- testApa(samp=samp,pr=pr,jtu=jtu,a=a,b=b,adjust.var=NULL,ncpu=ncpu,prs_gt=prs_gt)
	taadj <- testApa(samp=samp,pr=pr,jtu=jtu,a=a,b=b,adjust.var=adjust.var,prs_gt=prs_gt)
	
	list(ta=ta,taadj=taadj)

}
compnames <- paste0(comps$a,"_vs_",comps$b)
apa <- mclapply(1:nrow(comps),function(x) dotests(a=comps[x,]$a,b=comps[x,]$b,prs_gt=args$prs_gt),mc.cores=args$ncpu)
names(apa) <- compnames

# Merge down the results to one table going forward
handy::nsummary(apa)

# Put the adj info into the main res so get one ta object for each

tas <- lapply(apa,function(x){
	x$taadj$res$noadj_p <- x$ta$res$p
	x$taadj$res$noadj_padj <- x$ta$res$padj
	#x$taadj$res$noadj_sig <- x$ta$res$sig

	x$taadj$res$batchadj_p <- x$taadj$res$p
	x$taadj$res$batchadj_padj <- x$taadj$res$padj
	#x$ta$res$batchadj_sig <- x$taadj$res$sig

	x$taadj$res[,p:=NULL]
	x$taadj$res[,padj:=NULL]

	return(x$taadj)
})

apa <- tas

#save(apa,compress=T,file="output/apa.rd")
save(apa,compress=T,file=args$output_file)
