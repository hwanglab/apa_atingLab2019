# Also run t-test and edgeR and compare

library(edgeR)
library(argparse)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

parser <- ArgumentParser(description='testbatchadj')

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 2,
										help="number of cpus to utilize [2]")

parser$add_argument("-i", "--input_file", type="character", required=TRUE,
										dest="input_file",
										help="rscript workspace rd file")

parser$add_argument("-w", "--prepropa_rd", type="character", required=TRUE,
										dest="prepropa_rd",
										help="rscript workspace rd (prpropa_rd) file")

parser$add_argument("-o", "--output_file", type="character", required=TRUE,
										dest="output_file",
										help="output file")
args <- parser$parse_args()

#load("./output/prepropa.rd")
load(args$prepropa)
#load("./output/apa.rd")
load(args$input_file)

# Want to do the test and drop results right into the existing object
# Keep adding on cols with the results of each next test
testTTest <- function(myta)
{
	dat <- myta$use.frac
	# Found cases where both groups have zero variance and t-test failed on them, so have to flag this and pass 0 p-value
	t.p <- mclapply(1:nrow(dat),function(x){
		#message(x);
		avec <- dat[x,][,1:4]
		bvec <- dat[x,][,5:8]
		if((sum(avec==max(avec))==4)&(sum(bvec==max(bvec))==4))
		{
			return(0)
		} else
		{
			td <- t.test(x=avec,y=bvec)
			return(td$p.value)
		}
	},mc.cores=args$ncpu)
	t.p <- do.call(c,t.p)
	t.padj <- p.adjust(t.p,method="fdr")
	myta$res$ttest_p <- t.p
	myta$res$ttest_padj <- t.padj
	return(myta)
}

apa.tt <- lapply(apa,testTTest)
outdir <- dirname(args$output_file)
if (!dir.exists(outdir)){dir.create(outdir, showWarnings=FALSE)}
edger_prefix <- file.path(outdir,'edger_')

testEdgeR <- function(myta)
{
	mycomp <- paste0(unique(myta$samp.sub$group)[1],"_vs_",unique(myta$samp.sub$group)[2])
	message("Testing: ",mycomp)

	res <- myta$res
	ann <- data.frame(tu=res$tu,pr=res$pr)
	g <- as.character(myta$samp.sub$group)
	adj <- factor(myta$samp.sub$batch)
	cnt <- myta$counts.raw
	y.all <- DGEList(counts=cnt,genes=ann,group=g)
	y <- calcNormFactors(y.all)
	design <- model.matrix(~ g + adj)
	y <- estimateDisp(y, design, robust=TRUE)
	fit <- glmQLFit(y, design, robust=TRUE)

	#pdf(file=paste0("output/edger_",mycomp,".pdf"))
	pdf(file=paste0(edger_prefix,mycomp,".pdf"))
	plotBCV(y)
	plotQLDisp(fit)
	plotMDS(y)
	dev.off()

	qlf <- glmQLFTest(fit, coef=2)
	sp <- diffSpliceDGE(fit, coef=2, geneid="tu", exonid="pr")

	myta$edger <- sp

	myta$res$edger_exon_p <- sp$exon.p.value
	myta$res$edger_exon_padj <- p.adjust(sp$exon.p.value,method="fdr")

	ftest <- sp$gene.p.value
	ftest.adj <- ftest
	ftest.adj[,1] <- p.adjust(ftest[,1],method="fdr")
	simes <- sp$gene.Simes.p.value

	myta$res$edger_ftest_p <- ftest[match(res$tu,rownames(ftest))]
	myta$res$edger_ftest_padj <- ftest.adj[match(res$tu,rownames(ftest.adj))]

	names(simes) <- res[match(names(simes),res$pr),]$tu
	simes.adj <- p.adjust(simes,method="fdr")

	myta$res$edger_simes_p <- simes[match(res$tu,names(simes))]
	myta$res$edger_simes_padj <- simes.adj[match(res$tu,names(simes.adj))]


	return(myta)
}

apa.er <- lapply(apa.tt,testEdgeR)

apa.alt <- apa.er
#save(apa.alt,file="output/apa.alt.rd",compress=T)
save(apa.alt,file=args$output_file,compress=T)

