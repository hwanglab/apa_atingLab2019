# Plot venns between the different p-value calls
library(UpSetR)
library(gridExtra)
library(VennDiagram)
library(argparse)
library(data.table)

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
load(args$prepropa_rd)

#to get total pA sites
pr_mean_raw <- as.data.table(pr$means$raw)
rownames(pr_mean_raw) <- rownames(pr$means$raw)

T1<-list()
ret <- colSums(pr_mean_raw[,c('HCT116','DKO')]>=10)
T1[['All.pA.sites.PR.HCT116']] <- ret[['HCT116']]
T1[['All.pA.sites.PR.DKO']] <- ret[['DKO']]

T1[['All.pA.sites.PR.union']] <- pr_mean_raw[HCT116>=10|DKO>=10,.N,]

gene_tu <- red_ens$tu
prs <- unique(rownames(pr_mean_raw)[pr_mean_raw$HCT116>=10])
tus <- unique(jtu$join[jtu$join$pr %in% prs]$tu)
genes <- unique(gene_tu$gene.id[gene_tu$tu %in% tus])
T1[['All.pA.sites.gene.HCT116']] <- length(genes)

prs <- unique(rownames(pr_mean_raw)[pr_mean_raw$DKO>=10])
tus <- unique(jtu$join[jtu$join$pr %in% prs]$tu)
genes <- unique(gene_tu$gene.id[gene_tu$tu %in% tus])
T1[['All.pA.sites.gene.DKO']] <- length(genes)

prs <- unique(rownames(pr_mean_raw)[pr_mean_raw$HCT116>=10|pr_mean_raw$DKO>=10])
tus <- unique(jtu$join[jtu$join$pr %in% prs]$tu)
ugenes <- unique(gene_tu$gene.id[gene_tu$tu %in% tus])
T1[['All.pA.sites.gene.union']] <- length(ugenes)

pru <- unique(jtu$join[(tu %in% tus) & unique_tu==T,pr])
T1[['Unambig.Genes.PR']] <- length(pru)

tuu <- unique(jtu$join[(tu %in% tus) & unique_tu==T,tu])

T1[['Unambig.Genes.gene']] <- length(unique(gene_tu$gene.id[gene_tu$tu %in% tuu]))

dt2 <- jtu$join[(tu %in% tuu) & unique_tu==T,]
dt2pr <- dt2[,.(N=length(pr)),by=c('tu')]
T1[['Genes.w.multiplePA.PR']] <- length(which(dt2pr$N>1))
T1[['Genes.w.multiplePA.gene']] <- dt2pr[N>1,sum(N)]

# ---------------------------
#load("./output/apa.alt.rd")
load(args$input_file)

# Just running for the final "intersect" set that we already decided on, using the batchadj runs of APA calling
# Not generating a sig call column for the alternative tests we aren't using

# --------------------------------------------------------------------
# Call a sig for each test and criteria
# Perhaps we do two rounds, one for each?

callsig <- function(myta)
{
	# P-value and fraction delta
	myta$res$batchadj_delta_sig <- ((abs(myta$res$b_minus_a_frac)>=0.1)&(myta$res$batchadj_padj<0.0001))

	# P-value and fold change
	myta$res$batchadj_perfc_sig <- with(myta$res,(batchadj_padj<0.0001)&(abs(l2_b_over_a_frac)>log2(1.5))&((mean_a_frac>=0.05)|(mean_b_frac>=0.05)))

	# P-value and fold change and delta fraction (intersect set we want to use)
	myta$res$int_sig <- myta$res$batchadj_delta_sig & myta$res$batchadj_perfc_sig 

	return(myta)
}
apa.sig <- lapply(apa.alt,callsig)
#save(apa.sig,file="output/apa.sig.rd",compress=T)
save(apa.sig,file=args$output_file,compress=T)
# --------------------------------------------------------------------

dt <- apa.sig$HCT116_vs_DKO$res
dt <- dt[batchadj_padj<0.0001,]

T1[['padj1e-4.PR']] <- dim(dt)[1]
T1[['padj1e-4.gene']] <-dt[,length(unique(gene_name))]

dt <- dt[((mean_a_frac>=0.05)|(mean_b_frac>=0.05)),]
T1[['frac5pct.PR']] <- dim(dt)[1]
T1[['frac5pct.gene']] <-dt[,length(unique(gene_name))]

dt <- dt[abs(l2_b_over_a_frac)>log2(1.5),]
T1[['b_over_a_frac_1.5.PR']] <- dim(dt)[1]
T1[['b_over_a_frac_1.5.gene']]  <-dt[,length(unique(gene_name))]

dt <- dt[abs(dt$b_minus_a_frac)>=0.1,]
T1[['b_minus_a_frac_0.1.PR']] <- dim(dt)[1]
T1[['b_minus_a_frac_0.1.gene']]  <-dt[,length(unique(gene_name))]
T1 <- t(as.data.table(t(T1)))
colnames(T1) <- 'count'
metrics <- rownames(T1)
T1 <- as.data.table(T1)
rownames(T1) <- metrics

tsv_fpath <- file.path(dirname(args$output_file),'table01_pA_APA_counts.tsv')
fwrite(T1,file=tsv_fpath,row.names = T,sep='\t')
