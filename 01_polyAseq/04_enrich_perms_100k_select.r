library(goldmine)

source(file.path(Sys.getenv('R_UTIL_APA'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL_APA'),'lib_apps.R'))
mount_prefix <- get_mount_dir()

apa_sig_rd <- file.path(mount_prefix,'apa_atingLab2019','01_polyAseq','01_wkd','out','03_CallApa','output','apa.sig.rd')
load(apa_sig_rd)

summary(apa.sig)

getIntSet <- function(i)
{
	message(names(apa.sig)[i])
	myres <- apa.sig[[i]]$res
	myres$sig_intersect <- with(myres,(batchadj_perfc_sig==TRUE)&(batchadj_delta_sig==TRUE))
	sigtu <- unique(myres[sig_intersect==TRUE,]$tu)
	myres$tu_sig <- myres$tu %in% sigtu
	myres <- myres[myres$tu_sig,]
	spl <- split(myres,myres$tu)
	mytus <- rev(unique(myres[order(batchadj_padj,decreasing=T),]$tu))
	spl2 <- spl[mytus]
	out <- rbindlist(spl2)
}
apa.int <- lapply(1:length(apa.sig),getIntSet)


# Compare the DeltaWtMean with MostUpMostDown
callWtMean <- function(apa)
{
	calls <- apa[,list(gene_name=gene_name[1],chr=chr[1],strand=strand[1],wa=weighted.mean(x=end,w=mean_a_frac),wb=weighted.mean(x=end,w=mean_b_frac)),by="tu"]
	calls$wb_minus_wa <- calls$wb-calls$wa

	# Call direction based on strand and value
	cp <- calls[strand=="+",]
	cp$wcall <- ""
	cp[wb>wa,]$wcall <- "LongerInB"
	cp[wb<wa,]$wcall <- "ShorterInB"
	stopifnot(cp$wcall!="")
	cm <- calls[strand=="-",]
	cm$wcall <- ""
	cm[wb>wa,]$wcall <- "ShorterInB"
	cm[wb<wa,]$wcall <- "LongerInB"
	stopifnot(cm$wcall!="")
	calls <- rbind(cp,cm)

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
	upp[up_start<down_start,]$mcall <- "ShorterInB"
	upp[up_start>down_start,]$mcall <- "LongerInB"
	stopifnot(upp$mcall!="")
	
	upm <- ud[strand=="-",]
	upm$mcall <- ""
	upm[up_start<down_start,]$mcall <- "LongerInB"
	upm[up_start>down_start,]$mcall <- "ShorterInB"
	stopifnot(upm$mcall!="")

	udmeth <- rbind(upp,upm)

	calls <- cbind(calls[match(udmeth$tu,calls$tu),],udmeth[,c("MostUp","MostDown","up_start","up_end","down_start","down_end","mcall"),with=F])
	calls$agree <- calls$wcall==calls$mcall
	calls
}
apa.calls <- lapply(apa.int,callWtMean)

if (!dir.exists('./04_wkd')) {
  dir.create('./04_wkd')
}

#for(i in 1:length(names(apa.sig)))
#{
#	comp <- names(apa.sig)[i]
#	message(comp)
#	apa.calls[[i]]$comp <- comp
#}
#agg <- rbindlist(apa.calls)

#am <- agg[,list(n=length(tu)),by=c("comp","mcall")]
#setnames(am,"mcall","call")
#am$method <- "MostUpMostDown"
#aw <- agg[,list(n=length(tu)),by=c("comp","wcall")]
#setnames(aw,"wcall","call")
#aw$method <- "DeltaWtAvg"
#atab <- rbind(am,aw)

#atab <- atab[order(atab$call),]

#atab$comp <- factor(atab$comp,levels=names(apa.sig))

#pdf(file="04_wkd/compare_methods.pdf",width=10.5,height=4)
#ggplot(atab,aes(x=method,y=n,fill=call)) + geom_bar(stat="identity") + facet_wrap(~comp,scales="free_y",nrow=1) + handy::ggnice() + scale_fill_manual(values=rev(c("#4daf4a","#984ea3"))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#dev.off()

#disses <- agg[(agree==FALSE)&(comp=="HCT116_vs_DKO"),]
#write.csv(disses,row.names=F,file="04_wkd/HCT_DKO_method_disagree.csv")

# --------------------------------------------------------------------
# All ENCODE enrichment of changing PR flanks
genome <- "hg19"
cachedir <- "~/temp"

feat <- getFeatures("wgEncodeRegTfbsClusteredV3",genome=genome,cachedir=cachedir)[[1]]
enames <- feat$name
values(feat) <- NULL
feat$name <- enames

# Also append in the CTCF from HCT/DKO ChIP-seq
dr <- data.table(read.table("data/tf_binding_precomputed/diffreps/HCT_vs_DKO_report_win200_p0.0001.txt",comment.char="#",header=TRUE))
dr <- dr[abs(log2FC)>log2(1.5),]
dr.gr <- with(dr,GRanges(Chrom,IRanges(Start,End)))
dr.gr$name <- ""
dr.gr[dr$Event=="Down",]$name <- "CTCF_LostDko"
dr.gr[dr$Event=="Up",]$name <- "CTCF_GainedDko"
dr.gr <- dr.gr[seqnames(dr.gr) %in% handy::chrs()]

mac <- dir("data/tf_binding_precomputed/macs",pattern="_peaks.bed$",full.names=TRUE)
mac <- mac[str_detect(mac,"HCT|DKO")]
peaks <- lapply(mac,fread)
peak.gr <- lapply(peaks,function(x) with(x,GRanges(V1,IRanges(V2+1,V3))))
con <- reduce(do.call(c,peak.gr))
con <- con[seqnames(con) %in% handy::chrs()]
ovr <- lapply(peak.gr,function(x) con %over% x)
mat <- do.call(cbind,ovr)

any <- con
any <- any[!(any %over% dr.gr)]
any$name <- "CTCF_StaticAny"

all <- con[rowSums(mat)==ncol(mat)]
all <- all[!(all %over% dr.gr)]
all$name <- "CTCF_StaticAll"

feat <- c(feat,dr.gr,any,all)
table(feat$name)

etf <- split(feat,feat$name)
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Trying three-part region idea

mycalls <- apa.calls[[1]]
mycalls <- mycalls[mcall=="ShorterInB",]

btwp <- mycalls[strand=="+",list(tu=tu,strand=strand,chr=chr,start=up_start,end=down_end)]
btwm <- mycalls[strand=="-",list(tu=tu,strand=strand,chr=chr,start=down_start,end=up_end)]
btw <- rbind(btwp,btwm)
btw.gr <- makeGRanges(btw)

btw$size <- width(btw.gr)
stopifnot(btw$tu==mycalls$tu)
mycalls$btwsize <- btw$size

dsterm <- rbind(mycalls[strand=="+",list(tu=tu,
                                         strand=strand,
                                         chr=chr,
                                         start=down_end+1,
                                         end=down_end+btwsize)],
                mycalls[strand=="-",list(tu=tu,
                                         strand=strand,
                                         chr=chr,
                                         start=down_start-btwsize,
                                         end=down_start-1)])

usint <- rbind(mycalls[strand=="+",list(tu=tu,
                                        strand=strand,
                                        chr=chr,
                                        start=up_start-btwsize,
                                        end=up_start-1)],
               mycalls[strand=="-",list(tu=tu,
                                        strand=strand,
                                        chr=chr,
                                        start=up_start+1,
                                        end=up_start+btwsize)])

ds.gr <- makeGRanges(dsterm)
us.gr <- makeGRanges(usint)

summary(width(btw.gr))
summary(width(ds.gr))
summary(width(us.gr))

btw$name <- mycalls$gene_name
set.seed(20160124)
dgp.gr <- drawGenomePool(query=btw.gr,n=15,genome=genome,cachedir=cachedir)
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# PERMUTATION TESTING

# Want to know - how many genes can possibly produce 2 PRs drawn with distance <= max of real APA AND >= min of real PA
all <- apa.sig[[1]]$res
length(unique(all$tu))
all <- makeGRanges(all)
spl <- split(all,all$tu)
#d <- mcols(distanceToNearest(spl[[2]]))$distance
d <- mclapply(1:length(spl),function(x){
	message(x)
	mcols(distanceToNearest(spl[[x]]))$distance
},mc.cores=20,mc.preschedule=TRUE)

mymin <- min(width(btw.gr))
mymax <- max(width(btw.gr))

works <- sapply(d,function(x) any((x>=mymin)&(x<=mymax)))
table(works)
pool <- names(spl)[works]

# Now drawing the perms
set.seed(20160125)
ntu <- length(btw.gr)

# It makes sense to pre-compute which PR pairs are valid within each gene, so we don't need overhead for this each perm
all <- apa.sig[[1]]$res
length(unique(all$tu))
all <- makeGRanges(all,strand=TRUE)
spl <- split(all,all$tu)
spl <- spl[works]
d2 <- mclapply(1:length(spl),function(x){
	message(x)
	df <- data.table(as.data.frame(distanceToNearest(spl[[x]])))
	data.table(tu=spl[[x]]$tu[1],pr1=spl[[x]]$pr[df$queryHits],pr2=spl[[x]]$pr[df$subjectHits],distance=df$distance,pass=(df$distance>=mymin)&(df$distance<=mymax),chr=as.vector(seqnames(spl[[x]])[1]),up_start=as.vector(start(spl[[x]]))[df$subjectHits],up_end=as.vector(end(spl[[x]]))[df$subjectHits],down_start=as.vector(start(spl[[x]]))[df$queryHits],down_end=as.vector(end(spl[[x]]))[df$queryHits],strand=as.vector(strand(spl[[x]]))[1])
	
},mc.cores=20,mc.preschedule=TRUE)
pr_pool <- rbindlist(d2)
pr_pool <- pr_pool[pass==TRUE,]

doperm <- function()
{
	permtu <- sample(pool,ntu)
	permpr <- pr_pool[tu %in% permtu,]
	pspl <- split(permpr,permpr$tu)
	pbtw <- rbindlist(lapply(pspl,function(x) x[sample(1:nrow(x),1),]))
	
	pbtw.gr <- with(pbtw[up_start<down_end,],GRanges(chr,IRanges(up_start,down_end),mystrand=strand))
	pbtw.gr <- c(pbtw.gr,with(pbtw[up_start>down_end,],GRanges(chr,IRanges(down_start,up_end),mystrand=strand)))
	pbtw.gr$btwsize <- width(pbtw.gr)

	null_btw <- makeDT(pbtw.gr)
	flank1 <- GRanges(null_btw$chr,IRanges(null_btw$start-null_btw$btwsize,null_btw$end-1))
	flank2 <- GRanges(null_btw$chr,IRanges(null_btw$start+1,null_btw$end+null_btw$btwsize))

	null_ds.gr <- c(flank1[null_btw$mystrand=="-"],flank2[null_btw$mystrand=="+"])
	null_us.gr <- c(flank1[null_btw$mystrand=="+"],flank2[null_btw$mystrand=="-"])

	list(btw=pbtw.gr,us=null_us.gr,ds=null_ds.gr)
}
perms <- mclapply(1:100000,function(x){message(x);doperm();},mc.cores=20)

# Now that we have the perms, we need to compute enrichment accross them and the real data

# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Unfied enrichment functions for real and perms and application over all the perms
# x is GRanges of the factor of interest
# regions is a list of the 3-part regions
doEnrich <- function(x,regions)
{
	# Total count method
	#btw_n <- sum(regions$btw %over% x)		
	#ds_n <- sum(regions$ds %over% x)		
	#us_n <- sum(regions$us %over% x)		
	#btw_total <- length(regions$btw)	
	#ds_total <- length(regions$ds)
	#us_total <- length(regions$us)
	#btw_frac <- btw_n/btw_total
	#ds_frac <- ds_n/ds_total	
	#us_frac <- us_n/us_total	

	# Per-bp rate method
	btw_rate <- sum(width(intersect(regions$btw,x)))/sum(width(reduce(regions$btw)))
	ds_rate <- sum(width(intersect(regions$ds,x)))/sum(width(reduce(regions$ds)))
	us_rate <- sum(width(intersect(regions$us,x)))/sum(width(reduce(regions$us)))

	#data.table(factor=x$name[1],ds_n,btw_n,us_n,ds_total,btw_total,us_total,ds_frac,btw_frac,us_frac,ds_rate,btw_rate,us_rate)
	data.table(factor=x$name[1],us_rate,btw_rate,ds_rate)
}

reals <- list(btw=btw.gr,us=us.gr,ds=ds.gr)

realres <- rbindlist(mclapply(etf,function(x) doEnrich(x=x,regions=reals),mc.cores=20))

permres <- lapply(1:5,function(y){message(y);rbindlist(mclapply(etf,function(x) doEnrich(x=x,regions=perms[[y]]),mc.cores=20));})

# Could try just CTCF to go faster
fax <- etf[(str_detect(names(etf),"CTCF"))|(names(etf) %in% c("MYC","MBD4","SUZ12","HDAC2","GATA2","SMC3","EP300","MYBL2","YY1","GABPA","USF1"))]
realres <- rbindlist(mclapply(fax,function(x) doEnrich(x=x,regions=reals),mc.cores=20))
realmat <- as.matrix(realres[,-1,with=F])
rownames(realmat) <- realres$factor

permres <- mclapply(1:length(perms),function(y){message(y);rbindlist(lapply(fax,function(x) doEnrich(x=x,regions=perms[[y]])));},mc.cores=20)

save(perms,permres,compress=T,file="04_wkd/perm_100K_ctcf.rd")

resmat <- lapply(permres,function(x) as.matrix(x[,-1,with=F]))

permmeans <- apply(simplify2array(resmat),c(1,2),mean)
permsds <- apply(simplify2array(resmat),c(1,2),sd)
fc <- realmat/permmeans
l2fc <- log2(realmat/permmeans)
zscores <- (realmat-permmeans)/permsds


colnames(realmat) <- paste("real",colnames(realmat),sep="_")
colnames(permmeans) <- paste("permmean",colnames(permmeans),sep="_")
colnames(l2fc) <- paste("l2fc",colnames(l2fc),sep="_")
colnames(zscores) <- paste("z",colnames(zscores),sep="_")

out <- as.data.table(cbind(data.table(factor=rownames(realmat)),realmat,permmeans,l2fc,zscores))
out <- out[order(out$z_btw_rate,decreasing=T),]
write.csv(out,file="04_wkd/perm_100K_ctcf.csv",row.names=F)


dat <- rbindlist(permres)
dat <- melt(dat,id.vars=c("factor"))
pdf(file="04_wkd/Enrich_100K_ctcf.pdf",width=6,height=5)
for(i in out$factor)
{
	message(i)
	vdat <- data.frame(variable = c("us_rate", "btw_rate", "ds_rate"), xp = c(realres[factor==i,]$us_rate,realres[factor==i,]$btw_rate,realres[factor==i,]$ds_rate))
	print(ggplot(dat[factor==i,],aes(x=value,color=variable)) + geom_line(stat="density") + handy::ggnice() + labs(title=i) + geom_vline(data=vdat,aes(xintercept=xp,color=variable)) + scale_color_manual(values=c(us_rate="#4daf4a",btw_rate="#e41a1c",ds_rate="#377eb8")) + facet_grid(variable~.))
}
dev.off()


# --------------------------------------------------------------------


# Also want NoAPA genes, with the same scheme
sigtu <- unique(apa.int[[1]]$tu)
nstu <- unique(apa.sig[[1]]$res$tu)
nstu <- nstu[!(nstu %in% sigtu)]
length(nstu)
ns <- apa.sig[[1]]$res[tu %in% nstu,]
spl <- split(ns,ns$tu)
set.seed(20160124)
spl2 <- lapply(spl,function(x) x[sample(1:nrow(x),2),])
spl3 <- lapply(spl2,function(x) x[,list(tu=tu[1],gene_name=gene_name[1],chr=chr[1],strand=strand[1],MostUp=pr[1],MostDown=pr[2],up_start=start[1],up_end=end[1],down_start=start[2],down_end=end[2])])
nullcall <- rbindlist(spl3)

null_btw.gr <- with(nullcall[up_start<down_end,],GRanges(chr,IRanges(up_start,down_end),mystrand=strand))
null_btw.gr <- c(null_btw.gr,with(nullcall[up_start>down_end,],GRanges(chr,IRanges(down_start,up_end),mystrand=strand)))
null_btw.gr$btwsize <- width(null_btw.gr)
summary(null_btw.gr$btwsize)

null_btw <- makeDT(null_btw.gr)
flank1 <- GRanges(null_btw$chr,IRanges(null_btw$start-null_btw$btwsize,null_btw$end-1))
flank2 <- GRanges(null_btw$chr,IRanges(null_btw$start+1,null_btw$end+null_btw$btwsize))

null_ds.gr <- c(flank1[null_btw$mystrand=="-"],flank2[null_btw$mystrand=="+"])
null_us.gr <- c(flank1[null_btw$mystrand=="+"],flank2[null_btw$mystrand=="-"])

null_dgp.gr <- drawGenomePool(query=null_btw.gr,n=1,genome=genome,cachedir=cachedir)

doEnrich <- function(x)
{
	# Total count method
	btw_n <- sum(btw.gr %over% x)		
	ds_n <- sum(ds.gr %over% x)		
	us_n <- sum(us.gr %over% x)		
	null_n <- sum(dgp.gr %over% x)
	btw_total <- length(btw.gr)	
	ds_total <- length(ds.gr)
	us_total <- length(us.gr)
	null_total <- length(dgp.gr)	
	btw_frac <- btw_n/btw_total
	ds_frac <- ds_n/ds_total	
	us_frac <- us_n/us_total	
	null_frac <- null_n/null_total

	noapa_btw_n <- sum(null_btw.gr %over% x)		
	noapa_ds_n <- sum(null_ds.gr %over% x)		
	noapa_us_n <- sum(null_us.gr %over% x)		
	noapa_null_n <- sum(null_dgp.gr %over% x)
	noapa_btw_total <- length(null_btw.gr)	
	noapa_ds_total <- length(null_ds.gr)
	noapa_us_total <- length(null_us.gr)
	noapa_null_total <- length(null_dgp.gr)	
	noapa_btw_frac <- noapa_btw_n/noapa_btw_total
	noapa_ds_frac <- noapa_ds_n/noapa_ds_total	
	noapa_us_frac <- noapa_us_n/noapa_us_total	
	noapa_null_frac <- noapa_null_n/noapa_null_total

	# Per-bp rate method
	btw_rate <- sum(width(intersect(btw.gr,x)))/sum(width(reduce(btw.gr)))
	ds_rate <- sum(width(intersect(ds.gr,x)))/sum(width(reduce(ds.gr)))
	us_rate <- sum(width(intersect(us.gr,x)))/sum(width(reduce(us.gr)))
	null_rate <- sum(width(intersect(dgp.gr,x)))/sum(width(reduce(dgp.gr)))

	noapa_btw_rate <- sum(width(intersect(null_btw.gr,x)))/sum(width(reduce(null_btw.gr)))
	noapa_ds_rate <- sum(width(intersect(null_ds.gr,x)))/sum(width(reduce(null_ds.gr)))
	noapa_us_rate <- sum(width(intersect(null_us.gr,x)))/sum(width(reduce(null_us.gr)))
	noapa_null_rate <- sum(width(intersect(null_dgp.gr,x)))/sum(width(reduce(null_dgp.gr)))

	data.table(factor=x$name[1],ds_n,btw_n,us_n,null_n,ds_total,btw_total,us_total,null_total,ds_frac,btw_frac,us_frac,null_frac,ds_rate,btw_rate,us_rate,null_rate,noapa_ds_n,noapa_btw_n,noapa_us_n,noapa_null_n,noapa_ds_total,noapa_btw_total,noapa_us_total,noapa_null_total,noapa_ds_frac,noapa_btw_frac,noapa_us_frac,noapa_null_frac,noapa_ds_rate,noapa_btw_rate,noapa_us_rate,noapa_null_rate)
}
out <- mclapply(etf,doEnrich,mc.cores=20)

res <- rbindlist(out)

res$frac_diff <- res$btw_frac-res$null_frac
res <- res[order(res$frac_diff,decreasing=T),]
write.csv(res,file="04_wkd/enrich_3part_dyn_noapa.csv",row.names=F)

# Real APA
mat <- res[,list(factor=factor,us_fc=us_frac/null_frac,btw_fc=btw_frac/null_frac,ds_fc=ds_frac/null_frac,us_fc_bp=us_rate/null_rate,btw_fc_bp=btw_rate/null_rate,ds_fc_bp=ds_rate/null_rate)]

mat <- as.matrix(mat[,-1,with=F])
rownames(mat) <- res$factor

matl2 <- log2(mat)

# No APA
mat_noapa <- res[,list(factor=factor,noapa_us_fc=noapa_us_frac/noapa_null_frac,noapa_btw_fc=noapa_btw_frac/noapa_null_frac,noapa_ds_fc=noapa_ds_frac/noapa_null_frac,noapa_us_fc_bp=noapa_us_rate/noapa_null_rate,noapa_btw_fc_bp=noapa_btw_rate/noapa_null_rate,noapa_ds_fc_bp=noapa_ds_rate/noapa_null_rate)]

mat_noapa <- as.matrix(mat_noapa[,-1,with=F])
rownames(mat_noapa) <- res$factor

matl2_noapa <- log2(mat_noapa)

mat_diff <- matl2-matl2_noapa

mat_out <- cbind(mat,mat_noapa)
mat_out_l2 <- cbind(matl2,matl2_noapa,mat_diff)

write.csv(mat_out_l2,row.names=T,file="04_wkd/matl2_dyn_noapa.csv")
write.csv(mat_out,row.names=T,file="04_wkd/matfc_dyn_noapa.csv")

mat_vs <- res[,list(factor=factor,us_fc=us_frac/noapa_us_frac,btw_fc=btw_frac/noapa_btw_frac,ds_fc=ds_frac/noapa_ds_frac,us_fc_bp=us_rate/noapa_us_rate,btw_fc_bp=btw_rate/noapa_btw_rate,ds_fc_bp=ds_rate/noapa_ds_rate)]
mat_vs <- as.matrix(mat_vs[,-1,with=F])
rownames(mat_vs) <- res$factor

mat_vs_l2 <- log2(mat_vs)

write.csv(mat_vs_l2,row.names=T,file="04_wkd/matl2_dyn_noapa_VS.csv")
write.csv(mat_vs,row.names=T,file="04_wkd/matfc_dyn_noapa_VS.csv")

# --------------------------------------------------------------------
