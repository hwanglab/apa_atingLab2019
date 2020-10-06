# Save all the parts we'd like to have for the purposes of screening the APA calls
library(reshape)
library(matrixStats)
library(ggbio)
library(goldmine)
library(handy)
library(argparse)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

parser <- ArgumentParser(description='geneplots')

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 1,
										help="number of cpus to utilize [1]")

parser$add_argument("-i", "--input_file", type="character", required=TRUE,
										dest="input_file",
										help="rscript workspace rd file")

parser$add_argument("-w", "--prepropa_rd", type="character", required=TRUE,
										dest="prepropa_rd",
										help="rscript workspace rd (prpropa_rd) file")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output dir")
args <- parser$parse_args()

#load("./output/prepropa.rd")
load(args$input_file)
#load("./output/apa.sig.rd")
load(args$prepropa_rd)

# Make plots
plotOne <- function(mytu,myta)
{
	mydat <- myta$res[tu==mytu,]
	message(mydat$gene_name[1])
	wa <- weighted.mean(x=mydat$end,w=mydat$mean_a_frac)
	wb <- weighted.mean(x=mydat$end,w=mydat$mean_b_frac)
	wa.gr <- GRanges(mydat$chr[1],IRanges(wa,wa))
	wb.gr <- GRanges(mydat$chr[1],IRanges(wb,wb))
	#g1 <- ggplot(mydat,aes(x=end,y=mean_a_frac)) + geom_bar(stat="identity") + goldmine:::ggnice() + scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1))
	a <- myta$a
	b <- myta$b
	myuse <- myta$use.frac[myta$res$tu==mytu,]
	myuse$pr <- mydat$pr
	myuse$chr <- mydat$chr
	myuse$start <- mydat$start
	myuse$end <- mydat$end
	myuse <- data.table(myuse)
	myuse <- melt(myuse,id.vars=c("pr","chr","start","end"))
	mysamp <- myta$samp.sub
	myuse$group <- mysamp[match(myuse$variable,mysamp$sample),]$group
	myuse$color <- mysamp[match(myuse$variable,mysamp$sample),]$color
	g1 <- ggplot(myuse,aes(x=end,y=value,color=group)) + geom_point() + goldmine:::ggnice() + scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1)) + scale_color_manual(values=unique(myuse$color))

	myuse <- data.table(myuse)
	myuse2 <- data.table(myuse)
	myuse2 <- myuse2[,list(end=end[1],mean=mean(value),sem=sd(value)/sqrt(length(value)),ci_upper=mean(value)+(1.96*sd(value)/sqrt(length(value))),ci_lower=mean(value)-(1.96*sd(value)/sqrt(length(value)))),by=c("pr","group")]

	g2 <- ggplot(myuse2,aes(x=end,y=mean,color=group,ymin=ci_lower,ymax=ci_upper)) + geom_point() + geom_errorbar(width=35) + goldmine:::ggnice() + scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1)) + scale_color_manual(values=unique(myuse$color))

	diffd <- tryCatch({
		myuse[,list(end=end[1],diff=mean(value[group==b])-mean(value[group==a]),diff_ci_lower=t.test(y=value[group==a],x=value[group==b])$conf.int[1],diff_ci_upper=t.test(y=value[group==a],x=value[group==b])$conf.int[2]),by=c("pr")]
	}, error = function(e) {
		myuse[,list(end=end[1],diff=mean(value[group==b])-mean(value[group==a]),diff_ci_lower=mean(value[group==b])-mean(value[group==a]),diff_ci_upper=mean(value[group==b])-mean(value[group==a])),by=c("pr")]
	
	})

	ww1 <- autoplot(wa.gr) + goldmine:::ggnice()
	ww2 <- autoplot(wb.gr) + goldmine:::ggnice()
	
	diffd$type <- "other"
	diffd[diffd$diff==min(diffd$diff),]$type <- "minor"
	diffd[diffd$diff==max(diffd$diff),]$type <- "major"

	gdiff <- ggplot(diffd,aes(x=end,y=diff,ymin=diff_ci_lower,ymax=diff_ci_upper,color=type)) + geom_point() + geom_errorbar(width=35) + goldmine:::ggnice() + geom_hline(yintercept=0) + scale_color_manual(values=c(other="#000000",major="#ff7f00",minor="#6a3d9a"))

	fcd <- mydat
	fcd$type <- "other"	
	fcd[fcd$l2_b_over_a_frac==min(fcd[(mean_a_frac>=0.05)|(mean_b_frac>=0.05),]$l2_b_over_a_frac),]$type <- "minor"
	fcd[fcd$l2_b_over_a_frac==max(fcd[(mean_a_frac>=0.05)|(mean_b_frac>=0.05),]$l2_b_over_a_frac),]$type <- "major"

	gfc <- ggplot(fcd,aes(x=end,y=l2_b_over_a_frac,color=type)) + geom_point() + goldmine:::ggnice() + geom_hline(yintercept=0) + scale_color_manual(values=c(other="#000000",major="#ff7f00",minor="#6a3d9a"))

	print(tracks(g1,g2,ww1,ww2,gdiff,gfc,title=paste0(mydat$gene_name[1]," (",mydat$strand[1],"), width=",max(mydat$end)-min(mydat$start)," bp, TU=",mytu),heights=c(0.25,0.25,0.05,0.05,0.20,0.20)))
}

pdf_prefix = file.path(args$output_dir,'PlotSigTusIntersect_')
# Intersection set
for(i in 1:length(apa.sig))
{
	message(names(apa.sig)[i])
	
	theres <- apa.sig[[i]]$res[int_sig==TRUE,]
	# Sort by effect size so biggest come up first
	mytus <- unique(theres[order(abs(theres$b_minus_a_frac),decreasing=TRUE),]$tu)

	message("Sig TUs: ",length(mytus))
	#pdf(file=paste0("output/PlotSigTusIntersect_",names(apa.sig)[i],".pdf"),width=10.5,height=8)
	pdf(file=paste0(pdf_prefix,names(apa.sig)[i],".pdf"),width=10.5,height=8)
	null <- lapply(mytus,function(x) plotOne(mytu=x,myta=apa.sig[[i]]))
	dev.off()
}

# Save CSVs - only with APA genes
output_dir1 = sprintf('%s/',args$output_dir)
for(i in 1:length(apa.sig))
{
	message(names(apa.sig)[i])

	theres <- apa.sig[[i]]$res[int_sig==TRUE,]
	# Sort by effect size so biggest come up first
	mytus <- unique(theres[order(abs(theres$b_minus_a_frac),decreasing=TRUE),]$tu)

	myres <- theres[tu %in% mytus,]

	#write.csv(myres,row.names=F,file=paste0("output/",names(apa.sig)[i],"_APA_TUs.csv"))
	write.csv(myres,row.names=F,file=paste0(output_dir1,names(apa.sig)[i],"_APA_TUs.csv"))
}
