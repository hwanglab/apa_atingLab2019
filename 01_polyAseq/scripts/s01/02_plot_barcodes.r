#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(stringr)
library(ggplot2)
library(e1071)
library(matrixStats)
library(handy)

parser <- ArgumentParser(description='plot_barcodes')

parser$add_argument("-i", "--demux_dir", required=TRUE,
										dest="demux_dir",
										help="demultiplexed FASTQ files directory")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

hists <- dir(args$demux_dir,pattern="BarcodeCounts",full.names=T)

dat <- lapply(hists,fread)
names <- str_replace(str_replace(hists,".*/",""),"_Barcode.*","")
names(dat) <- names

for(i in 1:length(dat))
{
	dat[[i]]$File <- names(dat)[i]
}
dat <- rbindlist(dat)

# File-level summary of exact matches
agg <- dat[,list(n=sum(Count)),by=c("File","WantedBarcode")]
agg <- agg[,list(WantedBarcode=WantedBarcode,n=n,frac=n/sum(n)),by=c("File")]

write.csv(agg,row.names=F,file=sprintf("%s/SplitBarcodesExactFractions.csv",args$demux_dir))

# Mismatch Summary
# Annotate the dat with distances to each desired barcode
codes <- unique(dat$WantedBarcode)
codes <- codes[codes!="OTHER"]

for(c in codes)
{
	dat[,eval(paste0("MismatchFrom",c)):=sapply(dat$IncludedBarcode,function(x) hamming.distance(strsplit(c,"")[[1]],strsplit(x,"")[[1]]))]
}

spl <- split(dat,dat$File)

srt <- lapply(spl,function(x) x[order(x$Count,decreasing=T)])
srt <- rbindlist(srt)

srt <- srt[,list(WantedBarcode=WantedBarcode,IncludedBarcode=IncludedBarcode,Count=Count,Frac=Count/sum(Count),MismatchFromTGACCA=MismatchFromTGACCA,MismatchFromGCCAAT=MismatchFromGCCAAT,MismatchFromCTTGTA=MismatchFromCTTGTA,MismatchFromAAGTGC=MismatchFromAAGTGC),by="File"]
write.csv(srt,row.names=F,file=sprintf("%s/SplitBarcodesOtherFractions.csv",args$demux_dir))

# Check that 3 mismatch can uniquely call
mat <- as.matrix(srt[,str_detect(colnames(srt),"MismatchFrom"),with=F])

srt$n0 <- rowSums(mat==0)
srt$n1 <- rowSums(mat==1)
srt$n2 <- rowSums(mat==2)
srt$n3 <- rowSums(mat==3)
srt$n4 <- rowSums(mat==4)
srt$n5 <- rowSums(mat==5)
srt$n6 <- rowSums(mat==6)

apply(mat,1,function(x) data.table(sum(x==0)))


# Mismatch Distributions
srt$min <- rowMins(as.matrix(srt[,str_detect(colnames(srt),"MismatchFrom"),with=F]))


srt <- srt[,list(Freq=sum(Count)),by=c("File","min")]

#spl <- split(srt,srt$File)
#spl <- lapply(spl,function(x) as.data.frame(table(x$min)))
#for(i in 1:length(spl))
#{
#	spl[[i]]$File <- names(spl)[i]
#}
#mm <- rbindlist(spl)

pdf(file=sprintf("%s/SplitBarcodesMismatchFreq.pdf",args$demux_dir),width=10.5,height=8)
ggplot(srt,aes(x=factor(min),y=Freq)) + geom_bar(stat="identity") + handy::ggnice() + facet_wrap(~File,scales="free_x") + labs(x="Minimum Mismatches from a Barcode",y="Number of Sequences")
dev.off()

m2 <- srt[,list(MinimumMismatches=min,Freq=Freq,Frac=Freq/sum(Freq)),by="File"]

write.csv(m2,row.names=F,file=sprintf("%s/SplitBarcodesMinMismatches.csv",args$demux_dir))
