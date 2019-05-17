library(data.table)
library(ggplot2)
library(stringr)
library(handy)
library(argparse)

parser <- ArgumentParser(description='a_run_hist_plot')

parser$add_argument("-i", "--input_dir", type="character", required=TRUE,
										dest="input_dir",
										help="a directory of mismatched index count files")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output directory")
args <- parser$parse_args()

hists <- dir(args$input_dir,full.names=T)
dat <- lapply(hists,fread)
names <- str_replace(str_replace(hists,".*/",""),"_Barcode.*","")
names(dat) <- names

for(i in 1:length(dat))
{
	dat[[i]]$File <- names(dat)[i]
}
dat1 <- rbindlist(dat)

dat2 <- dat1[,list(Value=rep(LongestARunInSequence,Count)),by="File"]

boxes <- dat2[,list(ymin=quantile(Value,0.05),low=as.double(quantile(Value)[2]),mid=as.double(median(Value)),top=quantile(Value)[4],ymax=quantile(Value,0.95)),by="File"]
boxes$run <- str_split_fixed(boxes$File,"_R1_",2)[,1]
boxes$barcode <- sapply(strsplit(str_split_fixed(boxes$File,"_",4)[,3],""),function(x) x[1])
boxes$barcode <- str_split_fixed(boxes$File,"_",4)[,3]
#boxes$barcode <- factor(boxes$barcode,levels=c("","","","OTHER"))

pdf(file=sprintf("%s/LongestARunInSequenceHistograms_1Mismatch.pdf",args$output_dir),width=10.5,height=8)
ggplot(boxes,aes(x=barcode,ymin=ymin,lower=low,middle=mid,upper=top,ymax=ymax, fill=run)) + geom_boxplot(stat="identity") + handy::ggnice() + scale_y_continuous(breaks=round(seq(min(dat1$LongestARunInSequence), max(dat1$LongestARunInSequence), by = 2),1)) + labs(title="All Runs/Barcodes",x="Barcode",y="Longest A Run in Sequence")
for(i in 1:length(dat))
{
	message(i)
	print(ggplot(dat[[i]],aes(x=LongestARunInSequence,y=Count)) + geom_bar(stat="identity") + handy::ggnice() + labs(title=names(dat)[i],x="Longest A Run in Sequence",y="# of Sequences") + scale_x_continuous(breaks=round(seq(min(dat1$LongestARunInSequence), max(dat1$LongestARunInSequence), by = 2),1)))
}
dev.off()
