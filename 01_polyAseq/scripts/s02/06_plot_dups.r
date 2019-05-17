library(data.table)
library(stringr)
library(ggplot2)
library(parallel)
library(reshape)
library(handy)
library(argparse)

parser <- ArgumentParser(description='plot_dups')

parser$add_argument("-f", "--sample_sheet", type="character", required=TRUE,
										dest="sample_sheet",
										help="fastq_index.csv")

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 1,
										help="number of cpus to utilize [1]")

parser$add_argument("-i", "--input_dir", type="character", required=TRUE,
										dest="input_dir",
										help="dedupped BAM files directory")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output directory")
args <- parser$parse_args()

#samps <- fread("../01_ProcessSeq/input/fastq_index.csv")$sample
samps <- fread(args$sample_sheet)$sample

#files <- lapply(samps,function(x) dir("./output/03_dedup_bams",pattern=paste0(x,".*.csv"),full.names=T))
files <- lapply(samps,function(x) dir(args$input_dir,pattern=paste0(x,".*.csv"),full.names=T))
names(files) <- samps
dat <- lapply(files,function(x) lapply(x,fread))

files2 <- lapply(files,function(x) str_replace(str_replace(x,".*_",""),".csv",""))

for(i in 1:length(dat))
{
	names(dat[[i]]) <- files2[[i]]
}

lapply(dat[[1]],colnames)


# Dup count hist for each class by each sample
#pdf(file="output/PcrDuplicationLevelsByScheme.pdf",width=10.5,height=8)
pdf(file=file.path(args$output_dir, "PcrDuplicationLevelsByScheme.pdf"),width=10.5,height=8)

for(i in 1:length(dat))
{
	cell <- names(dat)[i]
	message(cell)
	myl <- dat[[i]]
	hist <- rbindlist(lapply(names(myl),function(x) {mine <- myl[[x]]; mine[,list(class=x,chr=chr,pos=pos,strand=strand,n=bamrow)];}))
	h1 <- hist[,list(n=1:10,reads=sapply(1:10,function(x) sum(n[n==x]))),by="class"]
	h2 <- hist[,list(n=">10",reads=sum(n[n>10])),by="class"]
	df <- rbind(h1,h2)
	df$n <- factor(df$n,levels=c(1:10,">10"))

	cols <- c(ByCPS="#a6cee3",ByCPLS="#1f78b4",ByCPST="#fb9a99",ByCPLST="#e31a1c")
	#ggplot(df,aes(x=n,y=reads,fill=class)) + geom_bar(stat="identity") + facet_wrap(~class,ncol=2) + handy::ggnice() + scale_fill_manual(values=cols)
	#ggplot(df,aes(x=n,y=reads,fill=class)) + geom_bar(stat="identity",position="dodge") + handy::ggnice() + scale_fill_manual(values=cols)
	print(ggplot(df,aes(x=n,y=reads,color=class,group=class)) + geom_line() + geom_point() + handy::ggnice() + scale_color_manual(values=cols) + labs(title=paste0("Duplication Levels for ",cell)))
}
dev.off()

# Table of remaining reads after each scheme
out <- mclapply(names(dat),function(x){
	myl <- dat[[x]]
	cell <- x
	message(cell)
	hist <- rbindlist(lapply(names(myl),function(x) {mine <- myl[[x]]; mine[,list(class=x,chr=chr,pos=pos,strand=strand,n=bamrow)];}))
	cast(hist[,list(sample=cell,n=length(n)),by="class"],formula="sample~class",value="n")
},mc.cores=args$ncpu)
out <- rbindlist(out)

#write.csv(out,row.names=F,file="output/ReadsAfterDeDupByScheme.csv")
write.csv(out,row.names=F,file=file.path(args$output_dir, "ReadsAfterDeDupByScheme.csv"))

#pdf(file="output/ReadsAfterDeDupByScheme.pdf",width=10.5,height=10.5)
pdf(file=file.path(args$output_dir, "ReadsAfterDeDupByScheme.pdf"),width=10.5,height=10.5)
df <- melt(out)
ggplot(df,aes(x=sample,y=value,fill=variable)) + geom_bar(stat="identity") + handy::ggnice() + scale_fill_manual(values=cols) + facet_wrap(~variable,ncol=2) + coord_flip() + labs(title="Reads Left After DeDuplication",y="# of Reads")
dev.off()

