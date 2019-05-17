library(data.table)
library(stringr)
library(ggplot2)
library(reshape)
library(handy)
library(argparse)

parser <- ArgumentParser(description='a_chop_plot')

parser$add_argument("-i", "--input_dir", type="character", required=TRUE,
										dest="input_dir",
										help="a directory of mismatched index count files")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output directory")
args <- parser$parse_args()

# Load histogram CSVs
#a_hists <- dir("./output/chop_a_tails",pattern="AsInWindow.csv",full.names=T)
a_hists <- dir(args$input_dir,pattern="AsInWindow.csv",full.names=T)
dat <- lapply(a_hists,fread)
names <- str_replace(str_replace(a_hists,".*/",""),"_AsInWindow.*","")
names(dat) <- names
a_dat <- dat

#t_hists <- dir("./output/chop_a_tails",pattern="TagLenghAfterACut.csv",full.names=T)
t_hists <- dir(args$input_dir,pattern="TagLenghAfterACut.csv",full.names=T)
dat <- lapply(t_hists,fread)
names <- str_replace(str_replace(t_hists,".*/",""),"_TagLenghAfter.*","")
names(dat) <- names
t_dat <- dat

# Compute total numbers table for each run
getClasses <- function(n,a,t)
{
	message(n)
	countNoA <- sum(a[(AsInWindow==-1)|(AsInWindow<9),]$Count)
	countAllA <- t[Length==0,]$Count
	countCutA <- sum(a[AsInWindow>=9,]$Count) - countAllA
	data.table(file=n,countCutA=countCutA,countAllA=countAllA,countNoA=countNoA)
}
summary <- rbindlist(lapply(names(a_dat),function(x) getClasses(n=x,a=a_dat[[x]],t=t_dat[[x]])))
#write.csv(summary,row.names=F,file="output/ChopATailsResults.csv")
write.csv(summary,row.names=F,file=sprintf("%s/ChopATailsResults.csv",args$output_dir))

# Plot the summary proportions
df <- melt(summary,id.vars="file")
df$file <- factor(df$file,levels=rev(unique(df$file)))
df <- df[,list(variable=variable,value=value,frac=value/sum(value)),by="file"]

#pdf(file="output/ChopATailsResults.pdf",width=10.5,height=8)
pdf(file=sprintf("%s/ChopATailsResults.pdf",args$output_dir),width=10.5,height=8)
ggplot(df,aes(x=file,y=value,fill=variable)) + geom_bar(stat="identity") + handy::ggnice() + coord_flip() + labs(title="Number of Reads with Each Cutting Decision",x="FASTQ File",y="Number of Reads") + scale_fill_manual(values=c("#4daf4a","#377eb8","#e41a1c"))
ggplot(df,aes(x=file,y=frac,fill=variable)) + geom_bar(stat="identity") + handy::ggnice() + coord_flip() + labs(title="Proportion of Reads with Each Cutting Decision",x="FASTQ File",y="Fraction of Total Reads") + scale_fill_manual(values=c("#4daf4a","#377eb8","#e41a1c")) + scale_y_continuous(breaks=seq(0,1,by=0.1))
dev.off()


# Plot the score distribution
for(i in 1:length(a_dat))
{
	a_dat[[i]]$File <- names(a_dat)[i]
}
adt <- rbindlist(a_dat)
adt <- adt[,list(AsInWindow=AsInWindow,Count=Count,Frac=Count/sum(Count)),by="File"]
adt$File <- factor(adt$File,levels=rev(unique(adt$File)))
#pdf(file="output/AsIn10bpWindow.pdf",width=10.5,height=8)
pdf(file=sprintf("%s/AsIn10bpWindow.pdf",args$output_dir),width=10.5,height=8)
ggplot(adt,aes(x=AsInWindow,y=File,fill=Frac)) + geom_tile() + handy::ggnice() + scale_x_continuous(breaks=seq(-1,10,by=1)) + scale_fill_gradient2(low="#3288bd", mid="#e6f598", high="#d53e4f") + labs(title="Distributions of Number of As in 10bp Window for all FASTQs")
for(i in 1:length(a_dat))
{
	print(ggplot(a_dat[[i]],aes(x=AsInWindow,y=Count)) + geom_bar(stat="identity") + handy::ggnice() + scale_x_continuous(breaks=seq(-1,10,by=1)) + labs(title=paste0("As in 10bp Window: ",names(a_dat)[i])))
}
dev.off()

# Plot the cut site density
for(i in 1:length(t_dat))
{
	t_dat[[i]]$File <- names(t_dat)[i]
}
tdt <- rbindlist(t_dat)
tdt <- tdt[,list(Length=Length,Count=Count,Frac=Count/sum(Count)),by="File"]
tdt$File <- factor(tdt$File,levels=rev(unique(tdt$File)))

#pdf(file="output/CutSiteDensity.pdf",width=10.5,height=8)
pdf(file=sprintf("%s/CutSiteDensity.pdf",args$output_dir),width=10.5,height=8)
ggplot(tdt,aes(x=Length,y=File,fill=Frac)) + geom_tile() + handy::ggnice() + scale_x_continuous(breaks=seq(min(tdt$Length),max(tdt$Length),by=5)) + scale_fill_gradient2(low="#3288bd", mid="#e6f598", high="#d53e4f") + labs(title="Distribution of Inferred Cleavage Site Location in Read for all FASTQs")
ggplot(tdt[!str_detect(tdt$File,"OTHER"),],aes(x=Length,y=File,fill=Frac)) + geom_tile() + handy::ggnice() + scale_x_continuous(breaks=seq(min(tdt$Length),max(tdt$Length),by=5)) + scale_fill_gradient2(low="#3288bd", mid="#e6f598", high="#d53e4f") + labs(title="Distribution of Inferred Cleavage Site Location in Read for non-other FASTQs")
for(i in 1:length(t_dat))
{
	print(ggplot(t_dat[[i]],aes(x=Length,y=Count)) + geom_point()  + geom_line() + handy::ggnice() + geom_vline(xintercept=18) + labs(title=paste0("Cleavage Site in Read Location Density: ",names(t_dat)[i])))
}
dev.off()
