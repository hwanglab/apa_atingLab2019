library(data.table)
library(stringr)
library(ggplot2)
library(parallel)
library(reshape)
library(handy)

# Get stats

getaln <- function(x)
{
	cmd <- paste0("samtools view -F 4 ",x," | wc -l")
	message(cmd)
	as.numeric(system(cmd,intern=TRUE))
}

getunaln <- function(x)
{
	cmd <- paste0("samtools view -f 4 ",x," | wc -l")
	message(cmd)
	as.numeric(system(cmd,intern=TRUE))
}

library(argparse)

parser <- ArgumentParser(description='plot_dups')

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 1,
										help="number of cpus to utilize [1]")

parser$add_argument("-i", "--input_dir", type="character", required=TRUE,
										dest="input_dir",
										help="dedupped BAM files directory")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="raw BAM files directory")


args <- parser$parse_args()

ncore <- args$ncpu

raw_bam_dir <- args$output_dir
message(raw_bam_dir)

#bams_sorted <- dir("./output",pattern="01_sort.bam$",full.names=T)
bams_sorted <- dir(raw_bam_dir,pattern="01_sort.bam$",full.names=T)
message(bams_sorted)

sorted_aln <- mclapply(bams_sorted,getaln,mc.cores=ncore)
sorted_un <- mclapply(bams_sorted,getunaln,mc.cores=ncore)
dt <- data.table(bam=bams_sorted,class="SortedAligned",nReads=do.call(c,sorted_aln))
dt <- rbind(dt,data.table(bam=bams_sorted,class="SortedUnAligned",nReads=do.call(c,sorted_un)))

#bams_mapq <- dir("./output",pattern="02_mapq20.bam$",full.names=T)
bams_mapq <- dir(raw_bam_dir,pattern="02_mapq20.bam$",full.names=T)

sorted_aln <- mclapply(bams_mapq,getaln,mc.cores=ncore)
dt <- rbind(dt,data.table(bam=bams_sorted,class="MapqGr20",nReads=do.call(c,sorted_aln)))

#bams_dedup <- dir("./output/03_dedup_bams",pattern=".bam$",full.names=T)
bams_dedup <- dir(args$input_dir,pattern=".bam$",full.names=T)

message(bams_dedup)

dedup_counts <- mclapply(bams_dedup,getaln,mc.cores=ncore)
dt <- rbind(dt,data.table(bam=bams_dedup,class="DePcrDup",nReads=do.call(c,dedup_counts)))

getauto <- function(x)
{
	cmd <- paste0("samtools idxstats ",x," | cut -f 1,3")
	message(cmd)
	idx <- str_split_fixed(system(cmd,intern=TRUE),"\t",2)
	idx <- data.table(idx)
	idx$V2 <- as.numeric(idx$V2)
	sum(idx[idx$V1 %in% handy::chrs(),]$V2)
} 
auto_counts <- mclapply(bams_dedup,getauto,mc.cores=ncore)

dt <- rbind(dt,data.table(bam=bams_dedup,class="AsmGenome",nReads=do.call(c,auto_counts)))

dt$bam <- str_replace(str_replace(dt$bam,".*/",""),"\\..*","")

dtc <- cast(dt,formula="bam~class",value="nReads")

out <- data.table(sample=dtc$bam,TotalToBowtie=dtc$SortedAligned+dtc$SortedUnAligned,Aligned=dtc$SortedAligned,UniqueMap=dtc$MapqGr20,AfterPcrDeDup=dtc$DePcrDup,MapToAsmChr=dtc$AsmGenome)

per <- data.table(apply(out[,c(-1),with=F],2,function(x) x/out$TotalToBowtie))
names(per) <- paste0(names(per),"_FracTotal")

tab <- cbind(out,per)
#write.csv(tab,row.names=F,file="output/AlignmentFullStats.csv")
write.csv(tab,row.names=F,file=file.path(raw_bam_dir,"AlignmentFullStats.csv"))

df <- melt(tab)
#df <- df[!str_detect(variable,"_FracTotal"),]

df$variable <- factor(df$variable,levels=unique(df$variable))

df$cell <- str_split_fixed(df$sample,"_",2)[,1]

#pdf(file="output/AlignmentFullStats.pdf",width=10.5,height=8)
pdf(file=file.path(raw_bam_dir,"AlignmentFullStats.pdf"),width=10.5,height=8)

for(lev in levels(df$variable))
{
	print(ggplot(df[variable==lev,],aes(x=sample,y=value,fill=cell)) + geom_bar(stat="identity") + handy::ggnice() + coord_flip() + labs(title=lev,y="# of Reads"))
}
dev.off()

df2 <- df[!str_detect(variable,"_FracTotal"),]
#pdf(file="output/AlignmentFullStatsByCell.pdf",width=10.5,height=8)
pdf(file=file.path(raw_bam_dir,"AlignmentFullStatsByCell.pdf"),width=10.5,height=8)
for(mycell in unique(df2$cell))
{
	print(ggplot(df2[cell==mycell,],aes(x=variable,y=value,fill=sample)) + geom_bar(stat="identity") + handy::ggnice() + coord_flip() + labs(title=mycell,y="# of Reads") + facet_wrap(~sample,ncol=2))
}
dev.off()
