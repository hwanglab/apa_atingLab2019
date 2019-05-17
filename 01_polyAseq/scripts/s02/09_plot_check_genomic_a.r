library(data.table)
library(stringr)
library(ggplot2)
library(reshape)
library(handy)
library(argparse)

parser <- ArgumentParser(description='plot_check_genomic_a')

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 1,
										help="number of cpus to utilize [1]")

parser$add_argument("-i", "--input_dir", type="character", required=TRUE,
										dest="input_dir",
										help="input directory")

parser$add_argument("-r", "--col_of_interest", type="character", required=FALSE, default = "7of10or6",
										dest="col_of_interest",
										help="column_of_interest[7of10or6]")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output directory")
args <- parser$parse_args()

gettab <- function(dirs,run)
{
	f <- lapply(dirs,fread)
	names <- dirs	
	for(i in 1:length(f))
	{
		f[[i]]$sample <- str_replace(str_replace(names[i],".*/",""),"\\..*","")
	}
	tab <- rbindlist(f)
	tab$run <- run
	return(tab)
}

tabs <- list()

tabs[[1]] <- gettab(dirs=dir(args$input_dir,pattern=".csv",full.names=T),run="7of10or6")
tab <- rbindlist(tabs)
tab$cell <- str_split_fixed(tab$sample,"_",2)[,1]

chk <- tab[,list(n=length(cell)),by="run"]$n
stopifnot(chk==max(chk))

# --------------------------------------------------------------------
# Table of genomicA vs noGenomicA
adf <- tab[,list(nReads=sum(count)),by=c("run","sample","group")]

adfc <- cast(adf,formula="run+sample~group",value="nReads")
adfc$frac_no <- adfc$noGenomicA/(adfc$genomicA+adfc$noGenomicA)

#write.csv(adfc,row.names=F,file="output/GenomicA_Fractions.csv")
write.csv(adfc,row.names=F,file=file.path(args$output_dir, "GenomicA_Fractions.csv"))

adfc$cell <- str_split_fixed(adfc$sample,"_",2)[,1]

#pdf(file="output/GenomicA_Fractions.pdf",width=10.5,height=8)
pdf(file=file.path(args$output_dir,"GenomicA_Fractions.pdf"),width=10.5,height=8)
ggplot(adfc[!str_detect(adfc$run,"flank100"),],aes(x=sample,y=frac_no,fill=run)) + geom_bar(stat="identity",position="dodge") + handy::ggnice() + coord_flip() + scale_y_continuous(breaks=seq(0.0,1.05,by=0.05))
dev.off()

#stat <- fread("./output/AlignmentFullStats.csv")
stat <- fread(file.path(args$output_dir,"AlignmentFullStats.csv"))
adfc$total <- adfc$genomicA+adfc$noGenomicA

stopifnot(adfc$total==stat[match(adfc$sample,stat$sample),]$MapToAsmChr)

# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Table of hexamer counts

hdf <- tab[,list(hexamer=hexamer,count=count,frac=count/sum(count)),by=c("run","sample","group")]

hdf$hashex <- FALSE
hdf[hexamer!="NONE",]$hashex <- TRUE

hdfa <- hdf[,list(count=sum(count),frac=sum(frac)),by=c("run","sample","group","hashex")]

#pdf(file="output/GenomicA_Hexamers.pdf",width=10.5,height=8)
pdf(file=file.path(args$output_dir,"GenomicA_Hexamers.pdf"),width=10.5,height=8)

# Hexamer/NoHexamer
for(r in unique(hdf$run))
{
	print(ggplot(hdfa[(run==r)&(hashex==TRUE),],aes(x=sample,y=frac,fill=group)) + geom_bar(stat="identity",position="dodge") + handy::ggnice() + coord_flip() + labs(title=r))
}
dev.off()

# Which Hexamer
cols <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#969696")

df <- hdf
df$hexamer <- factor(df$hexamer,levels=c("AATAAA","ATTAAA","AGTAAA","TATAAA","CATAAA","GATAAA","AATATA","AATACA","AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA","NONE"))
df <- df[order(df$hexamer),]

#pdf(file="output/GenomicA_HexamerFracs.pdf")
pdf(file=file.path(args$output_dir,"GenomicA_HexamerFracs.pdf"))
for(r in unique(hdf$run))
{
	print(ggplot(df[(hexamer!="NONE")&(run==r)&(group=="noGenomicA"),],aes(x=sample,y=frac,fill=hexamer)) + geom_bar(stat="identity") + handy::ggnice() + coord_flip() + labs(title=paste0("run: ",r,"; noGenomicA;")) + scale_fill_manual(values=cols))
}
dev.off()
# --------------------------------------------------------------------

