# Save all the parts we'd like to have for the purposes of screening the APA calls
library(reshape)
library(argparse)

if (TRUE) {

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

parser <- ArgumentParser(description='testbatchadj')

parser$add_argument("-i", "--input_file", type="character", required=TRUE,
										dest="input_file",
										help="rscript workspace rd file")

parser$add_argument("-w", "--prepropa_rd", type="character", required=TRUE,
										dest="prepropa_rd",
										help="rscript workspace rd (prpropa_rd) file")

parser$add_argument("-s", "--chrom_size_file", type="character", required=TRUE,
										dest="chrom_size_file",
										help="chrom_size_file")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output dir")
args <- parser$parse_args()
}

#load("output/prepropa.rd")
load(args$prepropa_rd)
#load("output/apa.sig.rd")
load(args$input_file)

# Save extra tracks
mypr <- jtu$prs
join <- jtu$join
mypr$has_tu <- mypr$pr %in% join$pr

# All TUs
tu <- red_ens$tu
tu$nopr <- !(tu$tu %in% join$tu)
tu$uniq <- tu$tu %in% unique(join[unique_tu==TRUE,]$tu)
tu$color <- ""
tu[tu$nopr==TRUE]$color <- "153,153,153"
tu[(tu$nopr==FALSE)&(tu$uniq==TRUE)]$color <- "77,175,74"
tu[(tu$nopr==FALSE)&(tu$uniq==FALSE)]$color <- "228,26,28"
table(tu$color)

bb <- data.table(chrom=as.vector(seqnames(tu)),chromStart=as.integer(start(tu)-1),chromEnd=as.integer(end(tu)),name=tu$tu,score=0,strand=as.vector(strand(tu)),thickStart=start(tu)-1,thickEnd=end(tu),itemRgb=tu$color)
#bb$name <- paste0(bb$strand,bb$name)
track_dir = file.path(args$output_dir,'tracks')
if (!dir.exists(track_dir)){dir.create(track_dir, showWarnings=FALSE)}

#write.table(bb,row.names=FALSE,file="output/tracks/TU.bed",sep="\t",quote=F,col.names=F)
write.table(bb,row.names=FALSE,file=file.path(track_dir,"TU.bed"),sep="\t",quote=F,col.names=F)

# sort -k1,1 -k2,2n TU.bed > TU.sort.bed
# bedToBigBed TU.sort.bed ~/lustre/Common/ChromSizes.hg19.txt TU.bb

# All PRs
mypr$uniq <- mypr$pr %in% join[unique_pr==TRUE,]$pr
mypr$color <- ""
mypr[mypr$has_tu==FALSE]$color <- "153,153,153"
mypr[(mypr$has_tu==TRUE)&(mypr$uniq==TRUE)]$color <- "77,175,74"
mypr[(mypr$has_tu==TRUE)&(mypr$uniq==FALSE)]$color <- "228,26,28"
table(tu$color)
bbp <- data.table(chrom=as.vector(seqnames(mypr)),chromStart=as.integer(start(mypr)-1),chromEnd=as.integer(end(mypr)),name=mypr$pr,score=0,strand=as.vector(strand(mypr)),thickStart=start(mypr)-1,thickEnd=end(mypr),itemRgb=mypr$color)
bbp$name <- paste0(bbp$strand,bbp$name)
#write.table(bbp,row.names=FALSE,file="output/tracks/PR.bed",sep="\t",quote=F,col.names=F)
write.table(bbp,row.names=FALSE,file=file.path(track_dir,"PR.bed"),sep="\t",quote=F,col.names=F)

# All ends
ends <- pr$ends$allends
bbe <- data.table(chrom=as.vector(seqnames(ends)),chromStart=as.integer(start(ends)-1),chromEnd=as.integer(end(ends)),name="",score=0,strand=as.vector(strand(ends)),thickStart=as.integer(start(ends)-1),thickEnd=as.integer(end(ends)),itemRgb="0,0,0")
bbe$name <- bbe$strand
#write.table(bbe,row.names=FALSE,file="output/tracks/AllEnds.bed",sep="\t",quote=F,col.names=F)
write.table(bbe,row.names=FALSE,file=file.path(track_dir,"AllEnds.bed"),sep="\t",quote=F,col.names=F)

apa.d <- apa.sig

# Usage BedGraphs
# Make big meta-matrix of all the % usages
# Can't because different rows
#hd <- apa.d[[1]]

#hd <- apa.d[[1]]$use.frac
#mel <- melt(as.data.frame.array(hd))

#hd <- cbind(hd,apa.d[[2]]$use.frac[,5:8]) #TODO/BUG: the number of columns do not match!!!

track_prefix = file.path(track_dir,"Frac")

# Start by making diff tracks for each comp. May want more detail but can at least start here
difftr <- lapply(1:length(apa.d),function(x){
	myres <- apa.d[[x]]$res
	a <- apa.d[[x]]$a
	b <- apa.d[[x]]$b
	rm <- myres[strand=="-",]
	rp <- myres[strand=="+",]
	
	bp <- with(rp,data.table(chr=chr,start=as.integer(start-1),end=as.integer(end),value=b_minus_a_frac))
	bm <- with(rm,data.table(chr=chr,start=as.integer(start-1),end=as.integer(end),value=b_minus_a_frac))

	bp <- bp[order(bp$chr,bp$start),]
	bm <- bm[order(bm$chr,bm$start),]

	#pp <- paste0("output/tracks/Frac",b,"-Frac",a,"_Plus.bg")
	pp <- paste0(track_prefix,b,"-Frac",a,"_Plus.bg")
	#pp2 <- paste0("output/tracks/Frac",b,"-Frac",a,"_Plus.bw")
	pp2 <- paste0(track_prefix,b,"-Frac",a,"_Plus.bw")
	#pm <- paste0("output/tracks/Frac",b,"-Frac",a,"_Minus.bg")
	pm <- paste0(track_prefix,b,"-Frac",a,"_Minus.bg")
	#pm2 <- paste0("output/tracks/Frac",b,"-Frac",a,"_Minus.bw")
	pm2 <- paste0(track_prefix,b,"-Frac",a,"_Minus.bw")

	write.table(bp,file=pp,row.names=F,col.names=F,quote=F)
	write.table(bm,file=pm,row.names=F,col.names=F,quote=F)

	message("Saving :",pp2)
	message("Saving :",pm2)

	#system(paste("bedGraphToBigWig",pp,"~/lustre/Common/ChromSizes.hg19.txt",pp2))
	system(paste("bedGraphToBigWig",pp,args$chrom_size_file,pp2))
	#system(paste("bedGraphToBigWig",pm,"~/lustre/Common/ChromSizes.hg19.txt",pm2))
	system(paste("bedGraphToBigWig",pm,args$chrom_size_file,pm2))
	NULL
})

outdir_base1 <- sprintf("%s/",dirname(args$output_dir))

# Save tables
summary(apa.d)
lapply(1:length(apa.d),function(x){
	#write.csv(apa.d[[x]]$res,row.names=F,file=paste0("output/",names(apa.d)[x],"_APA.csv"))
	write.csv(apa.d[[x]]$res,row.names=F,file=paste0(outdir_base1,names(apa.d)[x],"_APA.csv"))
})
