library(Rsamtools)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape)
library(BSgenome.Hsapiens.UCSC.hg19)
library(matrixStats)
library(edgeR)
library(DEXSeq)

# --------------------------------------------------------------------
# Append BAM file paths into Experiment Sheet
appendBamFilePath <- function(exp_sheet_fn,bam_dir,bam_suffix=".bam")
{
    exp_sheet <- fread(file=exp_sheet_fn,sep=',',header = TRUE)
    bams <- dir(bam_dir,pattern=sprintf('%s$',bam_suffix),full.names=T)
    samples <- str_replace(str_replace(bams,".*/",""),bam_suffix,"")
    bam_dt <- data.table(sample=samples,bam=bams)
    return(merge(exp_sheet,bam_dt,by='sample',all.x=TRUE,in_place=TRUE))
}

appendBamFilePath2 <- function(exp_sheet_fn,bam_dir,bam_suffix=".bam")
{
  
  exp_sheet <- fread(file=exp_sheet_fn,sep=',',header = TRUE)
  bams <- dir(bam_dir,pattern=sprintf('%s$',bam_suffix),full.names=T)
  
  exp_sheet$bam <- lapply(exp_sheet$sample,
         function(x){
           bams[which(grepl(x, bams))]}
         )
  return(exp_sheet)
}

# --------------------------------------------------------------------
# Load in Sample CSV
readSampleInfo <- function(file,colors=NULL)
{
	if(is.null(file)){stop("Must provide a path to a file to open.")}
	if(!file.exists(file)){stop(paste("Could not open: ",file,sep=""))}
	samplesheet <- read.csv(file, header=TRUE, stringsAsFactors=FALSE)
	if(sum(c("sample","group","bam") %in% colnames(samplesheet))!=3){stop("CSV must contain columns: sample, group, bam")}
	nsamp <- nrow(samplesheet)
	ngroup <- length(unique(samplesheet$group))
	message(paste("Found ", nsamp, " samples in ", ngroup, " groups",sep=""))
	grouporder <- unique(samplesheet$group)
	samplesheet$group <- factor(samplesheet$group,levels=grouporder,ordered=TRUE)
	message(paste("Group order detected as: ", paste(grouporder,collapse=", ", sep=""), sep=""))
	message("Number of samples in each group:")
	print(summary(samplesheet$group))
	samplesheet <- samplesheet[order(samplesheet$group,samplesheet$sample),]
	if(is.null(colors) & sum(colnames(samplesheet)=="color")==0)
	{
		message("Auto-picking group colors from RColorBrewer")
	if(ngroup<=9)
	{
		colors <- RColorBrewer::brewer.pal(ngroup,"Set1")
	} else
	{
		colors <- colorRampPalette(brewer.pal(9,"Set1"))(ngroup)
	}
	samplesheet$color <- colors[as.numeric(samplesheet$group)]
	} else if(!is.null(colors))
	{
		if(length(colors)!=ngroup){stop("colors vector must equal number of groups in file")}
		samplesheet$color <- colors[as.numeric(samplesheet$group)]
	}
	samplesheet	
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Load BAMs to RAM
getReads <- function(samp, chrs, ncore, is_10x=FALSE)
{
  message("Reading BAMs from disk with ",ncore," concurrent processes")
  # Validate BAM existence
  # Check that BAM files exist
  check <- sapply(samp$bam,file.exists)
  if(sum(check)!=nrow(samp)){stop(paste0("BAM file(s) could not be found: ",toString(samp$bam[check==F])))}
  
  # Check BAIs exist
  check <- sapply(paste0(samp$bam,".bai"),file.exists)
  if(sum(check)!=nrow(samp)){stop(paste0("Bam Index (.bai) file(s) could not be found for: ",toString(samp$bam[check==F])))}
  
  # Validate asked for chrs are in the BAM
  b1 <- Rsamtools::BamFile(samp$bam[1])
  sl <- seqlengths(b1)
  if(sum(chrs %in% seqlevels(b1)) != length(chrs)){stop(paste0("Could not find chrs: ",toString(chrs[!(chrs %in% seqlevels(b1))]), " in given BAM header"))}
  
  # Use the index so we don't bother reading in from the unaligned chrs
  which.gr <- GRanges(chrs,IRanges(1,seqlengths(b1)[chrs]))
  
  bam2gr <- function(bampath)
  {
    message("Reading BAM file: ",bampath)
    # What: fields to read in (column filtering)
    # Flag: records to read in (row filtering)
    # Which: what sequences must be overlaped (chr/pos filtering)
    if (is_10x){
      param <- Rsamtools::ScanBamParam(what=character(), 
                                       which=which.gr, 
                                       flag=Rsamtools::scanBamFlag(
                                         isUnmappedQuery=FALSE,
                                         isMinusStrand=is_minus_strand),
                                       tag=c("CB","CR","UB","UR","BC"))
      bam.ga <- GenomicAlignments::readGappedReads(bampath, param = param)
    } else {
      param <- Rsamtools::ScanBamParam(what=character(),
                                       which=which.gr,
                                       flag=Rsamtools::scanBamFlag(
                                         isUnmappedQuery=FALSE,
                                         isPaired=FALSE,
                                         isProperPair=FALSE))
      bam.ga <- GenomicAlignments::readGAlignments(bampath, param = param)
    }
    
    bam.gr <- as(bam.ga, "GRanges")
    seqlevels(bam.gr) <- chrs
    seqlengths(bam.gr) <- sl[chrs]
    return(bam.gr)
  }
  reads <- mclapply(samp$bam,bam2gr,mc.cores=ncore)
  #reads <- lapply(samp$bam,bam2gr)
  # browser() #debug
  names(reads) <- samp$sample
  
  message("Read GRanges Size in Memory=",format(object.size(reads),units="auto"))
  return(reads)
}


# --------------------------------------------------------------------
# Create genomic coverage WIG files
saveCov <- function(samp,reads,path=".",bigwig=FALSE,server=NULL,ncore=1, seq_suffix='PolyA', track_suffix=NULL)
{
	if(!dir.exists(path))
	{
		dir.create(path,recursive=TRUE)
	}

	if((!is.null(server))&(bigwig==TRUE))
	{
		if (!is.null(track_suffix)) {
			track_fn <- sprintf("tracks_%s.txt",track_suffix)
		} else {
			track_fn <- sprintf("tracks.txt")
		}
		trackpath <- file.path(path,track_fn)
		
		message("Saving tracklist to ",trackpath)
		rgbs <- str_replace_all(apply(col2rgb(samp$color),2,toString)," ","")
		trackfile <- paste0("track type=bigWig name=\"",paste0(rep(samp$sample,each=2),sprintf("_%s ",seq_suffix),c("(+)","(-)")),"\" visibility=full maxHeightPixels=50 color=",rep(rgbs,each=2)," autoScale=off viewLimits=",c("0:30","-30:0")," bigDataUrl=",server,paste0(rep(samp$sample,each=2),"_",c("Plus","Minus")),".bw")
		trackfile <- c(trackfile,paste0("track type=bigWig name=\"",paste0(rep(samp$sample,each=2),sprintf("_%s ",seq_suffix),c("(+, Depth Scaled)","(-, Depth Scaled)")),"\" visibility=full maxHeightPixels=50 color=",rep(rgbs,each=2)," autoScale=off viewLimits=",c("0:30","-30:0")," bigDataUrl=",server,paste0(rep(samp$sample,each=2),"_",c("Plus_Scaled","Minus_Scaled")),".bw"))
		writeLines(trackfile,trackpath)
	}

	message("Dividing stranded information")
	reads.p <- mclapply(reads,function(x) x[strand(x)=="+"],mc.cores=ncore)
	reads.m <- mclapply(reads,function(x) x[strand(x)=="-"],mc.cores=ncore)

	message("Computing genome-wide coverage vectors")
	cov.p <- mclapply(reads.p,coverage,mc.cores=ncore)
	cov.m <- mclapply(reads.m,coverage,mc.cores=ncore)

	message("Assembling to GRanges")
	covgr.p <- mclapply(cov.p,function(x) as(x,"GRanges"),mc.cores=ncore)
	covgr.m <- mclapply(cov.m,function(x) as(x,"GRanges"),mc.cores=ncore)

	covgr.p <- lapply(covgr.p,function(x) x[x$score>0])
	covgr.m <- lapply(covgr.m,function(x) x[x$score>0])

	# Normalize via sizefactors
	message("Performing sizeFactor normalization via DESeq2")
	libsizes <- sapply(reads,length)
	cnt <- matrix(libsizes,nrow=1)
	colnames(cnt) <- names(libsizes)
	sizefactors <- estimateSizeFactorsForMatrix(cnt)

	covgr.p.norm <- covgr.p
	for(i in 1:length(covgr.p.norm))
	{
		covgr.p.norm[[i]]$score <- covgr.p.norm[[i]]$score/sizefactors[i]
	}
	covgr.m.norm <- covgr.m
	for(i in 1:length(covgr.m.norm))
	{
		covgr.m.norm[[i]]$score <- covgr.m.norm[[i]]$score/sizefactors[i]*-1
	}

	# Cast all minus strand to show up as negatives for viz purposes
	for(i in 1:length(covgr.m))
	{
		covgr.m[[i]]$score <- covgr.m[[i]]$score*-1
	}

	chrs <- seqlevels(reads[[1]])
	chrlens <- seqlengths(reads[[1]])

	writeBed <- function(gr,name)
	{
		message(paste(name,": Creating BedGraph",sep=""))
		filename <- file.path(path,paste(name,".bed",sep=""))
		values <- gr$score
		bed <- data.frame(chr=seqnames(gr), start=as.integer(start(gr)-1), end=as.integer(end(gr)), value=values)
		bed <- bed[bed$value!=0,]
		write.table(bed, file=filename, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

		if(bigwig==TRUE)
		{
			message(paste(name,": Converting to BigWig",sep=""))
			cs <- tempfile(name)
			write.table(data.frame(chrs,chrlens),col.names=F,row.names=F,file=cs, quote=F)
			cmd <- paste("wigToBigWig",filename,cs,file.path(path,paste(name,".bw",sep="")),sep=" ")
			message(cmd)
			system(cmd)
			file.remove(cs)
		}
	}
	doit <- mclapply(1:length(reads),function(x) writeBed(gr=covgr.p[[x]],name=paste0(names(reads)[x],"_Plus")))
	doit <- mclapply(1:length(reads),function(x) writeBed(gr=covgr.m[[x]],name=paste0(names(reads)[x],"_Minus")))
	doit <- mclapply(1:length(reads),function(x) writeBed(gr=covgr.p.norm[[x]],name=paste0(names(reads)[x],"_Plus_Scaled")))
	doit <- mclapply(1:length(reads),function(x) writeBed(gr=covgr.m.norm[[x]],name=paste0(names(reads)[x],"_Minus_Scaled")))

	return (NULL)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Create genomic coverage WIG files
saveCov_unstranded <- function(samp,reads,path=".",bigwig=FALSE,server=NULL,ncore=1,seq_suffix='PolyA',track_suffix=NULL)
{
  if(!dir.exists(path))
  {
    dir.create(path,recursive=TRUE)
  }
  
  if((!is.null(server))&(bigwig==TRUE))
  {
        
    if (!is.null(track_suffix)) {
			track_fn <- sprintf("tracks_%s.txt",track_suffix)
		} else {
			track_fn <- sprintf("tracks.txt")
		}
		trackpath <- file.path(path,track_fn)
		
    message("Saving tracklist to ",trackpath)
    rgbs <- str_replace_all(apply(col2rgb(samp$color),2,toString)," ","")
    trackfile <- paste0("track type=bigWig name=\"",
                        sprintf("%s_%s ",samp$sample,seq_suffix),
                        "\" visibility=full maxHeightPixels=70 color=",
                        rgbs,
                        " autoScale=off viewLimits=",
                        "0:60",
                        " bigDataUrl=",
                        server,
                        samp$sample,
                        ".bw")
    
    trackfile <- c(trackfile,
                   paste0("track type=bigWig name=\"",
                          sprintf("%s_%s (Depth Scaled)",samp$sample,seq_suffix),
                          "\" visibility=full maxHeightPixels=70 color=",
                          rgbs,
                          " autoScale=off viewLimits=",
                          "0:60",
                          " bigDataUrl=",
                          server,
                          sprintf("%s_Scaled",samp$sample),
                          ".bw")
                   )
    
    writeLines(trackfile,trackpath)
  }
  
 
  message("Computing genome-wide coverage vectors")
  cov <- mclapply(reads,coverage,mc.cores=ncore)

  message("Assembling to GRanges")
  covgr <- mclapply(cov,function(x) as(x,"GRanges"),mc.cores=ncore)
  covgr <- lapply(covgr,function(x) x[x$score>0])

  # Normalize via sizefactors
  message("Performing sizeFactor normalization via DESeq2")
  libsizes <- sapply(reads,length)
  cnt <- matrix(libsizes,nrow=1)
  colnames(cnt) <- names(libsizes)
  sizefactors <- estimateSizeFactorsForMatrix(cnt)
  
  covgr.norm <- covgr
  for(i in 1:length(covgr.norm))
  {
    covgr.norm[[i]]$score <- covgr.norm[[i]]$score/sizefactors[i]
  }

  chrs <- seqlevels(reads[[1]])
  chrlens <- seqlengths(reads[[1]])
  
  writeBed <- function(gr,name)
  {
    message(paste(name,": Creating BedGraph",sep=""))
    filename <- file.path(path,paste(name,".bed",sep=""))
    values <- gr$score
    bed <- data.frame(chr=seqnames(gr), start=as.integer(start(gr)-1), end=as.integer(end(gr)), value=values)
    bed <- bed[bed$value!=0,]
    write.table(bed, file=filename, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    
    if(bigwig==TRUE)
    {
      message(paste(name,": Converting to BigWig",sep=""))
      cs <- tempfile(name)
      write.table(data.frame(chrs,chrlens),col.names=F,row.names=F,file=cs, quote=F)
      cmd <- paste("wigToBigWig",filename,cs,file.path(path,paste(name,".bw",sep="")),sep=" ")
      message(cmd)
      system(cmd)
      file.remove(cs)
    }
  }
  doit <- mclapply(1:length(reads),function(x) writeBed(gr=covgr[[x]],name=paste0(names(reads)[x])))
  
  doit <- mclapply(1:length(reads),function(x) writeBed(gr=covgr.norm[[x]],name=paste0(names(reads)[x],"_Scaled")))
  
	return(NULL)
}

# --------------------------------------------------------------------
writeBEDFromGRanges <- function(gr, file, name=NULL, offset=0)
{
	if(file.exists(file)){file.remove(file)}
	fileConn<-file(file)
	writeLines(c(paste("track name=\"",file,"\"",sep="")), fileConn)
	close(fileConn)
	df <- data.frame(chr=as.character(seqnames(gr)),start=start(gr)-1+offset, end=end(gr)-offset)
	if(!is.null(name))
	{
		df$name <- values(gr)[,name]
	}
	#df <- df[df$chr %in% c(sapply(seq(1,22),function(x) paste("chr",x,sep="")),"chrX","chrY"),]
	write.table(df, file=file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Hexamer counting for any set of genomic ranges
countHexes <- function(gr,us=50)
{
	hexes <- c("AATAAA","ATTAAA","AGTAAA","TATAAA","CATAAA","GATAAA","AATATA","AATACA","AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA")

	message("Generating flanks")

	# Also including the PR itself in addition to the flank here
	fl <- resize(gr,us+width(gr),fix="end")
	#fl <- flank(re,width=us)

	message("Pulling sequence data from hg19")
	seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,names=fl)

	message("Searching hexamer occurances")
	pd <- PDict(x=hexes,max.mismatch=0)
	vc <- vcountPDict(pd, seq)
	vc <- t(vc)

	message("Selecting matches")
	ap <- apply(vc,1,function(x) which(x>=1)[1])
	hp <- hexes[ap]
	gr$hex <- hp
	gr[rowSums(vc)==0]$hex <- "NONE"
	#table(gr$hex)

	gr$hex <- factor(gr$hex,levels=c(hexes,"NONE"))

	return(gr)
}
# --------------------------------------------------------------------
# Processing Region Generator
makeProcReg2 <- function(samp,read,dist_r1r2=50,dist=10,cov=10,ncore=1)
{
  
  message("Generating 3' End Map")
  ends <- resize(read,1,fix='end')
  end(ends) <- end(ends) + dist_r1r2
  
  ends.uniq <- unique(ends)
  allends <- unique(do.call(c,unname(ends.uniq)))
  
  message("Reducing to PRs")
  red <- reduce(allends,min.gapwidth=dist)
  
  message("Getting Sample-Level PR Counts")
  co <- mclapply(ends,function(x) countOverlaps(red,x),mc.cores=ncore)
  counts <- do.call(cbind,co)
  
  # Convert to RPKMs
  message("Converting to RPKMs")
  libsizes <- sapply(reads,length)
  rpkm <- rpkm(counts,gene.length=width(red),lib.sizes=libsizes)
  #r1 <- rpkm[,1]
  #c1 <- counts[,1]
  #r2 <- ((10^9)*c1)/(libsizes[1]*as.numeric(width(pr)))
  
  # Get group means for plotting
  message("Generating Group Means")
  rpkmmeans <- lapply(unique(samp$group),function(x) rowMeans(rpkm[,samp$group==x]))
  rpkmmeans <- do.call(cbind,rpkmmeans)
  colnames(rpkmmeans) <- unique(samp$group)
  
  countmeans <- lapply(unique(samp$group),function(x) rowMeans(counts[,samp$group==x]))
  countmeans <- do.call(cbind,countmeans)
  colnames(countmeans) <- unique(samp$group)
  
  # Build universe
  message("Applying Coverage Filter to Build Universe")
  keep <- rowMaxs(countmeans)>=cov
  pr <- red[keep]	
  pr$pr <- paste0("PR",1:length(pr))
  counts <- counts[keep,]
  rownames(counts) <- pr$pr
  rpkm <- rpkm[keep,]
  rownames(rpkm) <- pr$pr
  countmeans <- countmeans[keep,]
  rownames(countmeans) <- pr$pr
  rpkmmeans <- rpkmmeans[keep,]
  rownames(rpkmmeans) <- pr$pr
  
  return(list(pr=pr,counts=list(raw=counts,rpkm=rpkm),means=list(raw=countmeans,rpkm=rpkmmeans),ends=list(allends=allends,reduced=red)))
}

# --------------------------------------------------------------------
# Processing Region Generator
makeProcReg <- function(samp,reads,dist=10,cov=10,ncore=1)
{
	message("Generating 3' End Map")
	ends <- mclapply(reads,function(x) resize(x,1,fix="end"),mc.cores=ncore)
	ends.uniq <- mclapply(ends,unique,mc.cores=ncore)
	allends <- unique(do.call(c,unname(ends.uniq)))

	message("Reducing to PRs")
	red <- reduce(allends,min.gapwidth=dist)

	message("Getting Sample-Level PR Counts")
	co <- mclapply(ends,function(x) countOverlaps(red,x),mc.cores=ncore)
	counts <- do.call(cbind,co)

	# Convert to RPKMs
	message("Converting to RPKMs")
	libsizes <- sapply(reads,length)
	rpkm <- rpkm(counts,gene.length=width(red),lib.sizes=libsizes)
	#r1 <- rpkm[,1]
	#c1 <- counts[,1]
	#r2 <- ((10^9)*c1)/(libsizes[1]*as.numeric(width(pr)))

	# Get group means for plotting
	message("Generating Group Means")
	rpkmmeans <- lapply(unique(samp$group),function(x) rowMeans(rpkm[,samp$group==x]))
	rpkmmeans <- do.call(cbind,rpkmmeans)
	colnames(rpkmmeans) <- unique(samp$group)

	countmeans <- lapply(unique(samp$group),function(x) rowMeans(counts[,samp$group==x]))
	countmeans <- do.call(cbind,countmeans)
	colnames(countmeans) <- unique(samp$group)

	# Build universe
	message("Applying Coverage Filter to Build Universe")
	keep <- rowMaxs(countmeans)>=cov
	pr <- red[keep]	
	pr$pr <- paste0("PR",1:length(pr))
	counts <- counts[keep,]
	rownames(counts) <- pr$pr
	rpkm <- rpkm[keep,]
	rownames(rpkm) <- pr$pr
	countmeans <- countmeans[keep,]
	rownames(countmeans) <- pr$pr
	rpkmmeans <- rpkmmeans[keep,]
	rownames(rpkmmeans) <- pr$pr

	return(list(pr=pr,counts=list(raw=counts,rpkm=rpkm),means=list(raw=countmeans,rpkm=rpkmmeans),ends=list(allends=allends,reduced=red)))
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
reduceGenes <- function(genes,chrs,flank=5000,ncore=1)
{
  
	genes <- makeGRanges(genes,strand=T)
	genes <- genes[seqnames(genes) %in% chrs]
	genes$cdsStart <- genes$cdsStart+1
	genes$cdsEnd <- genes$cdsEnd+1
	stopifnot(genes$strand %in% c("+","-"))

	# Some genes may have disjoint txUnits, let's treat them as separate genes to make things easier
	message("Clustering TUs")
	df <- as.data.frame(table(genes$gene.id))
	genes.singles <- genes[genes$gene.id %in% df[df$Freq==1,]$Var1]
	genes.tored <- genes[genes$gene.id %in% df[df$Freq>1,]$Var1]
	spl <- split(genes.tored,genes.tored$gene.id)
	spl.red <- lapply(spl,function(x) reduce(x,with.revmap=FALSE))#collect union boundary of the same gene
	for(i in 1:length(spl.red))
	{
		spl.red[[i]]$gene.id <- names(spl)[i]
	}
	
	genes.red <- do.call(c,unname(spl.red))
	genes.red$name <- genes[match(genes.red$gene.id,genes$gene.id)]$name
	message("Single isoform annotated: ",length(unique(genes.singles$gene.id)))
	message("Multiple isoform annotated: ",length(unique(genes.tored$gene.id)))	
	dt <- makeDT(genes.red)
	genes.red.counts <- dt[,list(nPos=length(strand),nStrands=length(unique(strand))),by="gene.id"]
	message("Reduced to disjoint clusters on same strand only: ",nrow(genes.red.counts[(nPos>1)&(nStrands==1),]))
	message("Reduced to disjoint clusters on different strands: ",nrow(genes.red.counts[(nStrands>1),]))

	# Assign isoforms to TUs
	tu <- genes.singles
	values(tu) <- NULL
	tu$gene.id <- genes.singles$gene.id
	tu$name <- genes.singles$name
	tu <- c(tu,genes.red)
	seqlevels(tu) <- chrs
	tu <- tu[order(tu)]
	tu$tu <- paste0("TU",1:length(tu))

	# From now on, must aggregate by TU rather than gene.id, as there can be duplication in gene.id	
	message("Assigning transcripts to TUs")
	genes.singles$tu <- tu[match(genes.singles$gene.id,tu$gene.id),]$tu

	fo <- data.table(as.data.frame(findOverlaps(genes.tored,tu)))#to see how multiple-isoform genes are overlapped with the tu-assigned set; as a default, strand information is also included in the criteria
	fo$gene_gid <- genes.tored[fo$queryHits]$gene.id
	fo$tu_gid <- tu[fo$subjectHits]$gene.id
	fo$assign <- fo$gene_gid==fo$tu_gid
	fo <- fo[assign==TRUE,] #For the one w/ FALSE, a different gene is overlaped. The action with assign==TRUE still retains the one w/ overlapped with the other gene
	stopifnot(!duplicated(fo$queryHits))
	stopifnot(nrow(fo)==length(genes.tored))
	dt <- makeDT(genes.tored)
	dt[fo$queryHits,tu:= tu[fo$subjectHits]$tu] #bring tu id
	genes.tu <- c(genes.singles,makeGRanges(dt,strand=TRUE))

	stopifnot(length(unique(genes.tu$tu))==length(tu))
	stopifnot(length(genes.tu)==length(genes))
	stopifnot(length(unique(genes$gene.id))==length(unique(genes.tu$gene.id)))
	seqlevels(genes.tu) <- chrs
	genes.tu <- genes.tu[order(genes.tu)]

	# Assign TUs as coding
	message("Detecting coding TUs")
	genes.tu$coding <- genes.tu$cdsStart!=genes.tu$cdsEnd
	dt <- makeDT(genes.tu)
	dt <- dt[,list(coding=any(coding)),by="tu"]
	stopifnot(nrow(dt)==length(tu))
	tu$coding <- "new"
	tu[match(dt$tu,tu$tu)]$coding <- dt$coding

	# Get 3' end flanks (i.e., 5k bp segment right after TSE)
	message("Generating TU flanks")
	tu.flank <- flank(tu,width=flank,start=FALSE,both=FALSE,ignore.strand=FALSE)

	# Get 3' UTRs
	message("Clustering 3' UTRs for each coding TU")
	coding <- genes.tu[genes.tu$coding==TRUE]

	# If on (+), then want range between cdsEnd and txEnd
	# If on (-), then want range between txStart and cdsStart
	coding.p <- coding[strand(coding)=="+"]
	utr3.p <- GRanges(seqnames(coding.p),IRanges(coding.p$cdsEnd,end(coding.p)),strand=strand(coding.p),tu=coding.p$tu,gene.id=coding.p$gene.id,name=coding.p$name)
	coding.m <- coding[strand(coding)=="-"]
	utr3.m <- GRanges(seqnames(coding.m),IRanges(start(coding.m),coding.m$cdsStart),strand=strand(coding.m),tu=coding.m$tu,gene.id=coding.m$gene.id,name=coding.m$name)
	stopifnot((length(coding.p)+length(coding.m))==length(coding))

	redme <- function(myutrs)
	{			
		df <- as.data.frame(table(myutrs$tu))
		myutrs.singles <- myutrs[myutrs$tu %in% df[df$Freq==1,]$Var1]
		myutrs.tored <- myutrs[myutrs$tu %in% df[df$Freq>1,]$Var1]
		spl <- split(myutrs.tored,myutrs.tored$tu)
		spl.red <- lapply(spl,reduce)
		for(i in 1:length(spl.red))
		{
			spl.red[[i]]$tu <- names(spl)[i]
		}
		myutrs.red <- do.call(c,unname(spl.red))
		myutrs.red$gene.id <- myutrs[match(myutrs.red$tu,myutrs$tu)]$gene.id
		myutrs.red$name <- myutrs[match(myutrs.red$tu,myutrs$tu)]$name
		c(myutrs.singles,myutrs.red)
	}
	utr3.red.m <- redme(utr3.m)
	utr3.red.p <- redme(utr3.p)
	utr3 <- c(utr3.red.p,utr3.red.m)
	utr3 <- utr3[order(utr3)]
	utr3 <- utr3[width(utr3)>0]
	utr3$utr <- paste0("UTR",1:length(utr3))

	# Internal
	message("Computing coding gene internal regions via setdiff")	
	spl <- split(utr3,utr3$tu)
	coding <- tu[tu$coding==TRUE]
	# Remove coding genes that don't have 3' UTRs
	coding <- coding[coding$tu %in% names(spl)]
	diff <- mclapply(coding$tu,function(x) setdiff(coding[coding$tu==x],spl[[x]]),mc.cores=ncore)
	for(i in 1:length(diff))
	{
		diff[[i]]$tu <- coding$tu[i]
	}
	
	internal <- do.call(c,diff)	
	internal$gene.id <- tu[match(internal$tu,tu$tu)]$gene.id
	internal$name <- tu[match(internal$tu,tu$tu)]$name

	# No-UTR gene bodies
	coding <- tu[tu$coding==TRUE]
	noutr <- coding[!(coding$tu %in% names(spl))]

	# Between annotated 3' UTRs
	message("Computing between multi-UTR")
	dt <- makeDT(utr3)
	dt <- dt[,list(chr=chr[1],start=min(start),end=max(end),strand=strand[1],gene.id=gene.id[1],name=name[1],nUtr=length(utr)),by="tu"]
	multis <- dt[nUtr>1,]
	# Now have the maximal range for all the multis
	# If we setdiff out the real 3', then we have just the between ranges
	union <- makeGRanges(multis,strand=T)
	spl <- split(utr3,utr3$tu)
	diff <- mclapply(union$tu,function(x) setdiff(union[union$tu==x],spl[[x]]),mc.cores=ncore)
	for(i in 1:length(diff))
	{
		diff[[i]]$tu <- union$tu[i]
	}
	
	btw <- do.call(c,diff)	
	btw$gene.id <- tu[match(btw$tu,tu$tu)]$gene.id
	btw$name <- tu[match(btw$tu,tu$tu)]$name

	# Now remove the betweens from the internals
	internal <- internal[!(internal %in% btw)]
	
	ret <- list(tu=tu,genes.tu=genes.tu,flank=tu.flank,utr3=utr3,between=btw,internal=internal,noutr3=noutr)
	#writeBEDFromGRanges(ret$tu,file="output/tu_tu.bed",name="tu")
	#writeBEDFromGRanges(ret$flank,file="output/tu_flank.bed",name="tu")
	#writeBEDFromGRanges(ret$utr3,file="output/tu_utr3.bed",name="tu")
	#writeBEDFromGRanges(ret$between,file="output/tu_between.bed",name="tu")
	#writeBEDFromGRanges(ret$internal,file="output/tu_internal.bed",name="tu")
	#writeBEDFromGRanges(ret$noutr3,file="output/tu_noutr3.bed",name="tu")
	return(ret)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
joinTus <- function(pr,rg)
{
	prs <- pr$pr
  
  # tu_red <- reduce(rg$tu)
  # fodt <- data.table(data.frame(findOverlaps(rg$tu,tu_red)))
  # fodt[,cnt:=.N,by=subjectHits]
  # rg$tu[fodt[cnt>1,queryHits],]
  # tu_red[fodt[cnt>1,subjectHits],]
  # tsv <- '/media/tommy/cache/exp_out/ating/polyASeqs/04_AnnoApa/output/overlapped_tu.tsv'
  # fwrite(file=tsv,data.table(data.frame(rg$tu[fodt[cnt>1,queryHits],])),sep="\t",quote=FALSE,row.names=FALSE)
  
	fo_tu <- data.table(as.data.frame(findOverlaps(prs,rg$tu)))
	fo_flank <- data.table(as.data.frame(findOverlaps(prs,rg$flank)))

	fo_tu$set <- "tu"
	fo_flank$set <- "flank"
	fo_tu$tu <- rg$tu[fo_tu$subjectHits]$tu
	fo_tu$coding <- rg$tu[fo_tu$subjectHits]$coding
	fo_flank$tu <- rg$flank[fo_flank$subjectHits]$tu
	fo_flank$coding <- rg$flank[fo_flank$subjectHits]$coding

	# Remove cases where same PR links to same TU and flank of that TU
	fo_flank$key <- paste0(fo_flank$queryHits,"+",fo_flank$tu)
	fo_tu$key <- paste0(fo_tu$queryHits,"+",fo_tu$tu)
	fo_flank <- fo_flank[!(fo_flank$key %in% fo_tu$key),]

	ov <- rbind(fo_tu,fo_flank)

	# Join table that links each PR to a TU/flank
	join <- data.table(pr=prs[ov$queryHits]$pr,tu=ov$tu,type=ov$set,coding=ov$coding)
	
	# Count table for listing the contingencies
	agg <- join[,list(nTu=length(tu[type=="tu"]),nFlank=length(tu[type=="flank"])),by="pr"]

	# Reduced assignment table that assigns the uniques
	join <- join[,list(tu=tu,type=type,coding=coding,unique_pr=length(tu)==1,over_tus=toString(tu[type=="tu"]),flank_tus=toString(tu[type=="flank"])),by="pr"]
	join <- join[,list(pr=pr,type=type,coding=coding,unique_pr=unique_pr,unique_tu=all(unique_pr),over_tus=over_tus,flank_tus=flank_tus),by="tu"]

	#link <- join[,list(tus=toString(tu),uniq_pr=length(tu)==1),by="pr"]
	#link2 <- join[,list(tu=tu,uniq_pr=length(tu)==1),by="pr"]
	
	# Flag TUs where all PRs unique
	#utu <- link2[,list(uniq_tu=all(uniq_pr)),by="tu"]

	summary <- data.frame(group="intergenic",n=length(prs)-nrow(agg),stringsAsFactors=FALSE)
	summary <- rbind(summary,data.frame(group="unique_tu",n=sum((agg$nTu==1)&(agg$nFlank==0))))
	summary <- rbind(summary,data.frame(group="unique_flank",n=sum((agg$nTu==0)&(agg$nFlank==1))))
	summary <- rbind(summary,data.frame(group="multi_tu",n=sum((agg$nTu>1)&(agg$nFlank==0))))
	summary <- rbind(summary,data.frame(group="multi_flank",n=sum((agg$nTu==0)&(agg$nFlank>1))))
	summary <- rbind(summary,data.frame(group="unique_tu_multi_flank",n=sum((agg$nTu==1)&(agg$nFlank>0))))
	summary <- rbind(summary,data.frame(group="multi_tu_multi_flank",n=sum((agg$nTu>1)&(agg$nFlank>0))))
	summary$frac <- summary$n/sum(summary$n)
	stopifnot(sum(summary$n)==length(prs))
	
	ret <- list(prs=prs,join=join,agg=agg,summary=summary)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test for APA with DEXSeq
testApa <- function(samp,pr,jtu,a,b,cov=10,adjust.var=NULL,ncpu=1,prs_gt=1)
{
  #browser()
	message("Testing APA for ",a," vs ",b," using DEXSeq")
	stopifnot(a %in% samp$group)
	stopifnot(b %in% samp$group)

	message("Starting with ",length(unique(jtu$join$pr))," pA sites accross ",length(unique(jtu$join$tu))," TUs")

	# Require mean >= 10 in at least a or b before moving to testing
	matmean <- pr$means$raw
	pr.expr <- rownames(matmean)[(matmean[,a]>=cov)|(matmean[,b]>=cov)]
	pr.a <- rownames(matmean)[(matmean[,a]>=cov)]
	pr.b <- rownames(matmean)[(matmean[,b]>=cov)]

	# Grab unique TUs only (all PR of TU are unique)

	join.a <- jtu$join
	join.a <- join.a[(pr %in% pr.a)&(unique_tu==TRUE),]
	length(unique(join.a$tu))
	join.b <- jtu$join
	join.b <- join.b[(pr %in% pr.b)&(unique_tu==TRUE),]
	length(unique(join.b$tu))
	
	# browser()
	
	join.want <- jtu$join[((pr %in% pr.expr) & unique_tu==TRUE),]
	stopifnot(!str_detect(join.want$over_tus,","))
	stopifnot(!str_detect(join.want$flank_tus,","))
	stopifnot(all(join.want$unique_pr))
	stopifnot(length(unique(join.want$pr))==nrow(join.want))
	stopifnot(class(join.want$pr)=="character")
	stopifnot(join.want$pr %in% rownames(pr$counts$raw))
	join.want <- join.want[,list(pr=pr,
	                             type=type,
	                             coding=coding,
	                             unique_pr=unique_pr,
	                             unique_tu=unique_tu,
	                             over_tus=over_tus,
	                             flank_tus=flank_tus,
	                             num_pr=length(pr)),
	                       by="tu"]
	
	message(sprintf("Subsetting to TUs with > %d pA site",prs_gt))
	message("After subsetting to TUs with all PRs unique: ",length(unique(join.want$pr))," sites accross ",length(unique(join.want$tu))," TUs")
	join.want <- join.want[num_pr > prs_gt,] # ignore any tu where at most one pr appears only
	message(sprintf("After subsetting to > %d pA site TUs: ",prs_gt),length(unique(join.want$pr))," sites accross ",length(unique(join.want$tu))," TUs")

	# Get counts for them
	mat <- pr$counts$raw
	mat <- mat[join.want$pr,]
	stopifnot(rownames(mat)==join.want$pr)
	rownames(mat) <- paste(join.want$tu,join.want$pr,sep=":")

	#browser()
	
	# Build DEXSeq object
	message("Subsetting count matrix for groups 'a' and 'b' only")
	cd <- mat[,samp$group %in% c(a,b)]
	sd <- samp[samp$group %in% c(a,b),]
	stopifnot(samp$sample==colnames(mat))
	fr <- pr$pr[pr$pr$pr %in% join.want$pr]
	fr <- fr[match(join.want$pr,fr$pr)]
	stopifnot(fr$pr==join.want$pr)
	sd$group <- factor(sd$group,levels=c(a,b))
	rownames(cd) <- NULL
	sd.use <- data.frame(group=sd$group)
	sd.use$group <- factor(sd.use$group,ordered=F)

	# browser()
	
	message("Building DEXSeqDataSet object")
	if(is.null(adjust.var))
	{
		design <- "~ sample + exon + group:exon"
	} else
	{
		design <- paste0("~ sample + exon + adjust:exon + group:exon")
		sd.use$adjust <- factor(with(sd,get(adjust.var)))
	}
	message("Design Formula: ",design)

	exon_num <- join.want[,list(paste0("E",1:length(pr)),pr=pr),by="tu"]
	stopifnot(exon_num$tu==join.want$tu)
	stopifnot(exon_num$pr==join.want$pr)
	join.want$exon_num <- exon_num$V1
	join.want$key <- paste(join.want$tu,join.want$exon_num,sep=":")

	# browser()
	
	dxd <- DEXSeqDataSet(countData=cd, sampleData=sd.use, design=as.formula(design), featureID=join.want$exon_num, groupID=join.want$tu, featureRanges=fr)

	# Run tests
	message("Estimating size factors")
	dxd <- estimateSizeFactors(dxd)
	message("Estimating dispersions")
	dxd <- estimateDispersions(dxd)
	#pdf(file="output/dexseq_plotDispEsts.pdf")
	#plotDispEsts(dxd)
	#dev.off()
	message("Performing diffTesting")
	if(is.null(adjust.var))
	{
		dxd <- testForDEU(dxd)
	} else
	{
		message("Tesing w/ adjustment for variable: ", adjust.var)
		dxd <- testForDEU(dxd,fullModel=as.formula("~ sample + exon + adjust:exon + group:exon"),reducedModel=as.formula("~ sample + exon + adjust:exon"))
	}

	dxd <- estimateExonFoldChanges(dxd,fitExpToVar="group",BPPARAM=MulticoreParam(workers=ncpu))
	
	# Pull APA results table
	message("Extracting DEXSeqResults")
	dxr <- DEXSeqResults(dxd)

	# Fetch normalized counts
	counts.norm <- counts(dxd, normalized=TRUE)
	counts.norm <- counts.norm[,1:nrow(sd)]
	colnames(counts.norm) <- sd$sample
	stopifnot(rownames(counts.norm)==join.want$key)
	rownames(counts.norm) <- join.want$pr

	# Start building nice results list
	dt <- data.table(as.data.frame(dxr))
	out <- with(dt,data.table(tu=groupID,pr=join.want$pr,p=pvalue,padj=padj,dexl2fc_b_over_a=get(paste0("log2fold_",b,"_",a)),chr=genomicData.seqnames,start=genomicData.start,end=genomicData.end,strand=genomicData.strand))
	out$gene_name <- red_ens$tu[match(out$tu,red_ens$tu$tu)]$name
	#stopifnot(paste0(out$tu,":",out$pr)==rownames(counts.norm))

	#browser()
	
	# Make usage matrices
	c2 <- data.table(cd)
	c2$pr <- out$pr
	c2$tu <- out$tu
	cm <- melt(c2,id.vars=c("pr","tu"))
	stopifnot(!is.na(cm$value))
	# Within each TU and each sample, get each percent as percent of total
	#m2 <- cm[,list(pr=pr,raw_count=value,use_frac=ifelse(sum(value)==0,0,value/sum(value))),by=c("tu","variable")]
	m2 <- cm[,list(pr=pr,raw_count=value,use_frac=value/sum(value),tu_zero=sum(value)==0),by=c("tu","variable")]
	stopifnot(nrow(m2[is.na(use_frac),])==sum(m2$tu_zero))
	m2[is.na(use_frac),use_frac:=0]

	stopifnot(!is.na(m2$use_frac))
	m2$group <- samp[match(m2$variable,samp$sample),]$group
	stopifnot(!is.na(m2$group))

	# Groupwise usage means
	m3 <- m2[,list(tu=tu[1],mean_a=mean(use_frac[group==a]),mean_b=mean(use_frac[group==b])),by=c("pr")]
	m3$b_minus_a <- m3$mean_b-m3$mean_a
	m3$l2_b_over_a <- log2((m3$mean_b+0.0001)/(m3$mean_a+0.0001))
	stopifnot(out$tu==m3$tu)
	stopifnot(out$pr==m3$pr)
	setnames(m3,paste0(colnames(m3),"_frac"))

	cas <- cast(m2,formula="pr+tu~variable",value="use_frac")
	caskey <- paste(cas$tu,cas$pr,sep=":")
	outkey <- paste(out$tu,out$pr,sep=":")
	m <- match(outkey,caskey)
	stopifnot(!is.na(m))
	cas <- cas[m,]
	use.frac <- cas[,c(-1,-2)]
	stopifnot(nrow(use.frac)==nrow(out))
	rownames(use.frac) <- join.want$pr

	#browser()
	
	#myraw <- cd
	#colnames(myraw) <- paste0(colnames(cd),"_raw")
	myuse <- use.frac
	colnames(myuse) <- paste0(colnames(use.frac),"_frac")	
	mynorm <- counts.norm
	colnames(mynorm) <- paste0(colnames(mynorm),"_norm")

	nmeans_a <- rowMeans(counts.norm[,as.character(samp[samp$group==a,]$sample)])
	nmeans_b <- rowMeans(counts.norm[,as.character(samp[samp$group==b,]$sample)])
	dt <- data.table(mean_a=nmeans_a,mean_b=nmeans_b,b_minus_a=nmeans_b-nmeans_a,l2_b_over_a=log2((nmeans_b+0.0001)/(nmeans_a+0.0001)))
	setnames(dt,paste0(colnames(dt),"_norm"))

	res <- cbind(out,mynorm,dt,myuse,m3[,c("mean_a_frac","mean_b_frac","b_minus_a_frac","l2_b_over_a_frac"),with=F])

	rownames(cd) <- res$pr

	#browser()
	
	# Return
	list(dxd=dxd,dxr=dxr,res=res,counts.raw=cd,counts.norm=counts.norm,use.frac=use.frac,a=a,b=b,samp.sub=sd)
}
