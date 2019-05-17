library(data.table)

b37_chroms <- function () 
{
  paste0(c(1:22, "X", "Y", "MT"))
}

readGtf <- function(gtf_file,feature="exon")
{
	message("Reading GTF: ",gtf_file)
	if (summary(file(gtf_file))$class == "gzfile") {
      dt <- fread(cmd=sprintf("zgrep -v '^#!' %s",gtf_file))
  } else {
	    dt <- fread(cmd=sprintf("grep -v '^#!' %s",gtf_file))
  }

	dt <- dt[V3==feature,]

	message("Vectorized attribute parsing")
	metadat <- str_split(dt$V9,";")
	metadat <- lapply(metadat,function(x) x[1:(length(x)-1)])
	message("Extract column data")
	metadat <- lapply(metadat,function(x) str_replace_all(str_replace(x,"[^\"]*\"",""),"\"",""))
	
	message("Extract column names")
	metacols <- dt$V9
	metacols <- str_replace_all(metacols,"\"[^;]*\";","")
	metacols <- str_split(metacols," ")
	metacols <- lapply(metacols,function(x) x[x!=""])

	metaline <- lapply(1:length(metadat),function(x) rep(x,length(metadat[[x]])))

	stopifnot(length(metacols)==length(metadat))
	stopifnot(length(metacols)==length(metaline))
	stopifnot(sapply(metacols,length)==sapply(metadat,length))
	stopifnot(sapply(metacols,length)==sapply(metaline,length))

	message("Concat extracted data")
	attr <- data.table(line=do.call(c,metaline),name=do.call(c,metacols),dat=do.call(c,metadat))

	message("Casting attribute matrix")
	# all column names should be less than the nrow
	check <- as.data.frame(table(attr$name))
	check$check <- check$Freq<=nrow(dt)
	# tag seems to be the only offender
	# makes sense to concat these
	attr1 <- attr[name!="tag",]
	attr2 <- attr[name=="tag",]
	attr2 <- attr2[,list(name="tag",dat=toString(dat)),by="line"]
	attr <- rbind(attr1,attr2)
	check <- as.data.frame(table(attr$name))
	check$check <- check$Freq<=nrow(dt)
	stopifnot(check$check)

	cas <- dcast(attr,formula="line~name",value.var="dat")
	stopifnot(nrow(cas)==nrow(dt))
	stopifnot(cas$line==seq(1,nrow(dt)))
  
	
	out <- data.table(chr=dt$V1,source=dt$V2,feature=dt$V3,start=dt$V4,end=dt$V5,score=dt$V6,strand=dt$V7,frame=dt$V8)
	out <- cbind(out,cas)
}

convert_ucsc_asm_head <- function(gtf_dt){
  b37_chrs <- b37_chroms()
  
  gtf_dt$chr <- mapply(function(x){
    if (x %in% b37_chrs) {
      if (x == 'MT'){
        'chrM'
      } else {
        paste0('chr',x)
      }
    } else {
      NA
    }
  },
  gtf_dt$chr)
  
  return(gtf_dt)
  
}