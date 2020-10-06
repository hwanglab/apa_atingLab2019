# Collation of stats of read counts and where they went and why

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

load("../03_CallApa/output/prepro.rd")
load("../03_CallApa/output/apa.sig.rd")

stat <- fread("../02_ProcessBams/output/AlignmentFullStats.csv")

stat <- stat[,1:6,with=F]

# Add in the number of reads before A-chop (total that de-barcoded)
index <- fread("../01_ProcessSeq/input/fastq_index.csv")
index$file <- str_replace(str_replace(index$file,"fastq/",""),"\\.fastq\\.gz","")
index$key <- paste0(index$file,"_",index$barcode)
chop <- fread("../01_ProcessSeq/output/ChopATailsResults.csv")
stopifnot(index$key %in% chop$file)

sub <- chop[match(index$key,chop$file),]
sub$total <- rowSums(sub[,-1,with=F])

sub$sample <- index[match(sub$file,index$key),]$sample

stat$TotalDeBarcoded <- sub[match(stat$sample,sub$sample),]$total

stat <- stat[,c(1,ncol(stat),2:(ncol(stat)-1)),with=F]

# Add in number of reads left after removing genomic A
ga <- fread("../02_ProcessBams/output/GenomicA_Fractions.csv")
stopifnot(ga$sample==stat$sample)
stat$NoGenomicA <- ga$noGenomicA

# Save file
write.csv(stat,row.names=F,file="output/stat_table.csv")

