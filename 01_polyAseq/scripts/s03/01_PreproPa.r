library(handy)
library(argparse)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))
source(file.path(Sys.getenv('PADD_GIT'),'readGtf.r'))

if (T) { #debug
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))

  parser <- ArgumentParser(description='prepropa')
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                      dest="ncpu", default = 1,
                      help="number of cpus to utilize [1]")
  
  parser$add_argument("-i", "--input_file", type="character", required=TRUE,
                      dest="input_file",
                      help="sample_sheet file")
  
  parser$add_argument("-B", "--bam_dir", type="character", required=TRUE,
                      dest="bam_dir",
                      help="bam_dir")
  
  parser$add_argument("-f", "--ftp_serverl_url", type="character", required=TRUE,
                      dest="ftp_server_url",
                      help="ftp_server_url")
  
  parser$add_argument("-b", "--ensembl_biotype_table", type="character", required=TRUE,
                      dest="ensembl_biotype_table",
                      help="ensembl_biotype_table")
  
  parser$add_argument("-g", "--ensembl_gtf_file", type="character", required=TRUE,
                      dest="ensembl_gtf_file",
                      help="ensembl_gtf_file")
  
  parser$add_argument("-T", "--cache_dir", type="character", required=TRUE,
                      dest="cache_dir",
                      help="cache directory")
  
  parser$add_argument("-t", "--track_dir", type="character", required=TRUE,
                      dest="track_dir",
                      help="track directory")
  
  parser$add_argument("-o", "--output_file", type="character", required=TRUE,
                      dest="output_file",
                      help="output file to store workspace R variables")
  args <- parser$parse_args()
} else {
  # ----------------
  source('../../paddle-git/readGtf.r')
  args = DataFrame(input_file="",output_file="",track_dir="",bam_dir="",ftp_server_url="",ensembl_biotype_table="",ensembl_gtf_file="",cache_dir="")
  
  args$input_file="~/projects/apa_atingLab2019/01_polyAseq/01_wkd/fastq_index.csv"
  args$output_file="~/projects/apa_atingLab2019/01_polyAseq/01_wkd/out/03_CallApa/output/prepropa.rd"
  args$track_dir="~/projects/apa_atingLab2019/01_polyAseq/01_wkd/out/03_CallApa/output/tracks"
  args$bam_dir="/home/hongc2/projects/apa/paper_submission/polyAseq/01_wkd/out/02_ProcessBams/output/rm_int_a_7of10or6"
  args$ftp_server_utl="ftp://account:password@server_url/bw_ucsc_tracks/"
  args$ensembl_biotype_table="~/projects/apa_atingLab2019/resource/ensembl/ensembl_gene_biotypes.csv"
  args$ensembl_gtf_file="~/projects/apa_atingLab2019/resource/ensembl/Homo_sapiens.GRCh37.87.gtf"

  args$cache_dir="/tmp"
  args$ncpu = 6
}
# ----------------
# All data processing up until APA calling
# --------------------------------------------------------------------
# Preprocessing

# Settings
ncore <- args$ncpu
chrs <- handy::chrs()

# Samples
# Using the version that was already done, rather than the re-run version for now
# Separate sample sheets for each one
#samp <- readSampleInfo("input/samples-all-paper-run.csv")

#read sample sheet file and read bam dir, then, create a samples table

sample_dt <- appendBamFilePath(args$input_file,args$bam_dir,bam_suffix=".02_mapq20_DePcrDup_noGenomicA.bam")

exp_sheet_fn_tmp <- sprintf("%s.tmp",args$output_file)
write.table(sample_dt,row.names=FALSE,quote=FALSE,sep=",",col.names=TRUE,file=exp_sheet_fn_tmp)

samp <- readSampleInfo(exp_sheet_fn_tmp)

# Load reads from BAMs
reads <- getReads(samp=samp,chrs=chrs,ncore=ncore)

# Save bigWigs of coverages
saveCov(samp=samp,reads=reads,path=args$track_dir,bigwig=TRUE,server=args$ftp_server_url,ncore=ncore)

# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Make processing regions
dist <- 10
cov <- 10
pr <- makeProcReg(samp,reads,dist=dist,cov=cov,ncore=ncore)
handy::nsummary(pr)
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Grab gene database
library(goldmine)
genes <- getGenes("ensembl",genome="hg19",cachedir=args$cache_dir,sync=FALSE)
genes <- genes[chr %in% chrs,]

tab <- fread(args$ensembl_biotype_table)

bio.want <- tab[tab$group %in% c("Protein_coding","Long_noncoding"),]$Var1

gtf <- readGtf(args$ensembl_gtf_file,feature="gene")
stopifnot(genes$gene.id %in% gtf$gene_id)
stopifnot(!duplicated(gtf$gene_id))
genes$biotype <- gtf[match(genes$gene.id,gtf$gene_id),]$gene_biotype
stopifnot(!is.na(genes$biotype))
genes <- genes[biotype %in% bio.want,]
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Reduce genes to region sets for joining/annotation
red_ens <- reduceGenes(genes=genes,chrs=chrs,flank=5000,ncore=ncore)

# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Join PRs to TUs
jtu <- joinTus(pr=pr,rg=red_ens)
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Save the "preprocessing"
save(jtu,pr,red_ens,genes,samp,ncore,reads,file=args$output_file,compress=T)
# --------------------------------------------------------------------

if (file.exists(exp_sheet_fn_tmp)){
  unlink(exp_sheet_fn_tmp)
}
