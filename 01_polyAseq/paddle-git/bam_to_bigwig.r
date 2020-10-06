library(handy)
library(argparse)

source(file.path(Sys.getenv('R_UTIL_APA'),'lib_apps.R'))
source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))
debug2 <- F

if (!debug2){
parser <- ArgumentParser(description='bam_to_bigwig')

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                    dest="ncpu", default = 1,
                    help="number of cpus to utilize [1]")

parser$add_argument("-i", "--input_file", type="character", required=TRUE,
                    dest="input_file",
                    help="sample_sheet file")

parser$add_argument("-d", "--bam_dir", type="character", required=TRUE,
                    dest="bam_dir",
                    help="bam_dir")

parser$add_argument("-B", "--bam_suffix", type="character", required=FALSE,
                    dest="bam_suffix", default= '*.bam',
                    help="bam_suffix")
parser$add_argument("-u", "--add_chr", type="integer", required=FALSE,
                    dest="add_chr", default = 0,
                    help="[0]:No, 1:Yes Add chr in the contig header")

parser$add_argument("-s", "--read_pair_orient", type="integer", required=FALSE,
                    dest="read_pair_orient", default = 0,
                    help="[0]:NA, 1:Yes(1st pair mapped strand), 2:Yes(2nd pair mapped strand)")

parser$add_argument("-f", "--ftp_serverl_url", type="character", required=FALSE,
                    dest="ftp_server_url", default = "https://public_domain",
                    help="ftp_server_url")

parser$add_argument("-T", "--cache_dir", type="character", required=FALSE,
                    dest="cache_dir", default = "/tmp",
                    help="cache directory")

parser$add_argument("-t", "--track_dir", type="character", required=TRUE,
                    dest="track_dir",
                    help="track directory")

args <- parser$parse_args()
} else {
  args <- data.table(input_file="/home/hongc2/projects/apa/mrna/samples-all-mrna.csv",
                     bam_dir="/home/hongc2/projects/apa/mrna",
                     bam_suffix="_chr.bam",
                     add_chr=0,
                     read_pair_orient=2,
                     ncpu=4,
                     track_dir="/home/hongc2/projects/apa/mrna/mrnaseq_track")
                     
}
# Settings
ncore <- args$ncpu
chrs <- handy::chrs()

# Samples
# Using the version that was already done, rather than the re-run version for now
# Separate sample sheets for each one
#samp <- readSampleInfo("input/samples-all-paper-run.csv")

#read sample sheet file and read bam dir, then, create a samples table
sample_dt <- appendBamFilePath2(args$input_file,args$bam_dir,bam_suffix=args$bam_suffix)
if (args$add_chr==1){
	bam2 <- lapply(sample_dt$bam,add_chr_contig_header)
	sample_dt$bam <- bam2
}

stbin <- "samtools"

min_mapq <- 1

cnt <- 1
for (target in sample_dt$group2){
  if (cnt > 1) {break}
  cnt <- cnt + 1
  sample_dt_tmp <- sample_dt[group2%in%target,]
  
  exp_sheet_fn_tmp <- sprintf("%s.tmp",args$input_file)
  fwrite(sample_dt_tmp,file=exp_sheet_fn_tmp,quote=F,sep=',',col.names = T,row.names = F)
  samp <- readSampleInfo(exp_sheet_fn_tmp)
  
  # Load reads from BAMs
  if (args$read_pair_orient>0) {
    
    if (args$read_pair_orient == 1) {
      message('use 1st pair mapping orientation')
      sam_flag_comb <- list('fwd'=c('-f0x41 -G0x51','0x91'),'rev'=c('0x51','-f0x81 -G0x91'))
    } else if (args$read_pair_orient==2) {
      message('use 2nd pair mapping orientation')
      sam_flag_comb <- list('fwd'=c('-f0x51','-f0xA1 -G0x91'),'rev'=c('-f0x61 -G0x51','-f0x91'))
    }
    
    orig_bam_fns <- samp$bam
    reads2 <- lapply(1:nrow(samp),function(j) {
      orig_bam_fn <- orig_bam_fns[[j]]
      
      adj_bam_fns <- sapply(1:length(sam_flag_comb),function(i){
        sflag <- unlist(sam_flag_comb[i])
        sflag_name <- names(sam_flag_comb[i])
        
        merged_bam_fn <- sprintf('%s.12_%s.bam',orig_bam_fn,sflag_name)
        if (!file.exists(merged_bam_fn)) {
          out_bam1_fn <- sprintf('%s.1_%s.bam',orig_bam_fn,sflag_name)
          cmd <- sprintf("%s view -@%d -b %s -q%d -o %s %s\n",stbin,ncore,sflag[1],min_mapq,out_bam1_fn,orig_bam_fn)
          message(cmd)
          system(cmd)
          out_bam2_fn <- sprintf('%s.2_%s.bam',orig_bam_fn,sflag_name)
          cmd <- sprintf("%s view -@%d -b %s -q%d -o %s %s\n",stbin,ncore,sflag[2],min_mapq,out_bam2_fn,orig_bam_fn)
          message(cmd)
          system(cmd)
  
          cmd <- sprintf("%s merge -@%d %s %s %s\n",stbin,ncore,merged_bam_fn,out_bam1_fn,out_bam2_fn)
          message(cmd)
          system(cmd)
        }
        
        merged_bai_fn <- sprintf('%s.12_%s.bam.bai',orig_bam_fn,sflag_name)
        if (!file.exists(merged_bai_fn)) {
          cmd <- sprintf("%s index %s",stbin,merged_bam_fn)
          system(cmd)
        }
        return(merged_bam_fn)
      })
      
      samp_temp <- data.table(sample=c('fwd','rev'),bam=adj_bam_fns)
      read <- getReads(samp=samp_temp,chrs=chrs,ncore=ncore)
      read <- lapply(read,unique)
      strand(read$fwd)<-'+'
      strand(read$rev)<-'-'
      read2 <- c(read$fwd,read$rev)
      read2 <- sortSeqlevels(read2)
      read2 <- sort(read2)
      return(read2)
    })
    names(reads2)<-samp$sample
    saveCov(samp=samp,reads=reads2,path=args$track_dir,bigwig=TRUE,server=args$ftp_server_url,ncore=ncore,track_suffix=target)
  } else {
    message('no strand mode ...')
  }
}
