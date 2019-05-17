library(data.table)
library(argparse)

if (T) {
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  
  parser <- ArgumentParser(description='deseq2')
  
  parser$add_argument("-c", "--ctrl_bam", type="character", required=TRUE,
                      dest="ctrl_bam",
                      help="control sample bam file path")
  
  parser$add_argument("-e", "--expr_bam", type="character", required=TRUE,
                      dest="expr_bam",
                      help="control sample bam file path")

  parser$add_argument("-t", "--comp_tag", type="character", required=FALSE,
                      dest="comp_tag", default="HCT116_DKO",
                      help="comparison tag")
  
  parser$add_argument("-o", "--out_dir", type="character", required=TRUE,
                      dest="out_dir",
                      help="output directory")
  
  args <- parser$parse_args()
} else {
  args <- NA
}

library(handy)
library(openxlsx)
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))

run_macs2_deseq2 <- function(comp12,ctrl_bam,expr_bam,out_dir) {
  prog <- file.path(Sys.getenv('R_UTIL'),'macs2_deseq2.r')
  cmd <- sprintf("Rscript %s -c %s -t %s -p broad -m merged -q 0 -o %s",
								 prog,
                 ctrl_bam,
								 expr_bam,
                 out_dir)
  message(cmd)
  system(cmd) #debug

  return(sprintf("%s/deseq2.tsv",out_dir))
}

extract_diff <-function(deseq2_report_file,min_raw_count=3,min_abs_fc=1.){

  deseq2_report_file_filt <- sprintf('%s.rc%d_fc%d.tsv',deseq2_report_file,min_raw_count,min_abs_fc)

  mbdif <- fread(deseq2_report_file)
  mbdif <- mbdif[(raw_control>min_raw_count|raw_experiment>min_raw_count),]
  
  mbdif[,'Event':='TBA']
  mbdif[, 'log2FoldChange_raw' := mapply(log2fc,mbdif$raw_control,mbdif$raw_experiment)]
  mbdif[log2FoldChange_raw >= min_abs_fc, Event:='Up']
  mbdif[log2FoldChange_raw <= -min_abs_fc, Event:='Down']
  
  mbdif[,c("chrom","start","end"):=transpose(setDT(strsplit(mbdif$region, "[:-]+")))]
  mbdif[, start := as.numeric(start)]
  mbdif[, end := as.numeric(end)]
  mbdif[, peak_width := (end - start)]
  
  mbdif[,c('region','baseMean','stat','padj'):=NULL]
  chrs <- handy::chrs()
  mbdif <- mbdif[((chrom %in% chrs) & (Event!='TBA')),]
  
  setorder(mbdif,Event,chrom,start,end)
  
  setcolorder(mbdif,c('Event','chrom','start','end','peak_width','raw_control','raw_experiment','log2FoldChange_raw','norm_control','norm_experiment','log2FoldChange','lfcSE','pvalue'))

  
  fwrite(mbdif,file=deseq2_report_file_filt,sep='\t',col.names = T,quote = F) #debug

  return(deseq2_report_file_filt)
  
}


deseq2_report <- run_macs2_deseq2(args$ctrl_bam,args$expr_bam,args$out_dir)

deseq2_reports_filt <- extract_diff(deseq2_report,min_raw_count=3,min_abs_fc=1.)

wb <- createWorkbook()

message(sprintf('adding %s',deseq2_reports_filt))
dt <- fread(deseq2_reports_filt,header=T)

sheet_name <- args$comp_tag
addWorksheet(wb,sheet_name)
writeData(wb, sheet=sheet_name,x=dt)
message(sprintf('[%s] added',sheet_name))

saveWorkbook(wb, file.path(args$out_dir,"Supplementary_Table_S4_mbd_bt2.xlsx"), overwrite = TRUE)
