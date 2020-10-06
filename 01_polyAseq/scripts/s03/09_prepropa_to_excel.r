#objectives: take PR1_vs_PR2.csv file from 01_Prepropa.r and convert the csv file to excel file for supplementary table 2

library(data.table)
library(openxlsx)
library(argparse)

args_tmp <- commandArgs(trailingOnly = F)
scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))

if (T) {

  parser <- ArgumentParser(description='prepropa')
  
  parser$add_argument("-i", "--comp_dir", type="character", required=TRUE,
                      dest="comp_dir",
                      help="input directory")
  
  parser$add_argument("-c", "--control_tag", type="character", required=FALSE, 
                      dest="control_tag", default="HCT116",
                      help="control group tag [HCT116]")
  
  parser$add_argument("-e", "--exp_tag", type="character", required=FALSE, 
                      dest="exp_tag",default="DKO",
                      help="experimental group tag [DKO]")

  args <- parser$parse_args()
} else {
  args <- data.table(comp_dir='../../01_wkd/out/03_CallApa',
                     control_tag='HCT116',
                     exp_tag='DKO')
}

prepropa_to_excel <- function(sample_tags,comp_dir) {
  
  # browser()
  ctrl_tag <- sample_tags[1]
  expr_tag <- sample_tags[2]
  
  polya_dex_file <- sprintf('%s/%s_vs_%s_APA.csv',comp_dir,ctrl_tag,expr_tag)
  
  message(polya_dex_file)
  stopifnot(file.exists(polya_dex_file))
  
  out_file <- sprintf('%s/%s_vs_%s_APA.xlsx',comp_dir,ctrl_tag,expr_tag)
  
  dt <- fread(polya_dex_file,header = T,sep = ",")
  cn <- colnames(dt)
  setnames(dt,"chr","pr_chr")
  setnames(dt,"start","pr_start")
  setnames(dt,"end","pr_end")
  setnames(dt,"strand","pr_strand")
  setnames(dt,"mean_a_norm",sprintf("mean_%s_norm",ctrl_tag))
  setnames(dt,"mean_b_norm",sprintf("mean_%s_norm",expr_tag))
  
  setnames(dt,"b_minus_a_norm",sprintf("%s_minus_%s_norm",expr_tag,ctrl_tag))
  setnames(dt,"l2_b_over_a_norm",sprintf("l2_%s_over_%s_norm",expr_tag,ctrl_tag))
  
  setnames(dt,"mean_a_frac",sprintf("mean_%s_frac",ctrl_tag))
  setnames(dt,"mean_b_frac",sprintf("mean_%s_frac",expr_tag))
  
  setnames(dt,"b_minus_a_frac",sprintf("%s_minus_%s_frac",expr_tag,ctrl_tag))
  setnames(dt,"l2_b_over_a_frac",sprintf("l2_%s_over_%s_frac",expr_tag,ctrl_tag))
  
  setnames(dt,"batchadj_padj","dex_padj")
  setnames(dt,"batchadj_delta_sig","delta_sig")
  setnames(dt,"batchadj_perfc_sig","fracfc_sig")
  setnames(dt,"int_sig","both_sig")

  dt[,c('dexl2fc_b_over_a','noadj_p','noadj_padj','batchadj_p','ttest_p','ttest_padj','edger_exon_p','edger_exon_padj','edger_ftest_p','edger_ftest_padj','edger_simes_p','edger_simes_padj'):=NULL]
  
  setcolorder(dt, c('gene_name',
                    'tu',
                    'pr',
                    'pr_chr',
                    'pr_start',
                    'pr_end',
                    'pr_strand',
                     sprintf('%s_1_norm',ctrl_tag),
                     sprintf('%s_2_norm',ctrl_tag),
                     sprintf('%s_3_norm',ctrl_tag),
                     sprintf('%s_4_norm',ctrl_tag),
                     sprintf('%s_1_norm',expr_tag),
                     sprintf('%s_2_norm',expr_tag),
                     sprintf('%s_3_norm',expr_tag),
                     sprintf('%s_4_norm',expr_tag),
                     sprintf('mean_%s_norm',ctrl_tag),
                     sprintf('mean_%s_norm',expr_tag),
                    sprintf('%s_minus_%s_norm',expr_tag,ctrl_tag),
                  sprintf('l2_%s_over_%s_norm',expr_tag,ctrl_tag),
                       sprintf('%s_1_frac',ctrl_tag),
                       sprintf('%s_2_frac',ctrl_tag),
                       sprintf('%s_3_frac',ctrl_tag),
                       sprintf('%s_4_frac',ctrl_tag),
                       sprintf('%s_1_frac',expr_tag),
                       sprintf('%s_2_frac',expr_tag),
                       sprintf('%s_3_frac',expr_tag),
                       sprintf('%s_4_frac',expr_tag),
                       sprintf('mean_%s_frac',ctrl_tag),
                       sprintf('mean_%s_frac',expr_tag),
                       sprintf('%s_minus_%s_frac',expr_tag,ctrl_tag),
                       sprintf('l2_%s_over_%s_frac',expr_tag,ctrl_tag),
                       'dex_padj',
                       'delta_sig',
                       'fracfc_sig',
                       'both_sig'))
  
  write.xlsx(dt, file = out_file, colNames = TRUE, quote=F, rowNames = F)
  message(sprintf('wrote [%s]',out_file))
  return(out_file)
}

if (F) {
  comp_groups <- list(c('HCT116','DKO'),
                      c('PrEC','LNCaP'),
                      c('PrEC','DU145'),
                      c('LNCaP','DU145'))
} else {
  xlsx_file <- prepropa_to_excel(c(args$control_tag,args$exp_tag),
                                 args$comp_dir)
  message(sprintf('Generated the output file:[%s]',xlsx_file))
}
