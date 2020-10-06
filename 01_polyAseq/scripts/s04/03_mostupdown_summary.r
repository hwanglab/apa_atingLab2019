# Annotate the APA results with the MostUp/MostDown information
library(argparse)
library(data.table)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))
if (T) {
  parser <- ArgumentParser(description='summary_apa_table')
  debug <- T
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                      dest="ncpu", default = 1,
                      help="number of cpus to utilize [1]")
  
  parser$add_argument("-p", "--pr_table", type="character", required=TRUE,
                      dest="input_file",
                      help="PR xlsx table generated from s03/09_prepropa_to_excel.r")
  
  parser$add_argument("-a", "--apa_table", type="character", required=TRUE,
                      dest="apa_table",
                      help="APA xlsx table generated from s04/02_mostupdown_all_samples.r")
  
  parser$add_argument("-o", "--out_xlsx_file", type="character", required=TRUE,
                      dest="out_xlsx_file",
                      help="output excel file")
  
  args <- parser$parse_args()
} else {
  args <- data.table(ncpu=1,
                     pr_table='../../01_wkd/out/03_CallApa/HCT116_vs_DKO_APA.xlsx',
                     apa_table='../../01_wkd/out/04_AnnoApa/output/apa.ann_HCT116_vs_DKO.xlsx',
                     out_xlsx_file='../../01_wkd/out/04_AnnoApa/output/apa.ann_HCT116_vs_DKO_summary.xlsx')}

library(openxlsx)

mean_cvg <- 0.
pr_tb <- as.data.table(read.xlsx(args$pr_table))
pr_tb <- pr_tb[mean_HCT116_norm>mean_cvg | mean_DKO_norm>mean_cvg,]
apa_tb <- as.data.table(read.xlsx(args$apa_table))

ret <- apa_tb[,.N,by=c('call','tu_strand')]
ret$call2 <- ""
ret[call=='ShorterInB',call2:='distal_in_HCT116']
ret[call=='LongerInB',call2:='distal_in_DKO']

summary_tab <- data.table(pr=c(pr_tb[mean_HCT116_norm>mean_cvg,.N],pr_tb[mean_DKO_norm>mean_cvg,.N]),
                          genes=c(pr_tb[mean_HCT116_norm>mean_cvg,length(unique(gene_name))],pr_tb[mean_DKO_norm>mean_cvg,length(unique(gene_name))]),
                          pr_per_gene=c(mean(pr_tb[mean_HCT116_norm>0., list(pr_cnt = .N), by='gene_name']$pr_cnt),mean(pr_tb[mean_DKO_norm>0., list(pr_cnt = .N), by='gene_name']$pr_cnt)),
                          distal_watson=ret[tu_strand=='+',N],
                          distal_crick=ret[tu_strand=='-',N])

fields <- colnames(summary_tab)
samples <- c("HCT116","DKO")
summary_tab2 <- as.data.table(t(summary_tab))
colnames(summary_tab2) <- samples
summary_tab2$fields <- fields

message(sprintf("writing [%s] ...",args$out_xlsx_file))
write.xlsx(summary_tab2,file=args$out_xlsx_file)

sig_apa_dt <- apa_sig_res_annot[type!='1.NS',]
higher_pau_in_hct116 <- sig_apa_dt[,length(which(mean_a_frac>=mean_b_frac)),by='type']
higher_pau_in_hct116
higher_pau_in_dko <-sig_apa_dt[,length(which(mean_a_frac<mean_b_frac)),by='type']
higher_pau_in_dko
