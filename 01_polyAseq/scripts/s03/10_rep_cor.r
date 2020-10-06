# Annotate the APA results with the MostUp/MostDown information
library(argparse)
library(data.table)

args_tmp <- commandArgs(trailingOnly = F)
scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
source_root <- dirname(dirname(scriptPath))

if (T) {
  parser <- ArgumentParser(description='10_rep_cor')
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                      dest="ncpu", default = 1,
                      help="number of cpus to utilize [1]")
  
  parser$add_argument("-i1", "--pre_apa_rd", type="character", required=TRUE,
                      dest="pre_apa_rd",
                      help="rscript workspace rd file (prepropa.rd)")
  
  parser$add_argument("-i2", "--apa_sig_rd", type="character", required=TRUE,
                      dest="apa_sig_rd",
                      help="rscript workspace rd file (apa_sig.rd)")
  
  parser$add_argument("-c", "--control_tag", type="character", required=FALSE, default="HCT116",
                      dest="control_tag",
                      help="control group tag [HCT116]")
  
  parser$add_argument("-e", "--exp_tag", type="character", required=FALSE, default="DKO",
                      dest="exp_tag",
                      help="experimental group tag [DKO]")
  
  parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
                      dest="output_dir",
                      help="output directory")
  args <- parser$parse_args()
} else {
  args <- data.table(ncpu=1,
                     pre_apa_rd='../../01_wkd/out/03_CallApa/output/prepropa.rd',
                     apa_sig_rd='../../01_wkd/out/03_CallApa/output/apa.sig.rd',
                     control_tag='HCT116',
                     exp_tag='DKO',
                     output_dir='../../01_wkd/out/03_CallApa/output')
}


# Save all the parts we'd like to have for the purposes of screening the APA calls

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))
load(args$pre_apa_rd)
load(args$apa_sig_rd)

library(reshape)
library(matrixStats)
library(methylaction)
#library(ggbio)

comp_tag <- sprintf('%s_vs_%s',args$control_tag,args$exp_tag)
message(sprintf('Start to analyze %s ...',comp_tag))

# Want XY correlation plots among replicates for HCT/DKO
mat <- apa.sig[[comp_tag]]$counts.norm
#mel <- melt(mat)
#colnames(mel) <- c("pr","sample","count")

mat <- data.table(mat)

library(GGally)

M <- cor(mat)
png_file <- file.path(args$output_dir,"RepCor.png")
message(sprintf('printing a corr plot among polyA tech replicates [%s] ...',png_file))

png(filename=png_file,width=8,height=8,units="in",res=600)
#pdf(file=pdf_file,width=5,height=5)
#ggplot(mat,aes(x=HCT116_1,y=DKO_2)) + geom_point() + handy::ggnice()
pm <- ggpairs(mat)
print(pm)
dev.off()

library(corrplot)
pdf_file <- file.path(args$output_dir,"RepCorHeat.pdf")
message(sprintf('printing a corr plot among polyA tech replicates [%s]...',pdf_file))
pdf(file=pdf_file,width=5,height=5)
corrplot(M, method="color")
corrplot(M, method="number")
corrplot(M, type="lower", method="number")
dev.off()

