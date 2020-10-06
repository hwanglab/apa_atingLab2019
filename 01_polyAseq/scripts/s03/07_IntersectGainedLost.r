# Save all the parts we'd like to have for the purposes of screening the APA calls
library(reshape)
library(matrixStats)

library(handy)
library(argparse)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

parser <- ArgumentParser(description='geneplots')

parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
										dest="ncpu", default = 1,
										help="number of cpus to utilize [1]")

parser$add_argument("-i", "--input_file", type="character", required=TRUE,
										dest="input_file",
										help="rscript workspace rd file")

parser$add_argument("-c", "--control_tag", type="character", required=TRUE,
										dest="control_tag",
										help="control sample tag")

parser$add_argument("-e", "--exp_tag", type="character", required=TRUE,
										dest="exp_tag",
										help="experiment sample tag")

parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
										dest="output_dir",
										help="output dir")
args <- parser$parse_args()

#load("output/apa.sig.rd")
load(args$input_file)

# Save CSVs
for(i in 1:length(apa.sig))
{
	message(names(apa.sig)[i])
	myres <- apa.sig[[i]]$res
	myres$type <- "NS"
	myres[(batchadj_perfc_sig==TRUE)&(batchadj_delta_sig==FALSE)]$type <- "pr_perfc_uniq"
	myres[(batchadj_perfc_sig==FALSE)&(batchadj_delta_sig==TRUE)]$type <- "pr_delta_uniq"
	myres[int_sig==TRUE,]$type <- "int_sig"
	myres <- myres[type!="NS",]

	#pdf(file="output/IntersectXY.pdf")
	pdf(file=file.path(args$output_dir,sprintf("IntersectXY_%02d.pdf",i)))
	#print(ggplot(myres,aes(x=mean_a_frac,y=mean_b_frac,color=type)) + handy::ggnice() + geom_point() + labs(x="HCT Usage Fraction",y="DKO Usage Fraction"))
	print(ggplot(myres,aes(x=mean_a_frac,y=mean_b_frac,color=type)) + handy::ggnice() + geom_point() + labs(x=sprintf("%s Usage Fraction",args$control_tag),y=sprintf("%s Usage Fraction",args$exp_tag)))
	#print(ggplot(myres,aes(x=b_minus_a_frac,y=l2_b_over_a_frac,color=type)) + handy::ggnice() + geom_point() + labs(x="DKO Fraction-HCT Fraction",y="log2 DKO Frac/HCT Frac"))
	print(ggplot(myres,aes(x=b_minus_a_frac,y=l2_b_over_a_frac,color=type)) + handy::ggnice() + geom_point() + labs(x=sprintf("%s Fraction-%s Fraction",args$exp_tag,args$control_tag),y=sprintf("log2 %s Frac/%s Frac",args$exp_tag,args$control_tag)))
	dev.off()
}
