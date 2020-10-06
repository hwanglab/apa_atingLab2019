# Annotate the APA results with the MostUp/MostDown information
library(argparse)
library(data.table)
library(ggplot2)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))

if (T) {
  parser <- ArgumentParser(description='mostupdown')
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
  										dest="ncpu", default = 1,
  										help="number of cpus to utilize [1]")
  
  parser$add_argument("-i", "--input_file", type="character", required=TRUE,
  										dest="input_file",
  										help="rscript workspace rd file")
   parser$add_argument("-c", "--control_tag", type="character", required=FALSE,
   										default="HCT116", dest="control_tag",
  										help="control group tag [HCT116]")
   
   parser$add_argument("-e", "--exp_tag", type="character", required=FALSE, 
   										default="DKO", dest="exp_tag",
  										help="experimental group tag [DKO]")

  parser$add_argument("-o", "--output_file", type="character", required=TRUE,
  										dest="output_file",
  										help="output file")
  args <- parser$parse_args()
} else {
  args <- data.table(ncpu=1,
                     input_file='../../01_wkd/out/03_CallApa/output/apa.sig.rd',
                     control_tag='HCT116',
                     exp_tag='DKO',
                     output_file='../../01_wkd/out/04_AnnoApa/output/apa.ann.rd')
}
#load("../03_CallApa/output/apa.sig.rd")
load(args$input_file)

# --------------------------------------------------------------------
# Function to take apa.sig object and call MostUp/Down (MUD)
# This function actually does both wtmean (wcall) and MUD (mcall)
callMUD <- function(apa)
{
	# wt mean method
	calls <- apa[,list(gene_name=gene_name[1],chr=chr[1],strand=strand[1],wa=weighted.mean(x=end,w=mean_a_frac),wb=weighted.mean(x=end,w=mean_b_frac)),by="tu"]
	calls$wb_minus_wa <- calls$wb-calls$wa

	# Call direction based on strand and value
	cp <- calls[strand=="+",]
	cp$wcall <- ""
	cp[wb>wa,]$wcall <- "LongerInB"
	cp[wb<wa,]$wcall <- "ShorterInB"
	stopifnot(cp$wcall!="")
	cm <- calls[strand=="-",]
	cm$wcall <- ""
	cm[wb>wa,]$wcall <- "ShorterInB"
	cm[wb<wa,]$wcall <- "LongerInB"
	stopifnot(cm$wcall!="")
	calls <- rbind(cp,cm)
	#calls[wcall=="",]$wcall <- "Error"

	# mud method
	ud <- apa[,list(MostUp=pr[b_minus_a_frac==max(b_minus_a_frac)],MostDown=pr[b_minus_a_frac==min(b_minus_a_frac)]),by="tu"]
	stopifnot(nrow(ud)==length(unique(apa$tu)))

	ud$chr <- apa[match(ud$MostUp,apa$pr),]$chr
	ud$up_start <- apa[match(ud$MostUp,apa$pr),]$start
	ud$up_end <- apa[match(ud$MostUp,apa$pr),]$end
	ud$down_start <- apa[match(ud$MostDown,apa$pr),]$start
	ud$down_end <- apa[match(ud$MostDown,apa$pr),]$end
	ud$strand <- apa[match(ud$MostDown,apa$pr),]$strand

	upp <- ud[strand=="+",]
	upp$mcall <- ""
	upp[up_start<down_start,]$mcall <- "ShorterInB"
	upp[up_start>down_start,]$mcall <- "LongerInB"
	stopifnot(upp$mcall!="")
	
	upm <- ud[strand=="-",]
	upm$mcall <- ""
	upm[up_start<down_start,]$mcall <- "LongerInB"
	upm[up_start>down_start,]$mcall <- "ShorterInB"
	stopifnot(upm$mcall!="")

	udmeth <- rbind(upp,upm)

	calls <- cbind(calls[match(udmeth$tu,calls$tu),],udmeth[,c("MostUp","MostDown","up_start","up_end","down_start","down_end","mcall"),with=F])
	calls$agree <- calls$wcall==calls$mcall
	calls
}

# Function to return the "res" table, keeping all PRs (rows) for TUs that have one or more int set sig PR
getSigSet <- function(myres)
{
	sigtu <- unique(myres[int_sig==TRUE,]$tu)
	myres$tu_sig <- myres$tu %in% sigtu
	myres <- myres[myres$tu_sig,]
	return(myres)
}
# --------------------------------------------------------------------

plot_pau_usage_comparison <- function(myres, apa.calls, control_tag,exp_tag,output_dir=NA) {
  # Save CSVs
  
  myres$type <- "NS"
  
  myres[int_sig==TRUE,]$type <- "significant"
  pdf_file <- file.path(output_dir,sprintf("pau_scatter_%s_vs_%s.pdf",control_tag,exp_tag))
  pdf(pdf_file)
  
  # browser()
  apa_calls <- apa.calls[,c("tu", "gene_name", "mcall")]
  myres_mcall <- merge(myres, apa_calls, by = c("tu","gene_name"), all.x=T, allow.cartesian=TRUE)
  
  myres_mcall$type2 <- "NS"
  myres_mcall[int_sig==TRUE & b_minus_a_frac<=0.,]$type2 <- control_tag
  myres_mcall[int_sig==TRUE & b_minus_a_frac>0.,]$type2 <- exp_tag
  
  # myres_mcall[int_sig==TRUE & mcall=="ShorterInB",]$type <- '2.Proximal'
  # myres_mcall[int_sig==TRUE & mcall=="LongerInB",]$type <- '3.Distal'

  ctrl_higher <- myres_mcall[type2==control_tag,gene_name]
  expr_higher <- myres_mcall[type2==exp_tag,gene_name]
  
  myres_mcall <- myres_mcall[order(type)]
  p <- ggplot(myres_mcall,aes(x=mean_a_frac,y=mean_b_frac)) + 
    geom_point(aes(color=type,alpha=factor(type)),size=1) + 
    scale_alpha_manual(values=c("NS"=0.01,"significant"=1.),guide="none") +
    scale_color_manual(values=c("#000000","#ff0000")) +
    labs(x=sprintf("%s Usage Fraction",control_tag),
         y=sprintf("%s Usage Fraction",exp_tag),
         title=sprintf('pA[%d]/genes[%d],%s\npA[%d]/genes[%d],%s',
                       length(ctrl_higher), #sig.proximal.in.DKO
                       length(unique(ctrl_higher)), #sig.proximal.uniq.genes.in.DKO
                       control_tag,
                       length(expr_higher), #sig.distal.in.DKO
                       length(unique(expr_higher)), #sig.distal.uniq.genes.in.DKO
                       exp_tag)
         ) +
    handy::ggnice()
  plot(p)
  dev.off()
  return(myres_mcall)
}

# --------------------------------------------------------------------
comp_tag <- sprintf('%s_vs_%s',args$control_tag,args$exp_tag)
message(sprintf('Start to analyze %s ...',comp_tag))

apa.int <- getSigSet(apa.sig[[comp_tag]]$res)
apa.calls <- callMUD(apa.int)
#save(apa.int,apa.calls,compress=T,file="output/apa.ann.rd")

save(apa.int,apa.calls,compress=T,file=args$output_file)

out_dir <- dirname(args$output_file)
apa_sig_res_annot <- plot_pau_usage_comparison(apa.sig[[comp_tag]]$res, 
                                               apa.calls, 
                                               args$control_tag, 
                                               args$exp_tag, 
                                               output_dir = out_dir)

pau_scatter_fpath <- file.path(out_dir,'plot_pau_usage_comp_scatter.tsv')

fwrite(apa_sig_res_annot,file=pau_scatter_fpath,sep='\t')
