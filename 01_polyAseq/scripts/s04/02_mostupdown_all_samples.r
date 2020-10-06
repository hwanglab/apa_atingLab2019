# Annotate the APA results with the MostUp/MostDown information
library(argparse)
library(data.table)

source(file.path(Sys.getenv('R_UTIL_APA'),'paddle.r'))
if (T) {

  parser <- ArgumentParser(description='mostupdown')

  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
  										dest="ncpu", default = 1,
  										help="number of cpus to utilize [1]")
  
  parser$add_argument("-i", "--input_file", type="character", required=TRUE,
  										dest="input_file",
  										help="rscript workspace rd file (apa.sig.rd) ")
  
  parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
  										dest="output_dir",
  										help="output dir to save apa.ann_*.rd files")
  args <- parser$parse_args()
} else {
  args <- data.table(ncpu=1,
                     input_file='../../01_wkd/out/03_CallApa/output/apa.sig.rd',
                     output_dir='../../01_wkd/out/04_AnnoApa/output')
}

library(openxlsx)

#load("../03_CallApa/output/apa.sig.rd")
load(args$input_file)

# --------------------------------------------------------------------
# Function to take apa.sig object and call MostUp/Down (MUD)
# This function actually does both wtmean (wcall) and MUD (mcall)
callMUD2 <- function(apa)
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
	
	ud$up_batchadj_pad <- apa[match(ud$MostUp,apa$pr),]$batchadj_padj
	ud$up_l2_b_over_a_frac <- apa[match(ud$MostUp,apa$pr),]$l2_b_over_a_frac
	
	ud$down_start <- apa[match(ud$MostDown,apa$pr),]$start
	ud$down_end <- apa[match(ud$MostDown,apa$pr),]$end
	ud$strand <- apa[match(ud$MostDown,apa$pr),]$strand
	ud$down_batchadj_pad <- apa[match(ud$MostDown,apa$pr),]$batchadj_padj
	ud$down_l2_b_over_a_frac <- apa[match(ud$MostDown,apa$pr),]$l2_b_over_a_frac
	
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

	calls <- cbind(calls[match(udmeth$tu,calls$tu),],udmeth[,c("MostUp","MostDown","up_start","up_end","up_batchadj_pad","up_l2_b_over_a_frac", "down_start","down_end","down_batchadj_pad","down_l2_b_over_a_frac","mcall"),with=F])
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

# --------------------------------------------------------------------
comp_names <- names(apa.sig)
for (i in 1:length(apa.sig)){
  apa.int <- getSigSet(apa.sig[[i]]$res)
  apa.calls <- callMUD2(apa.int)
  #save(apa.int,apa.calls,compress=T,file="output/apa.ann.rd")
  out_rd_file <- file.path(args$output_dir,sprintf('apa.ann_%s.rd',comp_names[[i]]))
  save(apa.int,apa.calls,compress=T,file=out_rd_file)
  out_excel_file <- file.path(args$output_dir,sprintf('apa.ann_%s.xlsx',comp_names[[i]]))
  #write.csv(apa.calls,row.names=F,file=csv)
  #reorder columns before printing
  
  setnames(apa.calls,"chr","tu_chr")
  setnames(apa.calls,"strand","tu_strand")
  setnames(apa.calls,"MostUp","pr_MostUp")
  setnames(apa.calls,"MostDown","pr_MostDown")
  setnames(apa.calls,"up_start","pr_up_start")
  setnames(apa.calls,"up_end","pr_up_end")
  setnames(apa.calls,"down_start","pr_down_start")
  setnames(apa.calls,"down_end","pr_down_end")
  setnames(apa.calls,"mcall","call")
  
  apa.calls[,c('wa','wb','wb_minus_wa','wcall','agree'):=NULL]
  setcolorder(apa.calls, c('gene_name','tu','tu_chr','tu_strand','pr_MostUp','pr_MostDown','pr_up_start','pr_up_end','up_batchadj_pad','up_l2_b_over_a_frac','pr_down_start','pr_down_end','down_batchadj_pad','down_l2_b_over_a_frac','call'))

  write.xlsx(apa.calls, file = out_excel_file, colNames = TRUE, quote=F, rowNames = F)
  message(sprintf('wrote [%s]',out_excel_file))
}
# --------------------------------------------------------------------
