library(argparse)
library(data.table)

debug <- T
if (!debug) {
  parser <- ArgumentParser(description='locate 546 genes of interest in MAnorm M density plot [hongc2@ccf.org]')

  parser$add_argument("-e", "--extbp", type="integer", required=FALSE,
                      dest="extbp", default = 5e3,
                      help="extension basepair from TU region [5000]")
  
  parser$add_argument("-w", "--prepropa", type="character", required=FALSE,
                      dest="prepropa", default = "../01_polyAseq/01_wkd/out/03_CallApa/output/prepropa.rd",
                      help="prepropa.rd [../01_polyAseq/01_wkd/out/03_CallApa/output/prepropa.rd]")
  
  parser$add_argument("-a", "--apa_ann", type="character", required=FALSE,
                      dest="manorm_xlsx", default = "../01_polyAseq/01_wkd/out/04_AnnoApa/output/apa.ann.rd",help="apa_ann.rd [../01_polyAseq/01_wkd/out/04_AnnoApa/output/apa.ann.rd]")
  
  parser$add_argument("-m", "--manorm_xlsx", type="character", required=FALSE,
                      dest="apa_ann", default = "fastq/kundaje_encode/Supplementary_Table_S5.xlsx",help="M_Val xlsx file path [fastq/kundaje_encode/manorm/Supplementary_Table_S5.xlsx]")
  
  
  parser$add_argument("-M", "--mode", type="character", required=FALSE,
                      dest="mode", default = "m_value",help="[m_value]")

  parser$add_argument("-g", "--gc2print", type="integer", required=FALSE,
                      dest="gc2print", default = 30,help="[30]")
  
  parser$add_argument("-o", "--outprefix", type="character", required=FALSE,
                      dest="outprefix", default="fastq/kundaje_encode/goi_apa_manorm", help="output prefix")

} else {
  args = data.table(extbp=5e3,
                    prepropa="../01_polyAseq/01_wkd/out/03_CallApa/output/prepropa.rd",
                    apa_ann="../01_polyAseq/01_wkd/out/04_AnnoApa/output/apa.ann.rd",
                    manorm_xlsx="fastq/kundaje_encode/manorm/Supplementary_Table_S5.xlsx",
                    outprefix = 'fastq/kundaje_encode/goi_apa_manorm',
                    mode = 'm_value',
                    gc2print=30)
}

# ----------------------
message('loading libraries ...')
library(openxlsx)
library(goldmine)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(gridExtra)

# ----------------------
message('load polyA-seq rd files ...')
load(args$prepropa)
load(args$apa_ann)

# ----------------------
message('defining extended TU regions of interest ...')

goi_ext_tu_gr <- red_ens$tu[red_ens$tu$tu %in% apa.calls$tu,]
start(goi_ext_tu_gr) <- start(goi_ext_tu_gr) - args$extbp
end(goi_ext_tu_gr) <- end(goi_ext_tu_gr) + args$extbp

goi_ext_tu_gr_rd <- paste0(args$outprefix,'.goi.rd')
save(goi_ext_tu_gr,file=goi_ext_tu_gr_rd,compress = T)
message(sprintf('[%s] saved',goi_ext_tu_gr_rd))

# ----------------------
message('loading MAnorm tables ...')
manorm <- list()

out_pdf <- paste0(args$outprefix,'.pdf')
out_xlsx <- paste0(args$outprefix,'.xlsx')
out_rd <- paste0(args$outprefix,'.rd')

if (!debug) {pdf(out_pdf)}

S <- 6 #total 6 ChIP-seq tags; the last will be column keys
manorms <- list()

wb <- createWorkbook("MAnorm_546goi_diff_binding_analysis")

for (s in 1:S) {
  
  message(sprintf('reading sample sheet [%d/%d] from [%s]...',s,S,args$manorm_xlsx))
  
  manorm$dt <- as.data.table(read.xlsx(args$manorm_xlsx,sheet=s))
  exp_group <- tstrsplit(names(manorm$dt)[9],'_')[[5]]
  protein <- tstrsplit(names(manorm$dt)[9],'_')[[6]]
  ctrl_group <- tstrsplit(names(manorm$dt)[10],'_')[[5]]
  comp_tag <- sprintf('%s_%s',ctrl_group,exp_group)
  
  sheet_name <- sprintf("%s_%s",comp_tag,protein)
  addWorksheet(wb,sheetName = sheet_name)
  message(sprintf('working on a sheet [%s]',sheet_name))
  message('converting manorm dt to gr format ...')
  manorm$gr <- makeGRanges(manorm$dt,strand=F)
  
  message('annotating gene name if the diff binding sites are overlapped with our 546 goi ...')
  qi_sj <- findOverlaps(manorm$gr,goi_ext_tu_gr,ignore.strand=T)
  manorm$dt$goi546 <- as.character()
  manorm$dt$goi546[from(qi_sj)] <- goi_ext_tu_gr$name[to(qi_sj)]

  #message('counting overlapped region between manorm peaks and genes of interest (TU regions) ...')
  #manorm$dt$overlapped_by_goi <- countOverlaps(manorm$gr,goi_ext_tu_gr,ignore.strand=T)
  
  manorm$dt$group <- 'not_affected'
  manorm$dt[!is.na(goi546),group:='APA_perturbed_by_CH3']
  
  message('preparing X and Y axis value ...')
  manorm$dt$abs_M <- abs(manorm$dt$M_value)
  manorm$dt$nlog10pv <- -log10(manorm$dt$P_value)

  manorm$dt <- manorm$dt[group=='APA_perturbed_by_CH3',]
  
  message('finding cutoffs for top 25 pct ...')
  manorm$dt$m_x_pv <- manorm$dt$abs_M * manorm$dt$nlog10pv
  
  manorms[[sheet_name]] <- manorm
  
  manorm$dt <- manorm$dt[order(-m_x_pv)]
  manorm$dt$ranking <- 1:nrow(manorm$dt)

  M <- dim(manorm$dt)[1]
  top25r <- as.integer(round(0.25 * M))

  message("applying logFold2Change cutoff and -log10(P-value) cutoff for plotting...")
  manorm$dt$label2 <- as.character()
  manorm$dt[1:top25r,label2:='top25_pct']

  message('preparing which dot is labeled by gene name ...')
  if (args$gc2print > top25r) {
    args$gc2print <- top25r
  }
  
  both_cond <- manorm$dt[label2=='top25_pct',]
  ranking_best_per_gene <- both_cond[,min(ranking),by=goi546][1:args$gc2print,]
  
  manorm$dt$print_flag <- FALSE
  j <- which(manorm$dt$ranking %in% ranking_best_per_gene$V1)
  manorm$dt[j,print_flag:=TRUE]
  
  bottom_printed_gene <- manorm$dt[j[length(j)],]
  
  message('plotting histograms ...')
  
  if (T) {
    if (args$mode == 'm_value') {
      
      bottom_printed_gene$abs_M

      p <- ggplot(manorm$dt,aes(x=M_value,y=nlog10pv,group=label2))
      p <- p + 
        geom_point(aes(shape = ".",color=label2)) +
        geom_rug(sides="b",alpha=0.25) +
        ggtitle(sprintf("[%s:%s]\nmin|M| x min|nlog10(pv)|=%g x %g = %g",comp_tag,protein,bottom_printed_gene$abs_M,bottom_printed_gene$nlog10pv,bottom_printed_gene$m_x_pv)) +
        geom_text_repel(aes(x=M_value,
                            y=nlog10pv,
                            label=ifelse(print_flag == T, goi546,""))) +
        labs(x='APA M (log2 Fold Change)',y='-log10(p-value)')

    } else {
      p <- ggplot(data=manorm$dt, aes(x=value,
                                      fill=group,
                                      color=group)) + 
        geom_density(alpha=0.2,lwd=0.8,adjust=0.5) +
        ggtitle(sprintf("-log10(P-value) histogram [%s:%s]\npct lines(%3.2f@25%%,%3.2f@50%%,%3.2f@75%%,%3.2f@90%%)",comp_tag,protein,qtile[1],qtile[2],qtile[3],qtile[4])) +
        geom_vline(xintercept =qtile, linetype="dotted", 
                   color = "black", size=0.5) +
        theme(legend.position="bottom") +
        #xlim(min(manorm$dt$value),max(manorm$dt$value)*0.15)
        xlim(min(manorm$dt$value),10)
    }
    print(p)
  }
  
  writeData(wb,sheet_name,manorm$dt,rowNames=TRUE)

}
saveWorkbook(wb,out_xlsx,overwrite = TRUE)
save(manorms,file=out_rd,compress=T)
if (!debug) {dev.off()}
