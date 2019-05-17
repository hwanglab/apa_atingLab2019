setwd("~/projects/apa/src/s12_tcga_summary")
closeAllConnections()
rm(list=ls())

library(argparse)
library(data.table)
library(reshape2)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))

rd_fn <- '../546goi_all/corr_gtop8000.rd'
load(rd_fn)

if (!dir.exists('./output')) {
  dir.create('./output')
}


cohort_names <- names(corr_mpts)

pccs <- list()
for (i in 1:length(cohort_names)){
  cohort <- cohort_names[[i]]
  pcc_i <- rbindlist(corr_mpts[[i]])
  
  pcc_i$cohort <- cohort
  pccs[[cohort]] <- pcc_i[,list(probe_pos,pearson.r2,pearson.pval,data_points,test_id,gene_name,cohort)]
}

pcc <- rbindlist(pccs)
rm(pccs)
rm(pcc_i)

filter_opts <- list(no_filter=c(1.1, 0,'original'), #p.val, abs(corr)
                    filter=c(0.05,0.1,'filtered'))

for (i in 1:length(filter_opts)) {
  
  pcc <- pcc[data_points>=100,]
  
  pcc_pval_co <- as.numeric(filter_opts[[i]][1])
  pcc_co <- as.numeric(filter_opts[[i]][2])
  tag <- filter_opts[[i]][3]
  
  pcc <- pcc[(pearson.pval<pcc_pval_co & abs(pearson.r2)>pcc_co),]
  
  max_cc_per_cohgen <- as.data.table(pcc[, .SD[which.max(abs(pearson.r2))], by=c('cohort','gene_name')])
  test_id_map <- unique(max_cc_per_cohgen[,list(test_id,gene_name)])
  
  sel_prob_pos <- unique(max_cc_per_cohgen[,list(test_id,probe_pos)])
  
  max_cc_per_cohgen <- pcc[(test_id %in% sel_prob_pos$test_id & probe_pos %in% sel_prob_pos$probe_pos),]
  
  sel_genes_l <- max_cc_per_cohgen
  
  max_cc_per_cohgen[cohort=='merged',cohort:='Vmerged']
  
  sel_genes_w <- as.data.table(dcast(max_cc_per_cohgen, test_id + probe_pos ~ cohort, value.var="pearson.r2"))
  sel_genes_w$gene_name <- test_id_map[match(sel_genes_w$test_id,test_id_map$test_id),gene_name]
  sel_genes_w$gene_test_id <- paste0(sel_genes_w$gene_name,'[',sel_genes_w$test_id,']')
  
  sel_genes_w$chr <- sapply(strsplit(sel_genes_w$probe_pos,':'),function(chr_pos){chr_pos[[1]]})
  sel_genes_w$start <- sapply(strsplit(sel_genes_w$probe_pos,':'),function(chr_pos){as.integer(chr_pos[[2]])})
  sel_genes_w <- sel_genes_w[order(chr,start,gene_test_id)]
  setnames(sel_genes_w,'Vmerged','merged')
  sel_genes_w$chr <- NULL 
  sel_genes_w$start <- NULL
  sel_genes_w$test_id <- NULL
  sel_genes_w$gene_name <- NULL
  
  sel_genes_w$num_valid_cohorts <- rowSums(!is.na(sel_genes_w[,2:13]))
  sel_genes_w[is.na(sel_genes_w)] <- 0.
  
  
  sel_genes_pval_w <- as.data.table(dcast(max_cc_per_cohgen, test_id + probe_pos ~ cohort, value.var="pearson.pval"))
  sel_genes_pval_w$gene_name <- test_id_map[match(sel_genes_pval_w$test_id,test_id_map$test_id),gene_name]
  sel_genes_pval_w$gene_test_id <- paste0(sel_genes_pval_w$gene_name,'[',sel_genes_pval_w$test_id,']')
  
  sel_genes_pval_w$chr <- sapply(strsplit(sel_genes_pval_w$probe_pos,':'),function(chr_pos){chr_pos[[1]]})
  sel_genes_pval_w$start <- sapply(strsplit(sel_genes_pval_w$probe_pos,':'),function(chr_pos){as.integer(chr_pos[[2]])})
  sel_genes_pval_w <- sel_genes_pval_w[order(chr,start,gene_test_id)]
  setnames(sel_genes_pval_w,'Vmerged','merged')
  sel_genes_pval_w$chr <- NULL 
  sel_genes_pval_w$start <- NULL
  sel_genes_pval_w$test_id <- NULL
  sel_genes_pval_w$gene_name <- NULL
  
  sel_genes_pval_w$num_valid_cohorts <- rowSums(!is.na(sel_genes_pval_w[,2:13]))
  sel_genes_pval_w[is.na(sel_genes_pval_w)] <- 0.
  
  summary_rd_fn <- sprintf('../546goi_all/top/corr_gtop8000_summary_%s.rd',tag)
  save(sel_genes_w,sel_genes_pval_w,sel_genes_l,file=summary_rd_fn,compress=T)
  
  wb <- createWorkbook()
  sheet_name <- sprintf('pearson_r2_%s',tag)
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet=sheet_name,x=sel_genes_w)
  
  sheet_name <- sprintf('pearson_r2_pval_%s',tag)
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet=sheet_name,x=sel_genes_pval_w)
  
  xlsx_file <- sprintf('output/corr_gtop8000_summary_max_per_cohort_%s_rep2.xlsx',tag)
  saveWorkbook(wb,xlsx_file,overwrite=TRUE)
}