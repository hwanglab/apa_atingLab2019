closeAllConnections()
rm(list=ls())
library(data.table)
library(openxlsx)
library(reshape2)
library(ggplot2)
library(Rtsne)
library(rgl)
library(scatterplot3d)
library(car)
library(randomcoloR)
library("extrafont")

source(file.path(Sys.getenv('R_UTIL'),'lib_clustering.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))

loadfonts()

gen_heatmap2_gene <- function(dt2,pdf_fn=NA) {
  
  # browser()
  
  mycol <- colorpanel(20, "red", "white", "blue")
  if (!is.na(pdf_fn)) {pdf(pdf_fn,width=14,height = 3)}
  
  ret2 <- t(as.data.table(strsplit(dt2$gene_test_id,"\\[")))
  dt2$gene <- unlist(ret2[,1])
  dt_by_gene <- split(dt2,by='gene')
  ret <- lapply(dt_by_gene,function(dbg) {
    return(colMeans(dbg[,2:13]))
  })
  
  gxs_mat <- as.matrix(as.data.table(ret))
  colnames(gxs_mat) <- names(ret)
  rownames(gxs_mat) <- colnames(dt_by_gene[[1]])[2:13]
  
  x <- 1
  if (T) {
    # heatmap.2(gxs_mat,
    #           Rowv=as.dendrogram(h_fe),
    #           Colv=as.dendrogram(h_sm),
    #           col=mycol,
    #           lmat=rbind(c(2),c(3),c(1),c(4)), 
    #           lhei=c(1,1,9,0), 
    #           lwid=c(1),
    #           scale="none",
    #           #density.info="none",
    #           trace="none",
    #           dendrogram = "none")
    
    ret_hm <- heatmap.2(gxs_mat,
                        Rowv=F,
                        Colv=T,
                        col=mycol,
                        lhei=c(3,9),
                        lwid=c(2,10),
                        # lmat=rbind(c(2),c(3),c(1),c(4)), 
                        # lhei=c(1,1,9,0), 
                        # lwid=c(1),
                        scale="none",
                        #density.info="none",
                        trace="none",
                        srtRow = 45,
                        srtCol = 45,
                        dendrogram = "none")
    
  }
  x <- 1
  # browser()
  if (!is.na(pdf_fn)) {dev.off()}
  
  return(list(M=gxs_mat,ret_hm=ret_hm))
}

summary_rd_fn <- '../546goi_all/top/corr_gtop8000_summary_filtered.rd'

load(summary_rd_fn) #sel_genes_w,sel_genes_l

set.seed(42)

# browser()

message('preparing matrix ...')
#sel_genes_w[,merged:=NULL]

M2 <- as.matrix(sel_genes_w[,2:13])
cohorts <- colnames(sel_genes_w[,2:13])
prob_pos <- sel_genes_w$probe_pos
M2t <- t(M2)

j_indicator <- apply(M2t, 2, var) != 0
M <- M2t[ , j_indicator]
rownames(M) <- cohorts
colnames(M) <- prob_pos[j_indicator]

sel_genes_w <- sel_genes_w[j_indicator,]


# -------------------------
jmessage('heatmap2 ....')

# pos_by_gene <- sel_genes_w[,list(probe_pos,gene_test_id)]

message("instruction run clusterings in the gene-level first, then collect cluster results from an x-axis and run gen_heatmap2_transfer_genecluster_to_bp2()")

ghm <- gen_heatmap2_gene(sel_genes_w,
                  pdf_fn=NA)

genes_order_to_print <- rownames(ghm$ret_hm$carpet)
fxs_col <- gen_heatmap2_transfer_genecluster_to_bp2(M,
                                         genes_order_to_print,
             # pdf_fn='./output/heatmaps/heatmap_cluster_gene_to_bp2.pdf',
             # zi_pdf_fn='./output/heatmaps/heatmap_cluster_gene_to_bp2_zi.pdf',
             dt=sel_genes_w)

cluster_sorted_w <- sel_genes_w[match(fxs_col$pos,sel_genes_w$probe_pos),]
cluster_sorted_pval_w <- sel_genes_pval_w[match(fxs_col$pos,sel_genes_w$probe_pos),]

tag <- 'filtered'
wb <- createWorkbook()
sheet_name <- sprintf('pearson_r2_%s',tag)
addWorksheet(wb,sheetName = sheet_name)
writeData(wb,sheet=sheet_name,x=cluster_sorted_w)

sheet_name <- sprintf('pearson_r2_pval_%s',tag)
addWorksheet(wb,sheetName = sheet_name)
writeData(wb,sheet=sheet_name,x=cluster_sorted_pval_w)

xlsx_file <- sprintf('./output/corr_gtop8000_summary_max_per_cohort_%s_rep2_cls.xlsx',tag)
saveWorkbook(wb,xlsx_file,overwrite=TRUE)

