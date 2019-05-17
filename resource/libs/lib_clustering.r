library(ggplot2)
library(cluster)    # clustering algorithms
library(dendextend) # for comparing two dendrograms
library(factoextra)
library(gridGraphics)
library(gplots)
library(grid)
library(gridExtra)

grab_grob <- function(){
  grid.echo()
  grid.grab()
}

gen_heatmap2_transfer_genecluster_to_bp2 <- function(fxs_mat,gene_order0,pdf_fn=NA,zi_pdf_fn=NA,dt2=NA) {
  
  # browser()
  message('read order of genes from gene clustering result ...')
  #gene_order <- fread('./output/heatmaps/gene_cluster_by_ColvT_order.tsv')
  
  gene_order <- data.table(gene=gene_order0,
                           order2=1:length(gene_order0))
  
  message('extracting genes from the corr table ...')
  dt2$gene <- unlist(str_extract_all(dt2$gene_test_id,"[:alnum:]+(?=\\[)"))
  dt2$test_id <- unlist(str_extract_all(dt2$gene_test_id,"(?<=\\[)[:digit:]+(?=\\])"))

  message('count the number of probes in each gene ...')
  probe_cnt <- dt2[,list(cnt=.N),by=gene]
  
  probe_cnt$order2 <- gene_order[match(probe_cnt$gene,gene_order$gene),order2]
  
  probe_cnt <- probe_cnt[order(order2)]
  
  message('reformatting gene to probe template ...')
  idx_order <- lapply(1:dim(probe_cnt)[1], function(i) {
    return(data.table(gene=rep(probe_cnt[i,gene],probe_cnt[i,cnt]),
                      order2=rep(probe_cnt[i,order2],probe_cnt[i,cnt])))
  })
  idx_order <- rbindlist(idx_order)
  rm(probe_cnt)
  
  message('preparing gene label to print to reduce the complexity in the axis')
  idx_order_by_gene <- split(idx_order,by='gene')
  
  genes <- names(idx_order_by_gene)
  idx_order_by_gene2 <- list()
  for (i in 1:length(genes)) {
    rj <- idx_order_by_gene[[i]]
    rj$gene2 <- '.'
    M <- as.integer(ceiling(dim(rj)[1]/2.))
    rj[M,gene2:=genes[[i]]]
    idx_order_by_gene2[[genes[[i]]]] <- rj
  }
  idx_order <- rbindlist(idx_order_by_gene2)
  
  
  message('processing fxs_mat columns ...')
  fxs_col <- data.table(pos=colnames(fxs_mat),
                        idx=1:dim(fxs_mat)[2],
                        test_id=0)
  
  #get pos - gene mapper
  probe_gene <- dt2[,list(probe_pos,gene,test_id)]
  
  #append gene name to fxs_col
  fxs_col$test_id <- probe_gene[match(fxs_col$pos, probe_gene$probe_pos),test_id]
  fxs_col$gene <- probe_gene[match(fxs_col$pos, probe_gene$probe_pos),gene]

  #append cluster order to fxs_col
  fxs_col$order2 <- gene_order[match(fxs_col$gene,gene_order$gene),order2]
  fxs_col <- fxs_col[order(order2,test_id)]

  idx_order <- idx_order[order(order2)]
  # browser()
  fxs_mat2 <- fxs_mat[,fxs_col$idx]
  colnames(fxs_mat2) <- idx_order$gene2

  fxs_zis <- list()
  zi_genes <- c('HEATR2', 'MBD1', 'FANCL', 'IDUA', 'FBLN1', 'AFAP1', 'AP2A2')
  zi_genes <- sort(zi_genes)
  # browser()
  
  col.pal <- colorpanel(20, "red", "white", "blue")
  col.breaks <- seq(-1,1,length=21)
  
  if (!is.na(zi_pdf_fn)) {pdf(zi_pdf_fn,width=14,height=14)}
  par(cex.main=0.5)
  # par(mar=c(1,1,1,1))
  colw <- 9/11
  gl <- lapply(zi_genes, function(zi_gene) {
    pidx <- which(fxs_col$gene==zi_gene)
    
    fxs_zi <- fxs_mat2[,pidx]
    colnames(fxs_zi) <- fxs_col[pidx,pos]
    N <- dim(fxs_zi)[2]
    # browser()
    heatmap.2(fxs_zi,
              Rowv=F,
              Colv=F,
              key = F,
              col = col.pal,
              breaks = col.breaks,
              scale="none",
              density.info="none",
              trace="none",
              dendrogram = "none",
              lhei=c(3, 9),
              lwid=c(12-N*colw, N*colw),
              #margins=c(12,12),
              cexRow=0.7,
              cexCol=0.7,
              srtRow = 45,
              srtCol = 45,
              main=zi_gene)
    grab_grob()
  })
  
  # grid.newpage()
  browser()
  grid.arrange(grobs=gl, ncol=7, clip=TRUE)
  if (!is.na(zi_pdf_fn)) {dev.off()}
  
  # browser()
  if (!is.na(pdf_fn)) {pdf(pdf_fn,width=14,height = 3)}
  
  if (F) {
    heatmap.2(fxs_mat2,
              Rowv=T,
              Colv=F,
              key = F,
              col=col.pal,
              scale="none",
              lmat=rbind(c(2),c(3),c(1),c(4)),
              lhei=c(1,1,9,0),
              lwid=c(1),
              #density.info="none",
              trace="none",
              dendrogram = "none")
  }
  
  #w/o a legend key
  heatmap.2(fxs_mat2,
            Rowv=F,
            Colv=F,
            key = F,
            col=col.pal,
            scale="none",
            density.info="none",
            trace="none",
            dendrogram = "none",
            lhei=c(1, 11),
            lwid=c(1, 11),
            srtRow = 45,
            srtCol = 45)
  
  x <- 1
  if (F) {
  #w/o a legend key
  heatmap.2(fxs_mat2,
            Rowv=F,
            Colv=F,
            key = F,
            lmat=rbind(c(2),c(3),c(1),c(4)), 
            lhei=c(1,1,9,0), 
            lwid=c(1),
            col=mycol,
            scale="none",
            density.info="none",
            trace="none",
            dendrogram = "none")
  }
  if (!is.na(pdf_fn)) {dev.off()}
  
  x <-1
  return(fxs_col)
}


gen_heatmap2 <- function(fxs_mat,sm_measurement="pearson",fe_measurment="pearson",pdf_fn=NA) {
  
  browser()
  
  sm_dist <- get_dist(t(fxs_mat), method = sm_measurement) #column-wise (between basepairs)
  fe_dist <- get_dist(fxs_mat, method = fe_measurment) #row-wise (between feature)
  
  h_sm <- hclust(sm_dist, method = "complete", members=NULL)
  h_fe <- hclust(fe_dist, method = "complete", members=NULL)
  
  mycol <- colorpanel(40, "darkblue", "yellow", "white")
  
  mycls <- cutree(h_fe, h=max(h_fe$height)/1.5)
  mycolh_sm <- rainbow(length(unique(mycls)), start=0.1, end=0.9)
  mycolh_sm <- mycolh_sm[as.vector(mycls)] 
  
  if (!is.na(pdf_fn)) {pdf(pdf_fn)}
  
  heatmap.2(fxs_mat, 
            Rowv=as.dendrogram(h_fe),
            Colv=as.dendrogram(h_sm),
            col=mycol,
            scale="row",
            density.info="none",
            trace="none", 
            RowSideColors=mycolh_sm)
  
  
  
}