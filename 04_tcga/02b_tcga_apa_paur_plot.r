# Rscript /home/hongc2/projects/apa/src/s11_tcga/03_tcga_analy_apa_param_opt3_combined2b_s10_4corrs_module.r -n 1 -t /home/hongc2/projects/apa/src/s11_tcga/debug -c BLCA,BRCA,COAD,ESCA,GBM,KIRC,LUAD,LUSC,PRAD,SKCM,STAD,UCEC -T 8000 -b /home/hongc2/projects/apa/src/s11_tcga/input/apa546_2018_11_16.xlsx -o /media/sammy/s11_tcga_debug/546goi/topk -s gtop8000 -i 2

closeAllConnections()
rm(list=ls())

library(argparse)
library(data.table)

debug <- F
if (!debug) {
  parser <- ArgumentParser(description='TCGA analysis [hongc2@ccf.org]')
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                      dest="ncpu", default = 1,
                      help="number of cpus to utilize [1]")
  
  parser$add_argument("-f", "--min_fpkm_uq", type="integer", required=FALSE,
                      dest="min_fpkm_uq", default = 0,
                      help="minimum fpkm-uq value to test PAU [0]")
  
  parser$add_argument("-t", "--tcga_rd_prefix", type="character", required=TRUE,
                      default='output/03_tcga_analy_apa/546goi_all', dest="tcga_rd_prefix",
                      help="TCGA rd prefix")
  
  parser$add_argument("-c", "--cohorts", type="character", required=FALSE,
                      default='BLCA,BRCA,COAD,ESCA,KIRC,LUAD,LUSC,PRAD,SKCM,STAD,UCEC', dest="cohorts",
                      help="TCGA cohorts in abbreviation, e.g., DLBC or DLBC,ESCA")
  
  parser$add_argument("-b", "--goi_file", type="character", required=TRUE,
                      dest="goi_file",
                      help="xlsx file of genomic region of interest in the format following, chr shared_start shared_end uniq_start uniq_end goi_name")
  
  parser$add_argument("-s", "--method", type="character", required=FALSE,
                      dest="method", default="gtop8000",
                      help="pau vs. beta combine method")
  
  parser$add_argument("-T", "--expr_topk", type="integer", required=FALSE,
                      default=8000, dest="expr_topk",
                      help="expr_topk [10000]: top [10000] expressed samples for each gene [10000], Set to 0 to use min_mrn_md instead.")
  
  parser$add_argument("-rd1", "--rd1", type="character", required=TRUE,
                      dest="prepropa_rd",
                      help="prepropa.rd from polyA pipeline output file")
  
  parser$add_argument("-rd2", "--rd2", type="character", required=TRUE,
                      dest="apa_ann_rd",
                      help="apa.ann.rd from polyA pipeline output file")
  
  parser$add_argument("-o", "--out_dir", type="character", required=FALSE,
                      default='./output', dest="out_dir",
                      help="output directory for pdf files")
  
  parser$add_argument("-r", "--step_from", type="integer", required=FALSE,
                      default=1, dest="step_from",
                      help="step_from [1]: from 1st step")
  
  
  parser$add_argument("-u", "--min_obs", type="integer", required=FALSE,
                      default=100, dest="min_obs",
                      help="min_obs")
  
  parser$add_argument("-d", "--debug", type="integer", required=FALSE,
                      default=0, dest="debug",
                      help="reuse [0]:no debug, 1:debug")
  
  args <- parser$parse_args()
} else {
  
  # args <- data.table(ncpu=8,
  #                    tcga_rd_prefix='/home/hongc2/projects/apa/src/s11_tcga/debug',
  #                    cohorts='BLCA,BRCA,COAD,DLBC,ESCA,GBM,KIRC,LUAD,LUSC,PRAD,SKCM,STAD,UCEC',
  #                    out_dir='./output/03_tcga_analy_apa/3goi_all/topk',
  #                    method='gtop5000',
  #                    goi_file='/home/hongc2/projects/apa/src/s11_tcga/input/pau_testing_nfya_heatr2_dnajb5_input_hct116_vs_dko.xlsx',
  #                    expr_topk=5000,
  #                    min_mrn_md=0,
  #                    fit_model='linear',
  #                    min_fpkm_uq=0,
  #                    step_from=10,
  #                    debug=1)
  # 
  args <- data.table(ncpu=8,
                     tcga_rd_prefix='./output/03_tcga_analy_apa/546goi',
                     cohorts='BLCA,BRCA,COAD,ESCA,KIRC,LUAD,LUSC,PRAD,SKCM,STAD,UCEC',
                     out_dir='./output/03_tcga_analy_apa/546goi_all',
                     method='gtop8000',
                     goi_file='/home/hongc2/projects/apa/src/s11_tcga/input/apa546_2018_11_16.xlsx',
                     expr_topk=8000,
                     min_mrn_md=0,
                     min_obs=100,
                     min_fpkm_uq=0,
                     step_from=4,
                     debug=1)
  
}

message('loading modules -----------------------------------------------------')
library(goldmine)
library(ggbio)
# library(SigFuge)
library(ggforce)
library(Rsamtools)
library(gridExtra)
library(grid)
library(matrixStats)
library(handy)
library(matrixStats)
library(GenomicRanges)
library(parallel)
library(rtracklayer)
#library(canova)
library(reshape2)
library(Gviz)
library(coMET)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(phastCons100way.UCSC.hg19)
library(org.Hs.eg.db)
library(GenomicScores)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
mount_prefix <- get_mount_dir()
source(file.path(Sys.getenv('R_UTIL'),'lib_chipseq.r'))
library(R.utils)
library(openxlsx)

source(file.path(Sys.getenv('R_UTIL'),'lib_tcga.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_PRs.r'))

cpgIslands_UCSC2 <- function(gen, chr, start, end, title="CpG Islands UCSC") 
{
  if (is.null(chr)) {
    stop("Invalid in function cpgIslandsUCSC:chr null:\n")
  }
  if (is.null(start)) {
    stop("Invalid in function cpgIslandsUCSC :start null:\n")
  }
  if (is.null(end)) {
    stop("Invalid in function cpgIslandsUCSC :end null:\n")
  }
  if (is.null(gen)) {
    stop("Invalid in function cpgIslandsUCSC :gen null:\n")
  }
  genTrunk <- gsub("\\..*", "", gen)
  return(UcscTrack(genome = genTrunk, chromosome = chr, track = "cpgIslandExt", 
                   from = start, to = end, trackType = "AnnotationTrack", 
                   start = "chromStart", end = "chromEnd", id = "name", 
                   shape = "box", fill = "#006400", name = title, 
                   stacking = "dense", col.line = NULL, col = NULL))
}

merge_by_gene_and_submitter_id <- function(goi_apa,pau_dt,methyl_dt,debug=F) {
  
  if(debug){browser()}
  
  message('pau_dt last two columns correspond to test_id and gene_name ...')
  
  paur <- pau_dt$ratio
  
  ti2gene_map <- unique(goi$apa[,list(test_id,gene_name)])
  paur$gene_name <- ti2gene_map[match(paur$test_id,ti2gene_map$test_id),gene_name]

  row_ids <- c("test_id","gene_name")
  
  pau_samples <- setdiff(colnames(paur),row_ids)
  ret2 <- melt(paur,id.vars=row_ids,measure.vars = pau_samples)
  ret2 <- setnames(ret2,c('variable','value'),c('cohort_submitter_id','pau'))
  pau_melt <- unique(ret2[(!is.na(ret2$pau) & ret2$pau<=1.0),])
  
  methyl_vs_pau <- merge(methyl_dt, pau_melt, by = c("gene_name","cohort_submitter_id"), all=F, allow.cartesian=TRUE)
  methyl_vs_pau <- methyl_vs_pau[(!is.na(pau) &!is.na(nbeta)),]
  
  return(methyl_vs_pau)
}


masking_low_expr_by_ranking <- function(goi,fpkm0_dt,apa_gene_map,topk=10000,debug=F) {
  
  if(debug){browser()}
  message('use gene expression ranking across samples ....')

  message('we are interested in only query goi among apa_gene ....')
  goi_map <- apa_gene_map[match(unique(goi$apa$gene_name),gene_name),]
  
  message('annotate a ranking # to the goi genes in each sample w.r.t. gexpr.....')
  S <- length(fpkm0_dt)
  samples <- names(fpkm0_dt)
  for (s in 1:S) {
    
    fpkm0_dt[[s]]$goi <- goi_map[match(fpkm0_dt[[s]]$ensembl_gene,goi_map$ensembl_gene),gene_name]
    fpkm0_dt[[s]] <- fpkm0_dt[[s]][order(-expr)] #sort by the gexpr within a sample
    fpkm0_dt[[s]]$expr_ranking <- 1:dim(fpkm0_dt[[s]])[1] #annotate with a ranking
    fpkm0_dt[[s]]$sample <- samples[s]
    fpkm0_dt[[s]] <- fpkm0_dt[[s]][!is.na(goi),]
    fpkm0_dt[[s]][expr_ranking>topk,expr:=0.]
  }
  fpkm1_dt <- rbindlist(fpkm0_dt) #concat all sample tables to row-wide
  rm(fpkm0_dt)
  
  message('dcasting the fpkm matrix and to produce a gene x sample matrix with gexpr. Note that gexpr 0 means is either no expressed or out of topk expressed within a sample')
  
  fpkm2 <- as.data.table(dcast(fpkm1_dt,goi~sample,value.var="expr"))
  rownames(fpkm2) <- fpkm2$goi
  fpkm2$goi <- NULL
  fpkm2 <- sort_column(fpkm2)
  
  paur <- goi$pau$ratio
  # stopifnot(colnames(goi$pau)[1:L2]==colnames(goi$fpkm_prof_mat))
  
  cns <- colnames(paur)
  scols <- cns[cns!='test_id']
  pau_val <- sort_column(paur[,..scols])

  ti2gene_map <- unique(goi$apa[,list(test_id,gene_name)])
  pau_genes <- ti2gene_map[match(paur$test_id,ti2gene_map$test_id),gene_name]

  pau_to_fpkm_ridx <- match(pau_genes,rownames(fpkm2))
  pau_to_fpkm_cidx <- match(colnames(pau_val),colnames(fpkm2))
  
  min_fpkm_uq <- 1. #since we set 0 to the gexpr not in top-k
  
  mask_indicator <- (fpkm2[pau_to_fpkm_ridx,pau_to_fpkm_cidx,with=F] < min_fpkm_uq)
  mask_indicator[is.na(mask_indicator)] <- F
  
  if(debug){browser()}
  pau_val[mask_indicator] <- NA
  goi$pau$ratio <- cbind(pau_val,paur$test_id) #bring back to original format
  colnames(goi$pau$ratio) <- c(colnames(pau_val),'test_id')
  return(goi)
}

masking_low_expr <- function(goi,min_fpkm_uq=1e4,debug=F) {
  
  message('use minimum gene expression ....')
  
  if(debug){browser()}
  
  L <- dim(goi$pau$ratio)[2]#last column is test_id
  # stopifnot(colnames(goi$pau)[1:L1]==colnames(goi$fpkm_prof_mat))
  
  fpkm2 <- goi$fpkm_prof_mat
  fpkm2 <- sort_column(fpkm2)
  
  for (i in 1:length(goi$pau)) { #for each component of mrna (ratio, common, unique)
    paur <- goi$pau[[i]]

    cns <- colnames(paur)
    scols <- cns[cns!='test_id']
    pau_val <- paur[,..scols]

    test_ids <- paur$test_id
    pau_val <- sort_column(pau_val)
    pau_genes <- unlist(goi$apa[match(unlist(test_ids),unlist(goi$apa$test_id)),'gene_name'])
    
    pau_to_fpkm_ridx <- match(pau_genes,rownames(fpkm2))
    pau_to_fpkm_cidx <- match(colnames(pau_val),colnames(fpkm2))
    
    mask_indicator <- fpkm2[pau_to_fpkm_ridx,pau_to_fpkm_cidx] < min_fpkm_uq
    mask_indicator[is.na(mask_indicator)] <- F
    pau_val[mask_indicator] <- NA
    
    goi$pau[[i]] <- cbind(test_ids,pau_val)
  }
  
  return(goi)
}

get_pau_on_goi <- function(rnaseq_dt,goi38_gr,ncpu=1,debug=F,mosdepth_bin="mosdepth"){
  
  # https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
  # we set min mapq to 5. TCGA BAM is created by STAR and the min value will consider only uniquely mapping
  
  #rd_fn_tmp <- paste0(rd_fn,'.tmp')
  
  if (debug) {
    message('lapply *******************************')
    # bam,apa_test_gr,apa_test_bed_fpath,minMapQuality=30, min_depth=0., reuse_bed_fpath=T,ncpu=4,mosdepth_bin="mosdepth", comp_ratio=T,debug=F
    
    goi_pau <- lapply(1:nrow(rnaseq_dt),function(i){
      compPauRatioFromBam_mosdepth(rnaseq_dt$file_path[i][[1]],
                                   goi38_gr,
                                   minMapQuality=5,
                                   ncpu=ncpu,
                                   mosdepth_bin=mosdepth_bin,
                                   debug=debug)
    }
    )
  } else {
    goi_pau <- mclapply(1:nrow(rnaseq_dt),function(i){
      compPauRatioFromBam_mosdepth(rnaseq_dt$file_path[i][[1]],
                                   goi38_gr,
                                   minMapQuality=5,
                                   ncpu=1,
                                   mosdepth_bin=mosdepth_bin,
                                   debug=debug)
    },mc.cores = ncpu, mc.preschedule = FALSE
    )
  }
  
  #save(goi_pau,file=rd_fn_tmp,compress=T)
  
  test_ids <- colnames(goi_pau[[1]])
  row_entries <- rownames(goi_pau[[1]])
  
  pau_list <- list()
  for (j in 1:length(row_entries)) {
    ret <- lapply(1:nrow(rnaseq_dt), function(i){
      return(goi_pau[[i]][j,])
    })
    ret <- as.data.table(ret)
    colnames(ret) <- rnaseq_dt$cohort_submitter_id
    rownames(ret) <- test_ids
    ret <- sort_column(ret)
    ret$test_id <- as.numeric(rownames(ret))
    pau_list[[j]] <- ret
  }
  
  names(pau_list) <- row_entries
  
  #goi_pau <- pau_list$ratio
  
  missing_test_ids <- setdiff(goi38_gr$test_id,pau_list$ratio$test_id)
  message('missing test_ids in PAU analysis ...\n',missing_test_ids)
  
  stopifnot(isEmpty(missing_test_ids))
  
  #goi_pau$gene_name <- goi38_gr$gene_name[match(goi_pau$test_id,goi38_gr$test_id)]
  #pau_list$ratio <- goi_pau
  
  return(pau_list)
}

norm_450k_methyl <- function(methyl_sample_dt,rd_fpath,ncpu=1,debug=F) {
  if (debug){browser()}
  message('normalize beta per each sample/cohort ...')
  methyl_sample_dt$norm_beta_fpath <- run_bmiq_on_cohort(methyl_sample_dt$file_path,ncpu,reuse_bmiq=T,debug=debug)
  return(methyl_sample_dt)
}

build_tcga_data_table <- function (tcga_rd_prefix,cohorts,debug=F){
  if (debug){browser()}
  stopifnot(file.exists(tcga_rd_prefix))
  data_cat_prefix <- c('TCGA-MethArray','TCGA-RNAseq','TCGA-FPKM-UQ')
  cat_fext_pat <- c('gdc_hg38.txt','bam','FPKM-UQ.txt.gz')
  tcga_cat_based <- file.path(tcga_rd_prefix,data_cat_prefix)
  
  sample_dts <- list()
  n <- 0
  for (i in 1:length(tcga_cat_based)){
    catd <- tcga_cat_based[[i]]
    for (j in 1:length(cohorts)) {
      cohortd <- file.path(catd,cohorts[j])
      uu_ids <- list.dirs(cohortd, recursive=FALSE, full.names = F)
      
      sample_fpaths <- list()
      for (k in 1:length(uu_ids)) {
        sampd <- file.path(cohortd,uu_ids[k])
        files <- list.files(sampd,pattern = sprintf('\\.%s$',cat_fext_pat[i]),full.names = T)
        #message(sprintf('sampd:%s, cat_fext_pat:%s',sampd,cat_fext_pat[i]))
        if (isEmpty(files)){
          sample_fpaths[[k]] <- NA
        } else {
          sample_fpaths[[k]] <- files[[1]]
        }
      }
      sample_dt <- data.table(uuid=uu_ids,
                              submitter_id=NA,
                              file_path=unlist(sample_fpaths),
                              cohort=cohorts[j],
                              data_cat=data_cat_prefix[[i]])
      
      sample_dt <- sample_dt[!is.na(file_path),]
      if (dim(sample_dt)[1] > 0){
        n <- n + 1
        sample_dts[[n]] <- sample_dt
      }
    }
  }
  
  sample_dt <- rbindlist(sample_dts)
  
  message('retrieving submitter_ids ------------------------------------------')
  sample_dt$submitter_id <- map_uuid_to_submitter_id(sample_dt$uuid)
  
  message('combining submitter_id with cohorts -------------------------------')
  sample_dt$cohort_submitter_id <- paste0(sample_dt$cohort,'_',sample_dt$submitter_id)
  return(sample_dt)
}

liftover_hg19_to_hg38 <- function(query_hg19_gr,debug=F){
  if (debug){browser()}
  loChain <- list()
  loChain$b37_b38 <- import.chain(sprintf('%s/projects/data/ref/liftOverChain/hg19ToHg38.over.chain',mount_prefix))
  
  query38 <- liftOver(query_hg19_gr,loChain$b37_b38)
  query38_dt <- as.data.table(query38)
  missing_genes <- setdiff(query_hg19_gr$gene_name,query38_dt$gene_name)
  message('missing genes during hg38 conversion process ...\n',missing_genes)
  stopifnot(isEmpty(missing_genes))
  
  query38_dt <- query38_dt[,list(seqnames,start,end,strand,gene_name,type,test_id)]
  setnames(query38_dt,old=c('seqnames'),new=c('chr'))
  query38_gr <- makeGRanges(query38_dt)
  
  return(query38_gr)
}

compute_correlation <- function(methyl_pau,ncpu=1,min_obs=100,debug=F){
  if(debug){browser()}
  test_id <- methyl_pau[1,test_id]
  cohort <- methyl_pau[1,cohort]
  
  methyl_pau_by_pos <- split(methyl_pau, by=c('bp_pos'), sorted=TRUE, drop=TRUE)
  
  N <- length(methyl_pau_by_pos)
  
  message(sprintf('computing r2 between pau and beta, cohort[%s], test_id[%d], N[%d].....',cohort,test_id,N))
  
  ret_dt0 <- data.table(r2p=0.,
                        r2p_pval=1.,
                        # r2s=0.,
                        # r2s_pval=1.,
                        # r2k=0.,
                        # r2k_pval=1.,
                        # R2=1.,
                        # R2_pval=1.,
                        data_points=0)
  
  corrs <- lapply(1:N,function(i){
    mp_by_pos <- methyl_pau_by_pos[[i]]
    
    ret_dt <- ret_dt0
    
    if (dim(mp_by_pos)[1] >= min_obs & 
        length(unique(mp_by_pos$nbeta))>2 & 
        length(unique(mp_by_pos$pau))>2) {
      
      p_corr <-cor.test(mp_by_pos$nbeta,
                        mp_by_pos$pau,
                        method='pearson')
      
      # s_corr <-cor.test(mp_by_pos$nbeta,
      #                     mp_by_pos$pau,
      #                     method='spearman')
      #       
      # k_corr <-cor.test(mp_by_pos$nbeta,
      #                     mp_by_pos$pau,
      #                     method='kendall')
      # 
      # regfit_result <- regression_fit(mp_by_pos,
      #                                 method=fit_model,
      #                                 debug=debug)
      
      ret_dt <- data.table(r2p=p_corr$estimate,
                           r2p_pval=p_corr$p.value,
                           # 
                           # r2s=s_corr$estimate,
                           # r2s_pval=s_corr$p.value,
                           # 
                           # r2k=k_corr$estimate,
                           # r2k_pval=k_corr$p.value,
                           # 
                           # R2=regfit_result$rsq,
                           # R2_pval=regfit_result$pval,
                           data_points=dim(mp_by_pos)[1])
    }
    return(ret_dt)
  }
  )
  
  corrs <- rbindlist(corrs)
  corrs$bp_ps <- names(methyl_pau_by_pos)
  colnames(corrs) <- c('pearson.r2','pearson.pval',
                       # 'spearman.r2','spearman.pval',
                       # 'kendall.r2','kendall.pval',
                       # 'R2','R2.pval',
                       'data_points','probe_pos')
  
  smallest_pval <- min(corrs[pearson.pval>0.,pearson.pval])
  if (smallest_pval < 1e-30) {
    smallest_pval <- 1e-30
  }
  corrs[pearson.pval==0.,pearson.pval:=smallest_pval]
  
  return(corrs)
}


methyl_vs_pau_scatter <- function(methyl_pau_dt,tu_gr_ext,out_dir,cohort,method_prefix,ncpu=1,min_obs=100,reuse_pdf=F,debug=F){
  
  if (debug){browser()}
  
  message('remove all cases with pau 1.0 which unique md > common md ...')
  methyl_pau_dt <- methyl_pau_dt[pau<=1.0,]
  
  methyl_pau_dt_by_test <- split(methyl_pau_dt, by=c('test_id'), sorted=TRUE, drop=TRUE)
  test_ids <- names(methyl_pau_dt_by_test)
  
  corr_mpt2 <- list()
  mpt2s <- list()
  
  for (i in 1:length(methyl_pau_dt_by_test)) {
    if(debug){message(i)}
    
    mpt2 <-methyl_pau_dt_by_test[[i]]
    
    test_id <- test_ids[[i]]
    
    # debug <- T
    # browser()
    
    gene_name <- mpt2$gene_name[1]
    
    mpt2$bp_pos <- paste0(mpt2$chr,':',mpt2$start)
    mpt2 <- mpt2[order(bp_pos)]
    
    corr_mpt2[[i]] <- compute_correlation(mpt2,
                                          ncpu=ncpu,
                                          #fit_model=fit_model,
                                          min_obs=min_obs,
                                          debug=debug)
    
    corr_mpt2[[i]]$test_id <- test_id
    corr_mpt2[[i]]$gene_name <- gene_name
    
    mpt2s[[i]] <- mpt2
  }
  
  
  names(mpt2s) <- names(methyl_pau_dt_by_test)
  methyl_pau_dt <- rbindlist(mpt2s)
  
  names(corr_mpt2) <- names(methyl_pau_dt_by_test)
  corr_mpt2 <- rbindlist(corr_mpt2)
  rm(mpt2s)
  message('computing adjust p-value using BH ...')
  corr_mpt2_by_gene_name <- split(corr_mpt2, by=c('gene_name'), sorted=TRUE, drop=TRUE)
  
  corr_mpt3 <- list()
  for (i in 1:length(corr_mpt2_by_gene_name)) {
    mpt2 <-corr_mpt2_by_gene_name[[i]]
    
    mpt2$pearson.pval <- p.adjust(mpt2$pearson.pval, method = "BH")
    #mpt2$spearman.pval <- p.adjust(mpt2$spearman.pval, method = "BH")
    #mpt2$kendall.pval <- p.adjust(mpt2$kendall.pval, method = "BH")
    #mpt2$R2.pval <- p.adjust(mpt2$R2.pval, method = "BH")
    
    corr_mpt3[[i]] <- mpt2
  }
  
  corr_mpt2 <- rbindlist(corr_mpt3)
  rm(corr_mpt3)
  
  message('printing scatter plot ...')
  
  #corr_mpt2 <- corr_mpt2[(pearson.pval<0.05 & abs(pearson.r2)>0.1),]
  
  corr_mpt2_by_test <- split(corr_mpt2, by=c('test_id'), sorted=TRUE, drop=TRUE)
  methyl_pau_dt_by_test <- split(methyl_pau_dt, by=c('test_id'), sorted=TRUE, drop=TRUE)
  
  test_ids <- names(methyl_pau_dt_by_test)
  
  for (i in 1:length(methyl_pau_dt_by_test)) {
    
    mpt2 <-methyl_pau_dt_by_test[[i]]
    test_id <- test_ids[[i]]
    
    corr_mpt2 <- corr_mpt2_by_test[[test_id]]
    
    mpt2$data_points<-corr_mpt2[match(mpt2$bp_pos,corr_mpt2$probe_pos),data_points]
    mpt2 <- mpt2[data_points>min_obs,]
    
    if (dim(mpt2)[1]>0) {
      message(sprintf('filtering by data_points, cohort[%s], test_id[%s]...',cohort,test_id))
      
      gene_name <- mpt2$gene_name[1]
      tu <- tu_gr_ext[tu_gr_ext$name==gene_name,]
      chrom <- as.character(seqnames(tu))
      
      if(debug){message(i)}
      
      pdf_file <- file.path(out_dir,sprintf('%s_%s_%s_corrs.pdf',cohort,test_id,gene_name))
      
      if (reuse_pdf & file.exists(pdf_file) & (file.size(pdf_file)>5e3)){
        message(sprintf('reuse prev pdf [%s]',pdf_file))
      } else {
        if (F) {
          labels2 <- sprintf('%s\ndata_pts:%d\nr(%3.2e,%3.2e)\nR(%3.2e,%3.2e)',
                             corr_mpt2$probe_pos,
                             corr_mpt2$data_points,
                             corr_mpt2$kendall.r2,
                             corr_mpt2$kendall.pval,
                             corr_mpt2$R2,
                             corr_mpt2$R2.pval)
          names(labels2) <- corr_mpt2$probe_pos
        } else {
          mpt2s <- split(mpt2,by=c('bp_pos'), sorted=FALSE,drop=FALSE)
          labels <- lapply(mpt2s,function(x){
            
            label <- sprintf('%s\ndata_pts:%d',
                             x$bp_pos[1],
                             x$data_points[1])
            return(label)
            
          })
          
          labels2 <- unlist(labels)
        }
        
        ncol2 <- 4
        nrow2 <- 3
        n_pages <- ceiling(
          length(unique(mpt2$bp_pos)) / (ncol2 * nrow2)
        )
        
        if(!debug){pdf(pdf_file)}
        
        p0 <- ggplot(mpt2,aes(x=nbeta,y=pau,color=cohort)) +
          geom_point(alpha=0.5,size=1,shape=20)
        
        for (j in 1:n_pages) {
          p <- p0 + 
            facet_wrap_paginate( ~bp_pos, ncol=ncol2, nrow=nrow2, page=j, scales = "fixed", labeller=labeller(bp_pos=labels2)) + 
            ggtitle(sprintf('(%s)%s[%s:%s-%s], %s, #%d/%d',
                            test_id,gene_name, chrom,comma1k(start(tu)),comma1k(end(tu)),method_prefix,j,n_pages)) +
            xlab(expression("DNA methylation level ("~beta~") value")) +
            ylab('predicted APA isoform usage ratio') +
            scale_x_continuous(breaks=c(0., 0.5, 1.0),limits = c(0.,1.)) +
            scale_y_continuous(breaks=c(0., 0.5, 1.0),limits = c(0.,1.)) +
            theme_bw() +
            # guides(color = guide_legend(override.aes = list(size=6))) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
          
          plot(p)
        }
        if(!debug) {dev.off()}
      }
    }
  }
  
  return(corr_mpt2_by_test)
}


methyl_vs_pau_gtrack <- function(corr_mpts,goi,out_dir,chipseq_dt,ncpu=1,cohort_str='NA',method_prefix,reuse_pdf=F,min_obs=100,debug=F) {
  if (debug){browser()}
  
  num_of_test <-length(corr_mpts) #number of test
  test_ids <- as.numeric(names(corr_mpts))
  min_r2 <- -1.05
  max_r2 <- 1.05
  max_r2_abs <- 1.05
  ugrid <- 50
  
  r2_cutoff <- 0.10
  pval_cutoff <- 0.05
  
  message("loading phastCons data ...")
  #phast <- getGScores("phastCons100way.UCSC.hg19")
  phast <- phastCons100way.UCSC.hg19
  
  message("loading UCSC gene model ...")
  txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
  chroms <- as.character(unique(seqnames(goi$tu_gr_ext)))
  
  i.trs <- lapply(1:length(chroms),function(i){
    message(sprintf('loading IdeogramTrack [%s]...',chroms[i]))
    IdeogramTrack(genome="hg19", chromosome=chroms[i])
  })
  names(i.trs) <- chroms
  
  g.trs <- lapply(1:length(chroms),function(i){
    message(sprintf('loading GenomeTrack [%s]...',chroms[i]))
    GenomeAxisTrack(genome="hg19", chromosome=chroms[i], littleTicks = TRUE)
  })
  names(g.trs) <- chroms
  
  lapply(1:num_of_test,function(i) {
    
    message(sprintf('**** [%d/%d] ****',i,num_of_test))
    corr_dt <- corr_mpts[[i]]
    corr_dt[is.na(pearson.r2),pearson.r2:=0.]
    corr_dt[is.na(pearson.pval),pearson.pval:=1.]
    
    any_signifi_corr <- corr_dt[data_points>min_obs & abs(pearson.r2)>r2_cutoff & pearson.pval<pval_cutoff,]
    test_id2 <- test_ids[[i]]
    
    apa <- goi$apa[test_id==test_id2,]
    chrom <- apa$chr[1]
    gene_name <- apa$gene_name[1]
    
    tu <- goi$tu_gr_ext[goi$tu_gr_ext$name==gene_name,]
    from1 <- start(tu)
    to2 <- end(tu)
    
    if (width(tu)>1e6) {
      
      corr_peak_pos <- as.numeric(tstrsplit(corr_dt[abs(pearson.r2) == max(abs(corr_dt$pearson.r2)),'probe_pos'][1],':')[[2]])
      bnd <- c(corr_peak_pos,apa$start,apa$end)
      
      from1 <- min(bnd) - 5e3
      if (from1 < start(tu)) {
        from1 <- start(tu)
      }
      
      to2 <- max(bnd) + 5e3
      if (end(tu) < to2) {
        to2 <- end(tu)
      }
    }
    
    if (dim(any_signifi_corr)[1]>0 & (to2 - from1)<1e6) {
      message(sprintf('proceeding_analysis on test_id [%d]',test_id2))
      
      pdf_file <- file.path(out_dir,sprintf('%s_%s_%s_gtracks.pdf',cohort_str,test_id2,gene_name))
      
      if ((reuse_pdf & file.exists(pdf_file) & (file.size(pdf_file)>5e3))) {
        message(sprintf('reuse prev result [%s]',pdf_file))
      } else {
        message(sprintf('printing test_id[%d],gene_name[%s],pdf_file[%s]',test_id2,gene_name,pdf_file))
        
        if(!debug){pdf(pdf_file,width=7,height=9.06)}
        
        #-----------------------
        apa$color2 <- as.character()
        apa[type=='unique',color2:='blue']
        apa[type=='shared',color2:='red']
        apa_gr <- makeGRanges(apa,strand = F)
        
        expr.tr <- AnnotationTrack(range=apa_gr,
                                   group=apa_gr$type,
                                   genome='hg19',
                                   name='TR',
                                   fill=apa_gr$color2)
        
        #-----------------------
        i.tr <- i.trs[chrom][[1]]
        
        #-----------------------
        g.tr <- g.trs[chrom][[1]]
        
        #-----------------------
        ij <- findOverlaps(goi$pr,tu,ignore.strand=T)
        
        pr.gr <- unique(goi$pr[from(ij),])
        pr.tr <- AnnotationTrack(start=start(pr.gr),
                                 width=width(pr.gr),
                                 chromosome=chrom,
                                 strand='*',
                                 genome='hg19',
                                 name='PR')
        
        #-----------------------
        grid.gr <- GRanges(seqnames=chrom, 
                           IRanges(start=seq(from1,to2-ugrid,by=ugrid),
                                   width=ugrid))
        
        cons.tr <- DataTrack(range=scores(phast, grid.gr),
                             genome="hg19",
                             type="s",
                             chromosome = chrom,
                             name="PhCon")
        
        #-----------------------
        cg.tr <- cpgIslands_UCSC2("hg19", chrom, from1, to2, title="C")
        
        # ----------------------
        gr.tr <- GeneRegionTrack(txdb_hg19,
                                 genome="hg19", 
                                 chromosome=chrom, 
                                 name="Genes",
                                 stacking="hide",
                                 start=from1,
                                 end=to2,
                                 thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5", "miRNA", "lincRNA"))
        
        # -----------------------
        # corr_dt <- corr_mpts[[i]]
        # corr_dt[is.na(kendall.r2),kendall.r2:=0.]
        # corr_dt[is.na(kendall.pval),kendall.pval:=1.]
        corr_dt[, c("chromosome", "start") := tstrsplit(probe_pos, ":", fixed=TRUE)]
        corr_dt$start <- as.numeric(corr_dt$start)
        corr_dt$end <- corr_dt$start + 1
        
        #handling NA
        corr_dt[is.na(pearson.r2),pearson.r2:=0.]
        corr_dt[is.na(pearson.pval),pearson.pval:=1.]
        
        #handling 0 pvalues
        min_pval <- min(corr_dt[pearson.pval>0.,pearson.pval])
        if (min_pval < 1e-30){
          corr_dt[pearson.pval==0.,pearson.pval := min_pval*0.1]
        } else {
          corr_dt[pearson.pval==0.,pearson.pval := 1e-30]
        }
        
        #-log10(pval)
        corr_dt$pval_nlog10 <- -log10(corr_dt$pearson.pval)
        
        # handling unusual -log10(pval)
        corr_dt[is.na(pval_nlog10),pval_nlog10:=0.]
        corr_dt[is.infinite(pval_nlog10),pval_nlog10:=30.]
        
        pval_nlog10_max <- max(corr_dt$pval_nlog10)
        if (pval_nlog10_max == 0.) {
          pval_nlog10_max <- 1.
        }
        
        corr_dt <- setnames(corr_dt,'pearson.r2','value')
        
        corr.tr <- DataTrack(range = corr_dt[,list(chromosome,start,end,value)],
                             type = c('h'),
                             genome = 'hg19',
                             col = "black",
                             ylim = c(-max_r2_abs,max_r2_abs),
                             name = 'r2')
        
        corr_pval.tr <- DataTrack(range = corr_dt[,list(chromosome,start,end,pval_nlog10)],
                                  type = c('h'),
                                  genome = 'hg19',
                                  col = "red",
                                  ylim = c(0,pval_nlog10_max),
                                  name = '-log10(p)')
        
        ali.trs <- list()
        #https://support.bioconductor.org/p/59280/
        
        for (j in 1:dim(chipseq_dt)[1]){
          ali.trs[[j]] <- AlignmentsTrack(chipseq_dt$fpath[j],
                                          isPaired = FALSE, 
                                          name = chipseq_dt$sample[j],
                                          genome = 'hg19',
                                          chromosome=chrom,
                                          type=c("coverage"),
                                          start=from1,
                                          end=to2,
                                          size=0.03,
                                          fill = chipseq_dt$color2[j],
                                          id=NULL,
                                          cigar=NULL,
                                          mapq=0,
                                          seqs=NULL,
                                          referenceSequence = NULL)
        }
        
        # ----------------------
        tr1s <- list(i.tr, g.tr, corr.tr, corr_pval.tr, pr.tr, cons.tr, expr.tr, cg.tr, gr.tr)
        tr2s <- c(tr1s,ali.trs)
        
        size_trs <- c(2,     4,      4,       3,          1,       2,      2,      1,     2,   2,2,2,2,2,2,2,2)
        
        title2 <- sprintf('(%s)%s[%s:%s-%s],cohort:%s,%s',test_id2,gene_name, chrom,comma1k(from1),comma1k(to2),cohort,method_prefix)
        
        message(sprintf('printing gtrack [%s] ...',title2))
        
        rm(tr1s)
        rm(i.tr)
        rm(corr.tr)
        rm(corr_pval.tr)
        rm(pr.tr)
        rm(cons.tr)
        rm(expr.tr)
        rm(cg.tr)
        rm(gr.tr)
        rm(ali.trs)
        
        plotTracks(tr2s,
                   chromosome = chrom,
                   from = from1,
                   to = to2,
                   collapse=FALSE,
                   main=title2,
                   cex.main=1,
                   sizes=size_trs)
        
        rm(g.tr)
        rm(pr.gr)
        rm(tr2s)
        
        if(!debug){dev.off()}
        
      }
    }
    else {
      message(sprintf('skipping_analysis on test_id [%d]',test_id2))
    }
  })
}

synchronize_samples <- function(sample_dt) {
  
  sample_dts <- split(sample_dt,by=c("cohort","data_cat"),sorted=T,drop=T)
  
  sample_dts2 <- lapply(sample_dts, function(dt) {
    dt <- dt[!duplicated(dt$cohort_submitter_id),]
    dt <- dt[order(cohort_submitter_id)]
    return(dt)
  })
  
  #450k Methyl and FPKM should be matched !
  #FPKM and mRNA should be matched !
  cohorts <- unique(sample_dt$cohort)
  for (i in 1:length(cohorts)) {
    meth_name <- sprintf('%s.TCGA-MethArray',cohorts[i])
    fpkm_name <- sprintf('%s.TCGA-FPKM-UQ',cohorts[i])
    mcsi <- sample_dts2[[meth_name]]$cohort_submitter_id
    fcsi <- sample_dts2[[fpkm_name]]$cohort_submitter_id
    comm_samples <- intersect(mcsi,fcsi)
    
    sample_dts2[[meth_name]] <- sample_dts2[[meth_name]][cohort_submitter_id %in% comm_samples,]
    sample_dts2[[fpkm_name]] <- sample_dts2[[fpkm_name]][cohort_submitter_id %in% comm_samples,]
  }
  
  #FPKM and mRNA should be matched !
  for (i in 1:length(cohorts)) {
    fpkm_name <- sprintf('%s.TCGA-FPKM-UQ',cohorts[i])
    mrna_name <- sprintf('%s.TCGA-RNAseq',cohorts[i])
    fcsi <- sample_dts2[[fpkm_name]]$cohort_submitter_id
    mcsi <- sample_dts2[[mrna_name]]$cohort_submitter_id
    comm_samples <- intersect(fcsi,mcsi)
    
    sample_dts2[[fpkm_name]] <- sample_dts2[[fpkm_name]][cohort_submitter_id %in% comm_samples,]
    sample_dts2[[mrna_name]] <- sample_dts2[[mrna_name]][cohort_submitter_id %in% comm_samples,]
  }
  
  sample_dt <- rbindlist(sample_dts2)
  
  methyl_dt <- sample_dt[data_cat=='TCGA-MethArray',][order(cohort_submitter_id)]
  fpkm_dt <- sample_dt[data_cat=='TCGA-FPKM-UQ',][order(cohort_submitter_id)]
  rnaseq_dt <- sample_dt[data_cat=='TCGA-RNAseq',][order(cohort_submitter_id)]
  rm(sample_dts)
  
  return(tcga_fpath=list(all=sample_dt,
                         methyl_dt=methyl_dt,
                         fpkm_dt=fpkm_dt,
                         rnaseq_dt=rnaseq_dt))
}

message('main() ==============================================================')
prog_lab <- '03_tcga_analy_apa'

message(sprintf('resume from the step[%d]/10',args$step_from))

step_from <- as.integer(args$step_from)

out_dir<-args$out_dir
if ( !dir.exists(out_dir) ) {
  dir.create(out_dir)
} 

message(sprintf('output base directory:%s',out_dir))
method_prefix <- args$method
print_outdir <- file.path(out_dir,method_prefix)
if ( !dir.exists(print_outdir) ) {
  dir.create(print_outdir)
}

message(sprintf('print output directory:%s',print_outdir))

debug2 <- F
if (args$debug==1) {debug2<-T}

stepi <- 0
rd_fn <- file.path(out_dir,'apa_gene_map.rd')
message(rd_fn)
stepi <- stepi + 1
if ((step_from>stepi) & (file.exists(rd_fn))) {
  load(rd_fn)
} else {
  message('1. attaching gene locus (TU) of goi ----------------------------------')
  
  apa_gene_map <- map_refgene_to_ensembl(args$prepropa_rd,args$apa_ann_rd,debug=debug2)
  save(apa_gene_map,file=rd_fn,compress =T)
}

message('2. reading methyl/PAU rd files in multiple TCGA cohorts--------------')
cohorts <- unlist(strsplit(args$cohorts,','))

rd_fpath <- file.path(out_dir,'tcga_methyl_pau.rd')

stepi <- stepi + 1
if ((step_from>stepi) & (file.exists(rd_fpath))) {
  load(rd_fpath)
} else {
  message('loading cohort sample rd files ...')
  tcga_fpaths <- lapply(cohorts, function(cohort) {
    rd_fn <- file.path(sprintf('%s_%s',args$tcga_rd_prefix,cohort),'sample_dt.rd')
    message(rd_fn)
    load(rd_fn)
    return(tcga_fpath)
  })
  
  tcga_fpath<-list()
  
  tcga_fpath$all <- rbindlist(lapply(tcga_fpaths,function(tcga_fpath){return(tcga_fpath$all)}))
  tcga_fpath$tcga_methyl <- rbindlist(lapply(tcga_fpaths,function(tcga_fpath){return(tcga_fpath$methyl_dt)}))
  tcga_fpath$fpkm <- rbindlist(lapply(tcga_fpaths,function(tcga_fpath){return(tcga_fpath$fpkm_dt)}))
  tcga_fpath$rnaseq <- rbindlist(lapply(tcga_fpaths,function(tcga_fpath){return(tcga_fpath$rnaseq_dt)}))
  rm(tcga_fpaths)
  
  message('loading mRNA-seq data --------------------------------------------')
  
  gois <- lapply(cohorts, function(cohort) {
    rd_fn <- file.path(sprintf('%s_%s',args$tcga_rd_prefix,cohort),'goi_rnaseq.rd')
    message(sprintf('loading rd_fn[%s]',cohort))
    load(rd_fn)
    return(goi)
  })
  
  goi <- list()
  goi$apa <- gois[[1]]$apa
  goi$apa_gr <- gois[[1]]$apa_gr
  goi$tu_gr <- gois[[1]]$tu_gr
  goi$tu_gr_ext <- gois[[1]]$tu_gr_ext
  goi$pr <- gois[[1]]$pr
  goi$apa38_gr <- gois[[1]]$apa38_gr
  
  goi$methyl_prof_dt <- rbindlist(lapply(gois,function(goi2){return(goi2$methyl_prof_dt)}))
  
  pau_all <- list()
  for (i in 1:length(gois)) {
    
    pau <- gois[[i]]$pau
    
    if (i == 1) {
      for (j in 1:length(pau)) {
        pau_all[[j]] <- pau[[j]][order(test_id)]
      }
    } else {
      for (j in 1:length(pau)) {
        pau_all[[j]] <- merge(pau_all[[j]], 
                              pau[[j]],
                              all.x=T,
                              by="test_id")
        
      }
    }
  }
  names(pau_all) <- names(pau)
  goi$pau <- pau_all
  rm(gois)
  rm(pau_all)
  save(goi,tcga_fpath,file=rd_fpath,compress=T)
}

message('3. loading FPKM-UQ and extracting goi genes -------------------------')
rd_fn <- file.path(out_dir,'fpkm.rd')
message(rd_fn)
stepi <- stepi + 1
if ((step_from>stepi) & (file.exists(rd_fn))) {
  load(rd_fn) #fpkm0_dt
} else {
  
  # to check overall fpkm profile across all the genes
  fpkm0_dt <- load_fpkm_rnaseq(tcga_fpath$fpkm,
                               ncpu=args$ncpu,
                               debug=debug2)
  
  fpkm_dt <- load_fpkm_rnaseq_on_goi(tcga_fpath$fpkm,
                                     ncpu=args$ncpu,
                                     debug=debug2)
  
  goi$fpkm_prof_mat <- fpkm_dt[match(goi$tu_gr$name,rownames(fpkm_dt)),]
  goi$fpkm_prof_mat <- sort_column(goi$fpkm_prof_mat)
  save(fpkm0_dt,goi,file=rd_fn,compress = T)
}

message('4. masking genes not confidently expressed --------------------------')
rd_fn <- file.path(out_dir,sprintf('masked_%s.rd',method_prefix))
stepi <- stepi + 1
if (args$expr_topk > 0) {
  message(rd_fn)
  message(sprintf('top-K [%d] is applied ...',args$expr_topk))
  if ((step_from>stepi) & (file.exists(rd_fn))) {
    load(rd_fn)
  } else {
    #update goi$pau$ratio
    goi <- masking_low_expr_by_ranking(goi,fpkm0_dt,apa_gene_map,topk=args$expr_topk,debug=debug2)
    save(goi,file=rd_fn,compress = T)
  }
} else {
  
  message(rd_fn)
  message(sprintf('minimum FPKM-UQ [%d] is applied ...',args$min_fpkm_uq))
  if ((step_from>stepi) & (file.exists(rd_fn))) {
    load(rd_fn)
  } else {
    goi <- masking_low_expr(goi,min_fpkm_uq=args$min_fpkm_uq,debug=debug2)
    save(goi,file=rd_fn,compress = T)
  }
}
rm(fpkm0_dt)

message('5. combining pau to beta --------------------------------------------')
rd_fn <- file.path(out_dir,sprintf('convol_%s.rd',method_prefix))
message(rd_fn)
stepi <- stepi + 1
if ((step_from>stepi) & (file.exists(rd_fn))) {
  load(rd_fn)
} else {
  methyl_vs_pau_list <- list()
  
  message('generate merged cohorts (pau vs. apa) ...')
  
  methyl_vs_pau_list$merged <- merge_by_gene_and_submitter_id(goi$apa,
                                                              goi$pau,
                                                              goi$methyl_prof_dt,
                                                              debug=debug2)
  
  C <- length(cohorts)
  if (C > 1) {
    message('appending each cohort (pau vs. apa) ...')
    ret <- split(methyl_vs_pau_list$merged,by=c("cohort"))
    methyl_vs_pau_list <- append(methyl_vs_pau_list,ret)
  }
  save(methyl_vs_pau_list, file=rd_fn, compress = T)
}
goi$pau <- NULL
goi$methyl_prof_dt <- NULL

hostname <- Sys.getenv('HOSTNAME')

message('6. generating plots from [%s]-------------------------',hostname)

stepi <- stepi + 1

cohorts <- names(methyl_vs_pau_list)
C <- length(cohorts)

chipseq_dt <- get_chipseq_bam(mount_prefix)

rd_fn <- file.path(out_dir,sprintf('corr_%s.rd',method_prefix))

corr_mpts <- list()
if ((step_from>stepi) & (file.exists(rd_fn))) {
  load(rd_fn)
} else {
  for (i in 1:C) {
    cohort <- cohorts[i]
    
    methyl_vs_pau <- methyl_vs_pau_list[[cohort]]
    
    message(sprintf('6a. scatter plot between PAU and methylation rate [%s]',cohort))
    
    corr_mpts[[i]] <- methyl_vs_pau_scatter(methyl_vs_pau,
                                            goi$tu_gr_ext,
                                            print_outdir,
                                            cohort,
                                            method_prefix,
                                            #fit_model=args$fit_model,
                                            ncpu=args$ncpu,
                                            min_obs=args$min_obs,
                                            reuse_pdf=T,
                                            debug=debug2)
    
  }
  names(corr_mpts) <- cohorts
  save(corr_mpts,file=rd_fn,compress = T)
  rm(methyl_vs_pau)
  rm(methyl_vs_pau_list)
}

message('7. plotting gtrack() from [%s]-------------------------',hostname)

stepi <- stepi + 1
if (stepi>=step_from) {
  for (i in 1:C) {
    cohort <- cohorts[i]
    message(sprintf('7a. genome tracks between PAU and methylation rate [%s]',cohort))
    methyl_vs_pau_gtrack(corr_mpts[[i]], 
                         goi, 
                         print_outdir, 
                         chipseq_dt,
                         ncpu=1, 
                         cohort,
                         method_prefix,
                         reuse_pdf=T,
                         min_obs=args$min_obs,
                         debug=debug2)
  }
}
