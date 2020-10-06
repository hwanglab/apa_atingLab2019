library(data.table)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(RPMM)
library(wateRmelon)
library(parallel)
source(file.path(Sys.getenv('R_UTIL_APA'),'lib_util.R'))

methyl450k_extract <- function(methyl_dt,tu_gr_ext,ncpu=1,debug=F){
  #Note that we use hg19 coordinates in 450k array in this function
  if (debug){browser()}
  
  mprobe_in_goi <- methyl450k_extract_beta_by_gr(tu_gr_ext,debug)
  
  methyl_prof_gois <- mclapply(1:nrow(methyl_dt),function(i){
    fpath <- methyl_dt$norm_beta_fpath[i][[1]]
    if (debug){browser()}
    if (file.exists(fpath)){
      msample_dt <- fread(fpath)
      msample_dt <- msample_dt[!is.na(nbeta),]
      sg_idx <- match(msample_dt$probe_id,mprobe_in_goi$probe_id)
      matched_sidx <- which(!is.na(sg_idx))
      matched_gidx <- sg_idx[matched_sidx]
      
      matched_probe <- mprobe_in_goi[matched_gidx,]
      matched_probe$nbeta <- msample_dt[matched_sidx,nbeta]
      methyl_prof_goi <- as.data.table(matched_probe)
      methyl_prof_goi$cohort <- methyl_dt$cohort[i]
      methyl_prof_goi$cohort_submitter_id <- methyl_dt$cohort_submitter_id[i]
      methyl_prof_goi$width <- NULL
      setnames(methyl_prof_goi, "seqnames", "chr")
      
      return(methyl_prof_goi)
    } else {
      stop(sprintf('sample 450k methylation file [%s] does not exist!',fpath))
    }
  },mc.cores = ncpu)
  return(rbindlist(methyl_prof_gois))
}

methyl450k_extract_beta_by_gr <- function(tu_gr,debug){
  if (debug){browser()}
  message('Retrieving 450k array hg19 coordinates and all other info----------')
  methyl_probe_gr <- GRanges(seqnames = Locations$chr,
                             IRanges(start = Locations$pos,
                                     end = Locations$pos + 1),
                             strand = Locations$strand,
                             probe_id = rownames(Locations))
  
  
  methyl_probe_gr <- sortSeqlevels(methyl_probe_gr)
  methyl_probe_gr <- sort(methyl_probe_gr)
  
  ij <- findOverlaps(tu_gr,methyl_probe_gr,ignore.strand=TRUE)
  
  methyl_probe_goi_gr <-methyl_probe_gr[to(ij),]
  
  methyl_probe_goi_gr$gene_name <- tu_gr$name[from(ij)]
  
  missing_probe_genes <- setdiff(tu_gr$name,tu_gr$name[from(ij)])
  message('the following genes do not have any methylation probes!\n',missing_probe_genes)

  return(methyl_probe_goi_gr)
}

# run_BMIQ <- function(Beta_value,probe_idx,input_fname){
#   out <- tryCatch(
#     {
#       BMIQ(Beta_value,probe_idx,plots = F)
#     },
#     error=function(cond) {
#       message(cond)
#       D <- as.integer(length(Beta_value) * 0.06)
#       ret <- BMIQ(Beta_value,probe_idx,plots = F,nfit=D)
#       return(ret)
#     },
#     warning=function(cond) {
#       message(cond)
#       D <- as.integer(length(Beta_value) * 0.06)
#       ret <- BMIQ(Beta_value,probe_idx,plots = F,nfit=D)
#       return(ret)
#     },
#     finally={
#       message(paste("Processed File:", input_fname))
#     }
#   )    
#   return(out)
# }

run_bmiq_on_cohort <- function(methy_450k_files,ncpu,reuse_bmiq=T,debug=F) {
  if(debug){browser()} #debug
  message(sprintf('normalizing beta values with BMIQ method on total [%d] files...',length(methy_450k_files)))
  nfit_default <- 10000
  nprobes <- 350000

  norm_beta_fpaths <- mclapply(methy_450k_files,function(methy_450k_file) {
    
    message(sprintf("Normalizing methyl_450k_file[%s] and writing the result ...",methy_450k_file))

    if (file.exists(methy_450k_file)) {
      outd <- file.path(dirname(methy_450k_file),'normalized')
      norm_beta_fn <- file.path(outd,'beta.tsv')
      if (reuse_bmiq & (file.exists(norm_beta_fn))) {
        message(sprintf('norm_beta_fn[%s] reused ',norm_beta_fn))
      } else {
        dir.create(file.path(outd), showWarnings = FALSE)
        
        samp450k <- fread(methy_450k_file,select=c('Composite Element REF','Beta_value'),header=T)
        setnames(samp450k,"Composite Element REF", "probe_id")
        samp450k <- samp450k[!is.na(Beta_value),]
  
        samp450k$probe_type <- Manifest[match(samp450k$probe_id,Manifest$Name),'Type']
        samp450k <- samp450k[!is.na(probe_type),]
        samp450k$probe_idx <- as.integer()
        samp450k[probe_type=='I',probe_idx:=1]
        samp450k[probe_type=='II',probe_idx:=2]
        
        D <- dim(samp450k)[1]
        nfit2 <- nfit_default
        if (D < nprobes){
          nfit2 <- D * 0.05
        }
        ret <- BMIQ(samp450k$Beta_value,samp450k$probe_idx,nfit=nfit2,plots=FALSE)
        samp450k$probe_type <- NULL
        samp450k$nbeta <- as.numeric()
        samp450k$nbeta <- ret$nbeta
        samp450k$probe_idx <- NULL
        fwrite(samp450k,sep='\t',quote=F, file=norm_beta_fn)
        message(sprintf('norm_beta_fn[%s] successfully written',norm_beta_fn))
      }
    } else {
      norm_beta_fn <- NA
      stop(sprintf('file[%s] does not exist!',methy_450k_file))
    }
    return(norm_beta_fn)
  },mc.cores = ncpu)

  return(norm_beta_fpaths)

}

load_fpkm_rnaseq <- function(fpkm_dt,ncpu=1,debug=F) {
  
  if (debug){browser()}
  
  
  source(file.path(Sys.getenv('R_UTIL_APA'),'lib_PRs.r'))
  
  message('load read count matrix ...')
  exprs <- mclapply(1:nrow(fpkm_dt),function(i){
    sed_cmd <- sprintf("zcat %s | sed \'s/\\.[0-9]\\+//1\'",fpkm_dt[i,file_path])
    ge_dt <- fread(cmd=sed_cmd,header=F,col.names = c('ensembl_gene','expr'))
    return(ge_dt)
  },mc.cores = ncpu)

  names(exprs) <- as.character(fpkm_dt$cohort_submitter_id)
  return(exprs)
}

load_fpkm_rnaseq_on_goi <- function(fpkm_dt,ncpu=1,debug=F) {
  
  if (debug){browser()}
  mount_prefix <- get_mount_dir()

  source(file.path(Sys.getenv('R_UTIL_APA'),'lib_PRs.r'))
  
  prepropa_rd <- file.path(mount_prefix,'apa_atingLab2019','01_polyAseq','01_wkd','out','03_CallApa','output','prepropa.rd')
  
  apa_ann_rd <- file.path(mount_prefix,'apa_atingLab2019','01_polyAseq','01_wkd','out','04_AnnoApa','output','apa.ann.rd')

  goi_map <- map_refgene_to_ensembl(prepropa_rd,apa_ann_rd)
  
  message('initialize read count matrix [goi x submitter]')
  goi_exprs <- mclapply(1:nrow(fpkm_dt),function(i){
    sed_cmd <- sprintf("zcat %s | sed \'s/\\.[0-9]\\+//1\'",fpkm_dt[i,file_path])
    ge_dt <- fread(cmd=sed_cmd,header=F,col.names = c('ensembl_gene','expr'))
    
    return(ge_dt[match(goi_map$ensembl_gene,ge_dt$ensembl_gene),expr])
  },mc.cores = ncpu)
  
  goi_expr_dt <- t(as.data.table(do.call(rbind,goi_exprs))) #goi_exprs is a list of vector and thus we need rbind instead of cbind
  
  colnames(goi_expr_dt) <- as.character(fpkm_dt$cohort_submitter_id)
  rownames(goi_expr_dt) <- goi_map$gene_name
  goi_expr_dt <- sort_column(goi_expr_dt)
  
  return(goi_expr_dt)
}

map_uuid_to_submitter_id <- function(uuids){
  source(file.path(Sys.getenv('R_UTIL_APA'),'tcga_id_map.R'))
  ret <- TCGAtranslateID(uuids)
  return(ret[match(uuids,ret$file_id),"submitter_id"])
}


sort_column <- function(mytab) {
  
  rowname_bkp <- rownames(mytab)
  so_snames <- sort(colnames(mytab))
  
  if (data.class(mytab) == 'data.table') {
    mytab <- mytab[,..so_snames]
  } else if (data.class(mytab) == 'data.table') {
    mytab <- mytab[,so_snames]
  }
  colnames(mytab) <- so_snames
  
  if (!is.null(rowname_bkp)){
    rownames(mytab) <- rowname_bkp
  }
  return(mytab)
}
