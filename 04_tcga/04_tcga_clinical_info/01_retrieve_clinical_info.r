library(data.table)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))

decompose_sample_ids <- function(cohort_submitter_ids) {
  ret <- as.data.table(strsplit(cohort_submitter_ids,split = '_'))
  ret <- as.data.table(t(ret))
  ret2 <- as.data.table(strsplit(ret$V2,split='-'))
  ret2 <- as.data.table(t(ret2))
  submitter_id <- paste0(ret2$V1,'-',ret2$V2,'-',ret2$V3)
  trm_sample_dt <- data.table(submitter_id=submitter_id,
                              cohort=ret$V1,
                              cohort_sample_id=cohort_submitter_ids,
                              sample_id=ret$V2,
                              sample_type=ret2$V4)
  rm(ret)
  rm(ret2)
  
  trm_sample <- list()
  trm_sample$sample <- trm_sample_dt
  trm_sample$submitter_id <- unique(trm_sample_dt$submitter_id)
  rm(trm_sample_dt)
  
  return(trm_sample)
}

load_tcga_clinical_info <- function(tcga_clinical_info_dir,cohorts) {
  
  patient_dts <- list()
  sample_dts <- list()
  
  for (cohort in cohorts) {
    outd <- file.path(tcga_clinical_info_dir,cohort)
    if (!file.exists(outd)) {
      dir.create(outd)
    }
    
    prefix2 <- c('clinical','biospecimen')
    prefix3 <- c('clinical','sample')
    
    tsv_files <- lapply(1:2,function(i) {
      prefix <- prefix2[i]
      file_prefix <- prefix3[i]
      my_fn <- sprintf('%s.%s.tar.gz',prefix,cohort)
      my_fpath <- file.path(tcga_clinical_info_dir,my_fn)
      
      clinic_tsv <- file.path(outd,sprintf('%s.tsv',file_prefix))
      
      if (!file.exists(clinic_tsv)){
        cmd <- sprintf('tar -zxvf %s -C %s/',my_fpath,outd)
        message(cmd)
        system(cmd)
      }
      return(clinic_tsv)
    })
    
    clinic_tsv <- file.path(outd,sprintf('%s.tsv',prefix3[1]))
    patient_dts[[cohort]] <- fread(clinic_tsv)
    
    sample_tsv <- file.path(outd,sprintf('%s.tsv',prefix3[2]))
    sample_dts[[cohort]] <- fread(sample_tsv)[,list(sample_submitter_id,case_submitter_id,project_id,sample_type_id,sample_type)]
    
  }
  
  clin <- list()
  clin$patient <- unique(rbindlist(patient_dts))
  clin$patient[tumor_stage=="--",tumor_stage:="not reported"]
  clin$sample <- unique(rbindlist(sample_dts))
  return(clin)
}

jmessage('main() --------------------')

if (!dir.exists('./output_01')) {
  dir.create('./output_01')
}

# -------------------------------
jmessage('locating masked_gtop8000 rd files ...')
tcga_masked_rd_fn <- '../546goi_all/masked_gtop8000.rd'
load(tcga_masked_rd_fn)

apa_rd_fn <- './output_01/546goi_topk_masked8000_apa.rd'
apa <- goi$apa
save(apa,file=apa_rd_fn,compress = T)

cohort_submitter_ids <- unique(goi$methyl_prof_dt$cohort_submitter_id)

trm_samples <- decompose_sample_ids(cohort_submitter_ids)
rm(cohort_submitter_ids)

# ----------------
jmessage('locating TCGA clinical info direcotry (tsv) ...')

tcga_clinical_info_dir <- '../../resource/tcga_clinical_data'
cohorts <- c("BLCA","BRCA","COAD","ESCA","KIRC","LUAD","LUSC","PRAD","SKCM","STAD","UCEC")

tcga_clin <- load_tcga_clinical_info(tcga_clinical_info_dir,cohorts)

rd2_fn <- './output_01/tcga_rm_clin_info.rd'
save(tcga_clin,trm_samples,file=rd2_fn,compress = T)
