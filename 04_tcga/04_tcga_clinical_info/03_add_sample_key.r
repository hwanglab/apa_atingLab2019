setwd("~/projects/apa/src/s14_tcga_clinical_info")
closeAllConnections()
rm(list=ls())

library(argparse)
library(data.table)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))

# ----------------------
jmessage('loading clinical information ...')
rd_fn <- './output_01/tcga_rm_clin_info.rd'
if (file.exists(rd_fn)) {
  load(rd_fn) #trm_samples (TCGA samples obtained from mRNA and 4500 array), tcga_clin (TCGA clinical information table)
} else {
  stop('run 01_retrieve_clinical_info.r first!')
}

head(trm_samples$sample)
head(trm_samples$submitter_id)

# ----------------------
jmessage('locating masked_gtop8000 rd files ...')
tcga_masked_rd_fn <- '../546goi_all/masked_gtop8000.rd'
load(tcga_masked_rd_fn)
test_id <- goi$pau$ratio$test_id
goi$pau$ratio$test_id <- NULL
sample_id <- colnames(goi$pau$ratio)
pau_dt <- as.data.table(t(goi$pau$ratio))
colnames(pau_dt) <- as.character(test_id)
rownames(pau_dt) <- sample_id

methyl_dt <- goi$methyl_prof_dt

if (!dir.exists('./output_03')) {
  dir.create('./output_03')
}

rd_fn <- './output_03/add_sample_key.rd'

pau_dt$cohort_sample_id <- rownames(pau_dt)

methyl_dt$submitter_id <- trm_samples$sample[match(methyl_dt$cohort_submitter_id,trm_samples$sample$cohort_sample_id),submitter_id]

methyl_dt$sample_id <- trm_samples$sample[match(methyl_dt$cohort_submitter_id,trm_samples$sample$cohort_sample_id),sample_id]

pau_dt$submitter_id <- trm_samples$sample[match(pau_dt$cohort_sample_id,trm_samples$sample$cohort_sample_id),submitter_id]

pau_dt$sample_id <- trm_samples$sample[match(pau_dt$cohort_sample_id,trm_samples$sample$cohort_sample_id),sample_id]

pau_dt$cohort <- trm_samples$sample[match(pau_dt$cohort_sample_id,trm_samples$sample$cohort_sample_id),cohort]

methyl_dt$project_id <- paste0('TCGA-',methyl_dt$cohort)
pau_dt$project_id <- paste0('TCGA-',pau_dt$cohort)
methyl_dt$cohort <- NULL
pau_dt$cohort <- NULL

save(pau_dt,methyl_dt,file=rd_fn,compress=T)

