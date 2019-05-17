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
  stop('run 01_retrieve_clinical_info.r !')
}

head(trm_samples$sample)
head(trm_samples$submitter_id)

rd_fn <- './output_03/add_sample_key.rd'
if (file.exists(rd_fn)) {
  load(rd_fn) #pau_dt, methyl_dt
} else {
  stop('run 03_add_sample_key.r !')
}

# ------------
if (!dir.exists('./output_04')) {
  dir.create('./output_04')
}

rd_fn <- './output_04/annotated_clinical_trm.rd'

jmessage('quick looking at the table format ...')
head(tcga_clin$patient)
head(tcga_clin$sample)
head(pau_dt[,list(cohort_sample_id,submitter_id,sample_id,project_id)])
head(methyl_dt)

jmessage('generaring a key to match with ...')
pau_dt$key_id <- paste0(pau_dt$project_id,':',pau_dt$submitter_id)
methyl_dt$key_id <- paste0(methyl_dt$project_id,':',methyl_dt$submitter_id)

tcga_clin$patient$key_id <- paste0(tcga_clin$patient$project_id,':',tcga_clin$patient$submitter_id)
tcga_clin$sample$key_id <- paste0(tcga_clin$sample$project_id,':',tcga_clin$sample$case_submitter_id)

jmessage('annotate patient info (gender,tumor_stage,) ...')
pau_dt$gender <- tcga_clin$patient[match(pau_dt$key_id,key_id),gender]
pau_dt$tumor_stage <- tcga_clin$patient[match(pau_dt$key_id,key_id),tumor_stage]
pau_dt$race <- tcga_clin$patient[match(pau_dt$key_id,key_id),race]
pau_dt$ethnicity <- tcga_clin$patient[match(pau_dt$key_id,key_id),ethnicity]
pau_dt$primary_diagnosis <- tcga_clin$patient[match(pau_dt$key_id,key_id),primary_diagnosis]
pau_dt$year_of_birth <- tcga_clin$patient[match(pau_dt$key_id,key_id),year_of_birth]
pau_dt$year_of_death <- tcga_clin$patient[match(pau_dt$key_id,key_id),year_of_death]

methyl_dt$gender <- tcga_clin$patient[match(methyl_dt$key_id,key_id),gender]
methyl_dt$tumor_stage <- tcga_clin$patient[match(methyl_dt$key_id,key_id),tumor_stage]
methyl_dt$race <- tcga_clin$patient[match(methyl_dt$key_id,key_id),race]
methyl_dt$ethnicity <- tcga_clin$patient[match(methyl_dt$key_id,key_id),ethnicity]
methyl_dt$primary_diagnosis <- tcga_clin$patient[match(methyl_dt$key_id,key_id),primary_diagnosis]
methyl_dt$year_of_birth <- tcga_clin$patient[match(methyl_dt$key_id,key_id),year_of_birth]
methyl_dt$year_of_death <- tcga_clin$patient[match(methyl_dt$key_id,key_id),year_of_death]

jmessage('annotate sample info (sample_type) ...')
pau_dt$sample_type_id <- tcga_clin$sample[match(pau_dt$key_id,key_id),sample_type_id]
pau_dt$sample_type <- tcga_clin$sample[match(pau_dt$key_id,key_id),sample_type]

methyl_dt$sample_type_id <- tcga_clin$sample[match(methyl_dt$key_id,key_id),sample_type_id]
methyl_dt$sample_type <- tcga_clin$sample[match(methyl_dt$key_id,key_id),sample_type]

save(pau_dt,methyl_dt,file=rd_fn,compress = T)
