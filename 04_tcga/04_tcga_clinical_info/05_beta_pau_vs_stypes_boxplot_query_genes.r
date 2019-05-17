closeAllConnections()
rm(list=ls())

library(argparse)
library(data.table)
library(reshape2)
library(parallel)
library(ggplot2)
library(grid)
library(gridExtra)
library(handy)

source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))

jmessage('loading a rd file ...')
rd_fn <- './output_04/annotated_clinical_trm.rd'
load(rd_fn) #pau_dt,methyl_dt

load('./output_01/546goi_topk_masked8000_apa.rd')
apa<-apa[type == 'shared',]

methyl_probe_rd <- '../546goi_all/methyl_probes.rd'
if (file.exists(methyl_probe_rd)) {
  load(methyl_probe_rd) #tcga_methyl_probes
} else {
  load('../546goi_all/tcga_methyl_pau.rd')
  tcga_methyl_probes <- unique(goi$methyl_prof_dt[,c('chr','start','end','strand','probe_id','gene_name')])
  tcga_methyl_probes <- tcga_methyl_probes[order(chr,start,end)]
  save(tcga_methyl_probes,file=methyl_probe_rd,compress=T)
}

jmessage('reformatting pau data table ...')

# ============================
total_test_regions <- 651
test_ids <- colnames(pau_dt)[1:total_test_regions]

pau_dt$cohort <- sapply(pau_dt$project_id,function(project_id) {
  return(strsplit(project_id,'-')[[1]][2])
})

cfeatures <- c("tumor_stage","gender","race","ethnicity","primary_diagnosis","year_of_birth","year_of_death")
pau_cols <- c("cohort","sample_type",cfeatures)

pau_ti_dts <- mclapply(test_ids,function(test_id){
  ret <- pau_dt[,c(pau_cols,test_id),with=F]
  setnames(ret,test_id,'pau')
  ret$test_id <- test_id
  return(ret)
},mc.cores = 16)


pau_dt <- rbindlist(pau_ti_dts)
rm(pau_ti_dts)

tumor_stages <- list()
tumor_stages$i <- c('stage i','stage ia','stage ib','stage i')
tumor_stages$ii <- c('stage ii','stage iia','stage iib','stage iic')
tumor_stages$iii <- c('stage iii','stage iiia','stage iiib','stage iiic')
tumor_stages$iv <- c('stage iv','stage iva','stage ivb','stage ivc')

tumor_stages_names <- names(tumor_stages)
jmessage('box plotting pau vs. sample type ...')


# =====================================
pau_dt$tumor_stage_group <- as.character()
#pau_dt[sample_type %in% c('Solid Tissue Normal'),sample_type_simple:='Normal']
#pau_dt[sample_type %in% c('Primary Tumor','Recurrent Tumor'),sample_type_simple:='Cancer']
pau_dt <- pau_dt[sample_type %in% c('Primary Tumor','Recurrent Tumor'),]
unique(pau_dt$tumor_stage)

for (i in 1:length(tumor_stages_names)) {
  pau_dt[tumor_stage %in% tumor_stages[[i]],tumor_stage_group:=tumor_stages_names[[i]]]
}

pau_dt <- pau_dt[(!is.na(pau) & pau <= 1.0 & !is.na(tumor_stage_group)),]
pau_dt$gene_name <- apa[match(pau_dt$test_id,apa$test_id),gene_name]


# =====================================
methyl_dt <- methyl_dt[sample_type %in% c('Primary Tumor','Recurrent Tumor'),]
methyl_dt$tumor_stage_group <- as.character()
for (i in 1:length(tumor_stages_names)) {
  methyl_dt[tumor_stage %in% tumor_stages[[i]],tumor_stage_group:=tumor_stages_names[[i]]]
}
methyl_dt <- methyl_dt[(!is.na(nbeta) & !is.na(tumor_stage_group)),]

cohort <- mclapply(methyl_dt$project_id,function(project_id) {
  return(strsplit(project_id,'-')[[1]][2])
},mc.cores = 16)
methyl_dt$cohort <- unlist(cohort)
rm(cohort)

methyl_dt <- methyl_dt[,c('probe_id','gene_name','nbeta','cohort','sample_type','tumor_stage_group',"tumor_stage","gender","race","ethnicity","primary_diagnosis","year_of_birth","year_of_death")]

methyl_dt$gene_name <- probe_gene_map[match(methyl_dt$probe_id,probe_gene_map$probe_id),gene_name]

methyl_dt_query <- methyl_dt[gene_name == 'HEATR2',]

midx <- match(methyl_dt_query$probe_id,tcga_methyl_probes$probe_id)

methyl_dt_query$chr <- tcga_methyl_probes[midx,chr]
methyl_dt_query$start <- tcga_methyl_probes[midx,start]
methyl_dt_query$end <- tcga_methyl_probes[midx,end]

query.chr <- 'chr7'
query.start <- 807595
query.end <- 810959

methyl_dt_query_sub <- methyl_dt_query[chr==query.chr & start>=query.start & end<=query.end,]

methyl_dt_query_probe <- split(methyl_dt_query_sub,by=c("probe_id"))

dodge <- position_dodge(width = 0.4)

P <- length(methyl_dt_query_probe)
figh <- list()
for (i in 1:P) {
  
  probe_id <- methyl_dt_query_probe[[i]]$probe_id[1]
  probe_pos <- sprintf('%s:%d',methyl_dt_query_probe[[i]]$chr,methyl_dt_query_probe[[i]]$start)
  
  p <- ggplot(data = methyl_dt_query_probe[[i]],
              aes(x = cohort, 
                  y = nbeta,
                  color = as.factor(tumor_stage_group))) +
    geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + 
    #geom_jitter(shape=16, position=dodge) +
    ggtitle(sprintf("%s/%s",probe_id,probe_pos)) +
    theme(legend.position="bottom") +
    handy::ggnice()
  
  figh[[probe_id]] <- p
  
}

mtitle <- sprintf('tumor stage vs. methylation rate in HEATR2 (%s:%d-%d)',query.chr,query.start,query.end)

grid.arrange(grobs=figh,
             nrow=P,
             top = textGrob(mtitle))
