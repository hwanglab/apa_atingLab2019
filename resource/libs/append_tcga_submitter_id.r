library(argparse)
library(data.table)
source(file.path(Sys.getenv('R_UTIL_APA'),'tcga_id_map.R'))

if (TRUE) { #debug
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  
  parser <- ArgumentParser(description='prepropa')
  
  parser$add_argument("-i", "--manifest_fn", type="character", required=TRUE,
                      dest="manifest_fn",
                      help="manifest file")
  
  parser$add_argument("-d", "--data_root_dir", type="character", required=TRUE,
                      dest="data_root_dir",
                      help="TCGA actual data root dir")
  
  parser$add_argument("-o", "--rd_fn", type="character", required=TRUE,
                      dest="rd_fn",
                      help="rd file to save the output")
  
  args <- parser$parse_args()
  
} else {
  args = data.table(manifest_fn='/home/hongc2/projects/data/TCGA/blca/methylArray/manifest/TCGA-BLCA-MethArray-gdc_manifest.2018-06-14.txt',
                    data_root_dir='/mnt/lustre/tinga/TCGA/TCGA-MethArray/BLCA',
                   rd_fn='/home/hongc2/projects/data/TCGA/blca/methylArray/manifest/TCGA-BLCA-MethArray-gdc_manifest.2018-06-14.rd')
}

message("reading manifest file ...")
manifest <- fread(args$manifest_fn,header=T)
res <- TCGAtranslateID(manifest$id)

manifest$submitter_id <- res[match(manifest$id,res$file_id),'submitter_id']

manifest$fpath <- mapply(function(uuid,fbase,tcga_data_dir){
  fpath <- file.path(tcga_data_dir,uuid,fbase)
  },
  manifest$id,
  manifest$filename,
  MoreArgs = list(args$data_root_dir)
)

save(manifest,file=args$rd_fn)
message(sprintf("rd file [%s] saved",args$rd_fn))

