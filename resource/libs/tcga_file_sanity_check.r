library(argparse)

parser <- ArgumentParser(description='perform MD5 matching for TCGA files downloaded')

parser$add_argument("-b", "--base_dir", type="character", required=TRUE,
                    dest="base_dir",
                    help="base directory")

parser$add_argument("-m", "--manifest_file", type="character", required=TRUE,
                    dest="manifest_file",
                    help="manifest file name")

# ----------
library(data.table)
library(tools)

args <- parser$parse_args()

dt <- fread(args$manifest_file,header=T)

fpaths <- file.path(args$base_dir,dt$id,dt$filename)

dt$md5_match_result <- md5sum(fpaths)

fwrite(dt,file=sprintf('%s.sanity_check_result.txt',args$manifiest_file),quote=F)
