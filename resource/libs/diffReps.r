# https://github.com/shenlab-sinai/diffreps

library(argparse)
library(data.table)
library(parallel)
library(Rsubread)
source(file.path(Sys.getenv('R_UTIL_APA'),'lib_apps.R'))

# diffReps - Detect Differential Sites from ChIP-seq with Biological Replicates.
# 
# Usage: diffReps.pl --treatment bed_file1...[bed_fileN] --control bed_file1...[bed_fileN] --report output_results \
# [--gname genome_name|--chrlen chrom_length_file]
# 
# **Chromosome lengths can be specified through a text file or given a genome name.
# **Currently built-in genomes: mm9, hg19, rn4.
# 
# Optional parameters(defaults in parentheses):
#   - Background samples(DNA input or IgG control):
#   --btr            Treatment group background: bed_file1...[bed_fileN].
# --bco            Control group background: bed_file1...[bed_fileN].
# Hint: If background is only specified for one group, it will automatically be used for both groups.
# 
# - Genomic region parameters:
#   --mode(peak)     Scanning mode: a selection implies a different window size.
# Set window and step size manually to override.
# (p)eak      (=1000)  Histone mark peak (Default).
# (n)ucleosome(=200)   Single nucleosome (+DNAlinker).
# (b)lock     (=10000) Large chromatin modification block.
# --window(1000)   Window size (default=Histone mark peak size).
# --step(1/10 win) Window moving step size.
# --gap(0)         Gap allowed between two consecutive windows.
# 
# - Background filtering using: mean + nsd*deviation.
# --std            Use standard estimation of mean and deviation (Default=Robust estimation).
# In robust estimation, median absolute deviation is used in place of standard deviation.
# --nsd(broad)     Z-score cutoff for low read count. Choose from two default modes or set your own.
# (b)road     (=2)   Broad peak such as H3K36me3.
# (s)harp     (=20)  Sharp peak such as H3K4me3 or Transcription factors.
# --alpha(0.05)    Alpha for right-trimmed mean, must be in: [0, 0.5).
# --bkg(0)         Use fold enrichment vs. background as filter instead. Set a float number such as 2.0 here.
# Default is to use the Z-score as filter.
# 
# - Statistical testing parameters:
#   --meth(nb)       Statistical test (nb=Negative binomial; gt=G-test; tt=T-test; cs=Chi-square test).
# --pval(0.0001)   P-value cutoff for significant windows.
# 
# - Normalization can be done externally and be supplied as a text file:
#   --nrpass(1)      Do normalization on bins pass nsd cutoff?
#   --norm           File name to specify pre-determined norm constants (Default=Estimate by diffReps).
# 
# - Misc. parameters:
#   --frag(100)      ChIP-seq library fragment size. Use to shift read positions.
# --nproc(1)       Number of processors to use.
# --noanno         Switch off genomic annotation for differential sites (Default=Do annotation).
# --nohs           Switch off looking for chromatin modification hotspots (Default=Find hotspots).


debug <- FALSE
if (!debug) {
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  source_root <- dirname(dirname(scriptPath))
  parser <- ArgumentParser(description='diffChipSeq_lite')
  
  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                      dest="ncpu", default = 1,
                      help="number of cpus to utilize [1]")

  parser$add_argument("-g", "--genome", type="character", default='hg19', required=FALSE,
                      dest="genome",
                      help="genome[hg19]")
  
  parser$add_argument("-c", "--bam_co", type="character", required=TRUE,
                      dest="bam_co",
                      help="bam control")
  
  parser$add_argument("-t", "--bam_tr", type="character", required=TRUE,
                      dest="bam_tr",
                      help="bam treatment")
  
  parser$add_argument("-m", "--test_model", type="character", default="gt", required=FALSE, 
                      dest="test_model",
                      help="[gt], nb, tt, cs")
  
  parser$add_argument("-w", "--window", type="integer", required=FALSE,
                      dest="window", default = 200,
                      help="window length")
  
  parser$add_argument("-f", "--frag", type="integer", required=FALSE,
                      dest="frag", default = 150,
                      help="fragment length [150]")
  
  parser$add_argument("-o", "--out_dir", type="character", required=TRUE,
                      dest="out_dir",
                      help="output directory")
  
  args <- parser$parse_args()
} else { #debug
  
  # ----------------
  args <- data.table(s1="",
                     si1="",
                     s2="",
                     si2="",
                     name="",
                     group="",
                     bind_mode="",
                     encode_analy_dir="")
  
  args$encode_analy_dir <- "/home/hongc2/projects/apa/activemotif_fastq/encode_analy_out"
  
  args$group <- "ctrl,exp"
  args$name <- "hct116_ctcf,dko_ctcf"
  args$bind_mode <- "TF"
  
  args$s1 <- "01_0307_00CWCCF_HCT116_AHT_CTCF_hs_i72"
  args$si1 <- "37_030J_00CWCCF_HCT116_AHT_Input_hs_i84"
  
  args$s2 <- "02_0308_00CWCCF_DKO_AHT_CTCF_hs_i73"
  args$si2 <- "38_030K_00CWCCF_DKO_AHT_Input_hs_i85"
}

bed_co <- sprintf('%s.bed',args$bam_co)
bedtools_bamtobed_qc(args$bam_co,bed_co)

bed_tr <- sprintf('%s.bed',args$bam_tr)
bedtools_bamtobed_qc(args$bam_tr,bed_tr)

diffReps(bed_co,bed_tr,
         outfile=args$outfile,
         args$genome,
         test_model=args$test_model,
         window=args$window,
         frag=args$frag)
