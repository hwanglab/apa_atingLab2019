#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018
#inputs:
#- activemotif factor_path active region excel file
#- maurano table w/ annotation by hong pipeline comparing with dko samples
#method:
#- activemotif excel file contains HCT and DKO ctcf binding site flag (0:not present, 1:present) at column AL and AM, respectively. first, compares the column AL with maurano table by checking if each region is overlapped and count how concordant the peak calling is each other

library(data.table)
library(stringr)

load_spr_region_tsv <- function(spr_region_fn,to_ucsc_contig=FALSE) {
	edger = fread(spr_region_fn,header = TRUE)
	names(edger)[1] = "region"

	#split string vector with more than one delimiter
	region = transpose(setDT(strsplit(edger$region, "[_,-]+")))

	#assign anew column names
	names(region)[1:3] = c("chrom","chromStart","chromEnd")
	
	if (to_ucsc_contig) {
		edger$chrom = paste0('chr',region$chrom)
	} else {
		edger$chrom = region$chrom
	}
	
	edger$chromStart = as.integer(region$chromStart)
	edger$chromEnd = as.integer(region$chromEnd)
	edger$region = NULL
	edger
}