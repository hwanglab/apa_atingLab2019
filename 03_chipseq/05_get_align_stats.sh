#!/bin/bash -l

Rscript ./get_encode_align_stats.r -i fastq/kundaje_encode -o fastq/kundaje_encode/alignment_stats
Rscript ./get_mbdseq_align_stats.r
