#!/bin/bash -l

java -jar picard.jar CollectInsertSizeMetrics I=star/HI.1159.005.Index_13.HCT116_Aligned.sortedByCoord.out.bam O=star/HI.1159.005.Index_13.HCT116_Aligned.sortedByCoord_insert_size_metrix.txt H=star/HI.1159.005.Index_13.HCT116_Aligned.sortedByCoord_insert_size_hist.pdf

java -jar picard.jar CollectInsertSizeMetrics I=star/HI.1159.005.Index_16.DKO_Aligned.sortedByCoord.out.bam O=star/HI.1159.005.Index_16.DKO_Aligned.sortedByCoord_insert_size_metrix.txt H=star/HI.1159.005.Index_16.DKO_Aligned.sortedByCoord_insert_size_hist.pdf
