require(Biostrings)
require(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(stringi)

fetch_fastas <- function(regions) {
	
	genome <- BSgenome.Hsapiens.UCSC.hg19
	
	#jmessage('Preparing GRanges object ...')
	gr <- GRanges(seqnames = regions$chrom,
								ranges = IRanges(start = regions$chromStart, 
																 end=regions$chromEnd,
																 names = regions$chrom)
								)
	#jmessage(sprintf('Fetching [%d] genome sequences from hg19',dim(regions)[1]))
	seqs <- Views(genome, gr)
	
	#jmessage('Converting to data.table format to return ...')
	seqs
}

write_fasta_file <- function(dt,fasta_out_gz){
	dt$header = NA
	dt$header = paste(dt$seqnames,dt$start,dt$end,sep = '_')
	fw = gzfile(fasta_out_gz,"w")
	jmessage('writing FASTA file [%s] ...',fasta_out_gz)
	for (i in 1:nrow(dt)){
		cat(sprintf(">%s\n%s\n",dt[i]$header,dt[i]$dna), file = fw)
	}
	close(fw)
	jmessage('Done')
}

best_match_freq <- function(BsgViewSet,target_freq) {
	dimer_freqs = dinucleotideFrequency(BsgViewSet)
	idx0 = which(rowSums(dimer_freqs)>0)
	j = which.min(rowSums(abs(dimer_freqs[idx0,] - target_freq)))
	return(idx0[j])
}

slide_window_seq <- function(aheader,astring,W,slide,fw) {
	S = str_length(astring)
	s1 = 1 - slide
	j = 0
	run_flag = TRUE
	while (run_flag) {
		s1 = s1 + slide
		s2 = s1 + W
		if (s2 > S) {
			astringj = str_sub(astring, start = s1, end = S)
			cat(sprintf(">%s_%d\n%s%s\n",aheader,j,astringj,stri_dup('N',s2-S)), file = fw)
			run_flag = FALSE
		} else {
			astringj = str_sub(astring, start = s1, end = s2)
			cat(sprintf(">%s_%d\n%s\n",aheader,j,astringj), file = fw)
		}
		j = j + 1
	}
}

slide_window_fasta <- function(in_fa_gz,W,slide,out_fa_gz) {
	
	jmessage(sprintf('slide %d window by %d on FASTA file [%s] ...',W,slide,in_fa_gz))
	inFa <- readDNAStringSet(in_fa_gz)
	fw = gzfile(out_fa_gz,"w")
	dt <- data.table(head=names(inFa), seq=paste(inFa))
	for (r in 1:nrow(dt)) {
		slide_window_seq(dt$head[r],dt$seq[r],W,slide,fw)
	}
	close(fw)
}