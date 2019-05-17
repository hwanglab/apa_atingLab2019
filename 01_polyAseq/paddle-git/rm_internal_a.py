#!/usr/bin/env python

import sys,os
# Needs this due to version weirdness on the HPC cluster, other users should remove
#sys.path.insert(0,"/home/bhasinj/local/lib/python2.7/site-packages")
import re
# import gzip
import argparse
# from numpy import array,in1d,savetxt
from collections import defaultdict
import csv
import sys
# import itertools
import pysam
# import pandas as pd
# import pybedtools
from pyfaidx import Fasta
from Bio.Seq import Seq

# --------------------------------------------------------------------
# Hexamer finding function
def checkhex(seq,dorev):
	# Will look on the revcomp of seq if dorev==True
	seq = Seq(seq)
	if dorev:
		seq = seq.reverse_complement()
	#print "Checking seq: " + seq
	
	hexes = ["AATAAA","ATTAAA","AGTAAA","TATAAA","CATAAA","GATAAA","AATATA","AATACA","AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA"]

	out = "NONE"
	for checkhex in hexes:
		if seq.find(checkhex)>-1:
			out = checkhex
			break
	#print "Found hex: " + out
	return out
# --------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Remove reads due to internal priming (genomic A).", epilog="Creates a new BAM file without internally primed reads.")
parser.add_argument('inpath', metavar='PATH_TO_BAM', help='Path to .bam file to use as input.')

parser.add_argument('--ref_genome_fa', dest='ref_genome_fa', required=True, help='Path to reference genome sequence file')

parser.add_argument('--outdir', dest='outdir', metavar='PATH_TO_OUTPUT_DIR', action='store', default=".", help='Path to place out BAM and summary CSV (default: current working directory)')
parser.add_argument('--log', dest='dolog', action='store_true', help='Redirect all output to a log file')
parser.add_argument('--print', dest='doprint', action='store_true', help='Print debug log of sequences used')
parser.add_argument('--abam', dest='doabam', action='store_true', help='Also save a second BAM file (hasGenA) in addition to the noGenA BAM.')
parser.add_argument('--win', dest='win', metavar='WINDOW_SIZE', action='store', default=8, help='Size of downstream window (with respect to genomic strand) in which to search for genomic As. (default: 8)')

parser.add_argument('--chr_header', dest='chr_header', metavar='CHR_HEADER', action='store', default=1, help='does bam file have chr in contig header?[1:Yes],0:No')

parser.add_argument('--mm', dest='mm', metavar='ALLOWED_MISMATCHES', action='store', default=1, help='Number of mismatches to allow when searching for genomic A in the window (default: 1)')
parser.add_argument('--flank', dest='flank', metavar='FLANK_SIZE', action='store', default=50, help='Size of flank regions in bp (default: 50 bp)')
parser.add_argument('--orfirst', dest='orfirst', metavar='OR_FIRST_BP', action='store', default=-1, help=' (default: disabled)')
parser.set_defaults(dolog=False)

args = parser.parse_args()

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

outbase = re.sub(r'\.bam$', '', args.inpath)
outbase = os.path.basename(outbase)
outbase = os.path.join(args.outdir,outbase)

if args.dolog==True:
	log_file = open(outbase + "_rm_internal_a.log","w")
	sys.stdout = log_file

print "Reading BAM file: " + args.inpath
print "Basename for output files: " + outbase
sys.stdout.flush()

if args.doprint:
	print_a = open(outbase + "_print_genomicA.txt",'w')
	print_noa = open(outbase + "_print_noGenomicA.txt",'w')

#fa = Fasta('/home/bhasinj/lustre/Common/hg19.fa',as_raw=True,sequence_always_upper=True)
fa = Fasta(args.ref_genome_fa,as_raw=True,sequence_always_upper=True)

# --------------------------------------------------------------------
# Step 1: Iterate over the BAM file and extract key values. Save them into PANDAS DF.
samfile = pysam.AlignmentFile(args.inpath, "rb")
newbam_noa = pysam.AlignmentFile(outbase + "_noGenomicA.bam", "wb", template=samfile)
newbam_a = pysam.AlignmentFile(outbase + "_genomicA.bam", "wb", template=samfile)

bamrow = []
chrs = []
starts = []
ends = []
isrevs = []
#seqs = []

print "Loading data from BAM to RAM"

want = []

if args.chr_header == 1:
	for i in xrange(1,23):
		want.append("chr" + str(i))
	want.append("chrX")
	want.append("chrY")
else:
	for i in xrange(1,23):
		want.append(str(i))
	want.append("X")
	want.append("Y")

want = set(want)

a_hex_hist = defaultdict(int)
noa_hex_hist = defaultdict(int)

recno=0
for read in samfile.fetch():
	my_chr = samfile.getrname(read.rname)
	my_pos = read.pos
	my_isrev = read.is_reverse
	my_len = len(read.seq)

	if my_chr in want:
		# Both coords 1-based
		my_start = my_pos+1
		my_end = my_pos+my_len
	
		# Pull seq
		read_seq = fa[my_chr][my_start-1:my_end]

		if my_isrev==False:
			win_start = my_end+1
			win_end = my_end+int(args.win)
		elif my_isrev==True:
			win_start = my_start-int(args.win)
			win_end = my_start-1

		win_seq = fa[my_chr][win_start-1:win_end]

		firsts = False

		if my_isrev==False:
			my_strand = "+"
			f1_start = my_start-int(args.flank)
			f1_end = my_start-1
			f2_start = win_end+1
			f2_end = win_end+int(args.flank)
			f1_seq = fa[my_chr][f1_start-1:f1_end]
			f2_seq = fa[my_chr][f2_start-1:f2_end]
			win_a = win_seq.count("A")

			myhex = checkhex(seq=f1_seq,dorev=False)

			if(int(args.orfirst)>0):
				first_seq = win_seq[0:int(args.orfirst)]
				if first_seq.count("A")==int(args.orfirst):
					firsts = True

			if (win_a >= (int(args.win)-int(args.mm)))|(firsts==True):
				newbam_a.write(read)
				a_hex_hist[myhex] += 1
				if args.doprint:
					print_a.write("[%s] f1:%s-%s read:%s-%s win:%s-%s f2:%s-%s hexamer:%s\n" % (my_chr,f1_start,f1_end,my_start,my_end,win_start,win_end,f2_start,f2_end,myhex))
					print_a.write("[%s:%s-%s (%s)] flank:%s read:%s win(%s A):%s flank:%s\n\n" % (my_chr,f1_start,f2_end,my_strand,f1_seq,read_seq,win_a,win_seq,f2_seq))
			else:
				newbam_noa.write(read)
				noa_hex_hist[myhex] += 1
				if args.doprint:
					print_noa.write("[%s] f1:%s-%s read:%s-%s win:%s-%s f2:%s-%s hexamer:%s\n" % (my_chr,f1_start,f1_end,my_start,my_end,win_start,win_end,f2_start,f2_end,myhex))
					print_noa.write("[%s:%s-%s (%s)] flank:%s read:%s win(%s A):%s flank:%s\n\n" % (my_chr,f1_start,f2_end,my_strand,f1_seq,read_seq,win_a,win_seq,f2_seq))


		elif my_isrev==True:
			my_strand = "-"
			f1_start = win_start-int(args.flank)
			f1_end = win_start-1
			f2_start = my_end+1
			f2_end = my_end+int(args.flank)
			f1_seq = fa[my_chr][f1_start-1:f1_end]
			f2_seq = fa[my_chr][f2_start-1:f2_end]
			win_a = win_seq.count("T")

			myhex = checkhex(seq=f2_seq,dorev=True)

			if(int(args.orfirst)>0):
				first_seq = win_seq[0:int(args.orfirst)]
				if first_seq.count("A")==int(args.orfirst):
					firsts = True

			if (win_a >= (int(args.win)-int(args.mm)))|(firsts==True):
				newbam_a.write(read)
				a_hex_hist[myhex] += 1
				if args.doprint:
					print_a.write("[%s] f1:%s-%s win:%s-%s read:%s-%s f2:%s-%s hexamer:%s\n" % (my_chr,f1_start,f1_end,win_start,win_end,my_start,my_end,f2_start,f2_end,myhex))
					print_a.write("[%s:%s-%s (%s)] flank:%s win(%s A):%s read:%s flank:%s\n\n" % (my_chr,f1_start,f2_end,my_strand,f1_seq,win_a,win_seq,read_seq,f2_seq))
			else:
				newbam_noa.write(read)
				noa_hex_hist[myhex] += 1
				if args.doprint:
					print_noa.write("[%s] f1:%s-%s win:%s-%s read:%s-%s f2:%s-%s hexamer:%s\n" % (my_chr,f1_start,f1_end,win_start,win_end,my_start,my_end,f2_start,f2_end,myhex))
					print_noa.write("[%s:%s-%s (%s)] flank:%s win(%s A):%s read:%s flank:%s\n\n" % (my_chr,f1_start,f2_end,my_strand,f1_seq,win_a,win_seq,read_seq,f2_seq))

		recno += 1
		if (recno % 100000 == 0):
			print "Done with Record Number: " + str(recno)
			sys.stdout.flush()
			#break

# Clean up and save CSV
samfile.close()
newbam_a.close()
newbam_noa.close()

#print a_hex_hist
#print noa_hex_hist

print "Saving hexamer histogram CSV"
csvpath = outbase + "_HexamerCounts.csv"
print "Writing histogram: " + csvpath
sys.stdout.flush()
writer = csv.writer(open(csvpath, 'w'))
writer.writerow(["group","hexamer","count"])
for key, value in a_hex_hist.items():
	writer.writerow(["genomicA", key, value])
for key, value in noa_hex_hist.items():
	writer.writerow(["noGenomicA", key, value])

if args.doprint:
	print_a.close()
	print_noa.close()

if args.dolog==True:
	log_file.close()

