#!/usr/bin/env python

import sys,os
# Needs this due to version weirdness on the HPC cluster, other users should remove
# sys.path.insert(0,"/home/bhasinj/local/lib/python2.7/site-packages")
import re
import gzip
import argparse
from numpy import array,in1d,savetxt
from collections import defaultdict
import csv
import sys
import itertools
import pysam
import pandas as pd

parser = argparse.ArgumentParser(description="Collapse PCR duplicates using 6N sequences.", epilog="Creates a new BAM file with all unique reads and all duplicates collapsed to a single read.")
parser.add_argument('inpath', metavar='PATH_TO_BAM', help='Path to .bam file to use as input.')
parser.add_argument('--outdir', dest='outdir', metavar='PATH_TO_OUTPUT_DIR', action='store', default=".", help='Path to place out BAM and summary CSV (default: current working directory)')
parser.add_argument('--log', dest='dolog', action='store_true', help='Redirect all output to a log file')
#parser.add_argument('--mm', dest='mm', metavar='ALLOWED_MISMATCHES', action='store', default=0, help='Number of mismatches to allow when finding the A-stretch (default: 0)')
parser.set_defaults(dolog=False)

args = parser.parse_args()

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

outbase = re.sub(r'\.bam$', '', args.inpath)
outbase = os.path.basename(outbase)
outbase = os.path.join(args.outdir,outbase)

if args.dolog==True:
	log_file = open(outbase + "_collapse_dups.log","w")
	sys.stdout = log_file

print "Reading BAM file: " + args.inpath
print "Basename for output files: " + outbase
sys.stdout.flush()

# --------------------------------------------------------------------
# Step 1: Iterate over the BAM file and extract key values. Save them into PANDAS DF.
samfile = pysam.AlignmentFile(args.inpath, "rb")

bamrow = []
tags = []
chrs = []
poss = []
strands = []
lens = []
#seqs = []

print "Loading data from BAM to RAM"

recno=0
for read in samfile.fetch():
	my_tag = read.qname.split("_")[-1]
	my_chr = samfile.getrname(read.rname)
	my_pos = read.pos
	my_strand = read.is_reverse
	my_len = len(read.seq)
	
	tags.append(my_tag)
	chrs.append(my_chr)
	poss.append(my_pos)
	strands.append(my_strand)
	#seqs.append(my_seq)
	lens.append(my_len)
	bamrow.append(recno)

	recno += 1
	if (recno % 1000000) == 0:
		print "Done with Record Number: " + str(recno)
		sys.stdout.flush()
		#break
#print tags
#print chrs
#print poss
#print strands
#print seqs

samfile.close()

print "Creading pandas DF"
df = pd.DataFrame(data={'bamrow': bamrow, 'tag': tags, 'chr': chrs, 'pos': poss, 'strand': strands, 'len': lens})
#df.sort(['chr','pos','strand','tag'])

#print "Sorting pandas DF"
#df = df.sort(["chr","pos","strand","tag"])

#print df.ix[:,"bamrow"]
#print df
#df.to_csv("tmp.csv")
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Step 2: Loop over the DF and flag the duplicates

print "Detect/Collapse Duplicates"

byCPLS = df.groupby(["chr","pos","len","strand"])
byCPLST = df.groupby(["chr","pos","len","strand","tag"])
byCPST = df.groupby(["chr","pos","strand","tag"])
byCPS = df.groupby(["chr","pos","strand"])

#bypos_sizes = array(bypos.agg(len).bamrow,dtype=int)
#bytag_sizes = array(bytag.agg(len).bamrow,dtype=int)
byCPLS_sizes = byCPLS.agg(len)
byCPLST_sizes = byCPLST.agg(len)
byCPST_sizes = byCPST.agg(len)
byCPS_sizes = byCPS.agg(len)

bytag_want = set(byCPLST.first().bamrow)

print "Saving Duplication Vectors"
#savetxt(outbase + "_ByPosSizes.csv", bypos_sizes, delimiter=",", fmt='%i')
#savetxt(outbase + "_ByTagSizes.csv", bytag_sizes, delimiter=",", fmt='%i')

byCPLS_sizes.to_csv(outbase + "_ByCPLS.csv")
byCPLST_sizes.to_csv(outbase + "_ByCPLST.csv")
byCPST_sizes.to_csv(outbase + "_ByCPST.csv")
byCPS_sizes.to_csv(outbase + "_ByCPS.csv")

#print "Mapping Reads to Keep"
#yn = in1d(array(df.ix[:,"bamrow"]),bytag_want)
#df['keep_read'] = pd.Series(yn, index=df.index)
#print df

#print "Writing CSV"
#df.to_csv("collapse_check.csv")
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Step 3: Write out a filtered BAM file

# Iterate over it again, and dump out to new BAM if the ID checks out
print "Writing new BAM with just the de-duplicated reads"
samfile = pysam.AlignmentFile(args.inpath, "rb")
newbam = pysam.AlignmentFile(outbase + "_DePcrDup.bam", "wb", template=samfile)

#bamrow = []
#print "Loading data from BAM to RAM"
recno=0
for read in samfile.fetch():
	if recno in bytag_want:
		newbam.write(read)

	#newbam.write(read)

	recno += 1
	if (recno % 1000000) == 0:
		print "Done with Record Number: " + str(recno)
		sys.stdout.flush()
		#break
#print tags
#print chrs
#print poss
#print strands
#print seqs

samfile.close()
newbam.close()

log_file.close()
# --------------------------------------------------------------------

