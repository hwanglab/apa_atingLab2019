#!/usr/bin/env python

import sys
import difflib
import gzip
import re
import os
import argparse
from numpy import array
from collections import defaultdict
import csv
import sys
#import distance
import itertools

def hamming1(str1, str2):
	 return sum(itertools.imap(str.__ne__, str1, str2))

def mymatcher(mycode):
	# return difflib.SequenceMatcher(None,myseq,mycode).ratio()
	# if myseq==mycode:
	#	return 1
	# else:
	#	return 0
	return hamming1(mycode, myseq)

parser = argparse.ArgumentParser(description="Split a FASTQ file based on barcodes at start of read.", epilog="Creates a new file for each barcode.")

parser.add_argument('-i','--inpath', dest='inpath', metavar='PATH_TO_GZIPPED_FASTQ', action='store', required=True, help='Path to .fastq(.gz) file to use as input.')

parser.add_argument('-o','--outdir', dest='outdir', metavar='PATH_TO_OUTPUT_DIR', action='store', default="", help='Path to place split .fastq.gz files in (default: current working directory)')
parser.add_argument('--mm', dest='mm', metavar='ALLOWED_MISMATCHES', action='store', default=0, help='Number of mismatches to allow when matching barcodes (default: 0)')

parser.add_argument('--barcode_fn', dest='barcode_fn', metavar='BARCODES', action='store', required=True, help='Specify a file containing barcodes')

args = parser.parse_args()

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)
if not os.path.exists(args.barcode_fn):
	raise IOError('barcode file [%s] does not exist'%args.barcode_fn)

outbase = re.sub(r'\.fastq', '', args.inpath)
outbase = re.sub(r'\.gz$', '', outbase)
outbase = os.path.basename(outbase)
outbase = os.path.join(args.outdir,outbase)

print "Reading FASTQ.GZ file: " + args.inpath
print "Basename for output files: " + outbase
sys.stdout.flush()

#barcodes = array(['TGACCA','GCCAAT','CTTGTA','AAGTGC'])
fp = open(args.barcode_fn,'r')
barcodes = [i.strip() for i in fp]
fp.close()
barcodes = array(barcodes)

outpaths = []

for b in barcodes:
	b = outbase + "_" + b + ".fastq.gz"
	outpaths.append(b)

#print "Paths for split FASTQs:"
#print outpaths

outs = []
for b in outpaths:
	outs.append(gzip.open(b,'w+'))
unk_path = outbase + "_" + "OTHER" + ".fastq.gz"
out_unk = gzip.open(unk_path, 'w+')
outs = array(outs)

# Set up for histogram
unk_hist = defaultdict(int)
outs_hist = []
for b in outs:
	outs_hist.append(defaultdict(int))
outs_hist = array(outs_hist)
#print outs_hist
record = []

if (args.inpath.endswith('.gz')):
	f = gzip.open(args.inpath,'r')
else:
	f = open(args.inpath, 'r')

lineno = 1
recno = 1
for line in f:
	record.append(line)
	# work in blocks of 4
	if lineno == 4:
		#print "got record"
		myseq = record[1][0:6]
		#print myseq
		scores = array(map(mymatcher,barcodes))
		#print scores
		#print min(scores)

		# Select the match we want
		matched = barcodes[scores==min(scores)]
		if(min(scores)>int(args.mm)):
			matched = array([])

		# Prep the output
		myflank = record[1][6:12]
		thename = record[0]
		theseq = record[1][12:len(record[1])]
		thesep = record[2]
		thequals = record[3][12:len(record[3])]

		thename = re.sub(r' ','_',thename)
		thename = thename.rstrip("\n")
		thename = thename + "_" + myseq + "_" + myflank + "\n"

		if matched.size > 1:
			raise Exception('Matched more than one barcode!')
		elif matched.size == 0:
			#print "Matched no barcodes!"
			unk_hist[myseq] += 1
			out_unk.write(thename)
			out_unk.write(theseq)
			out_unk.write(thesep)
			out_unk.write(thequals)
		elif matched.size == 1:
			#print "Matched Barcode: " + matched[0] + " for " + myseq + " (mismatches=" + str(min(scores)) + ")"
			outs[scores==min(scores)][0].write(thename)
			outs[scores==min(scores)][0].write(theseq)
			outs[scores==min(scores)][0].write(thesep)
			outs[scores==min(scores)][0].write(thequals)
			outs_hist[scores==min(scores)][0][myseq] += 1
		# reset iteration
		record = []
		lineno = 0
		if (recno % 1000000) == 0:
			print "Done with Record Number: " + str(recno)
			sys.stdout.flush()
		recno += 1
	lineno += 1
f.close()

for o in outs:
	o.close()
out_unk.close()

# Write the histograms
#print unk_hist
#print outs_hist

csvpath = outbase + "_BarcodeCounts.csv"
print "Writing histograms: " + csvpath
sys.stdout.flush()
writer = csv.writer(open(csvpath, 'w'))
writer.writerow(["WantedBarcode","IncludedBarcode","Count"])
for b in barcodes:
	#print b
	myhist=outs_hist[barcodes==b][0]
	for key, value in myhist.items():
		writer.writerow([b, key, value])

for key, value in unk_hist.items():
	writer.writerow(["OTHER", key, value])
