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
import regex

parser = argparse.ArgumentParser(description="Remove A tails from polyA-seq reads.", epilog="Creates a new file with polyA trimmed.")
parser.add_argument('inpath', metavar='PATH_TO_GZIPPED_FASTQ', help='Path to .fastq.gz file to use as input.')
parser.add_argument('--outdir', dest='outdir', metavar='PATH_TO_OUTPUT_DIR', action='store', default=".", help='Path to place split .fastq.gz files in (default: current working directory)')
parser.add_argument('--mm', dest='mm', metavar='ALLOWED_MISMATCHES', action='store', default=0, help='Number of mismatches to allow when finding the A-stretch (default: 0)')

args = parser.parse_args()

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

outbase = re.sub(r'\.fastq\.gz$', '', args.inpath)
outbase = os.path.basename(outbase)
outbase = os.path.join(args.outdir,outbase)

print "Reading FASTQ.GZ file: " + args.inpath
print "Basename for output files: " + outbase
sys.stdout.flush()

alen = []
hist = defaultdict(int)

with gzip.open(args.inpath,'r') as f:
	record = []
	lineno = 1
	recno = 1
	for line in f:
		record.append(line)
		# work in blocks of 4
		if lineno == 4:
			#print "got record"
			#print record[1]
			seq = record[1].rstrip()
			find = regex.findall(r'(A*){s<=' + str(args.mm) + '}', seq)
			#print find
			lens = map(len,find)
			mymax = max(lens)
			#alen.append(mymax)
			hist[mymax] += 1

			# reset iteration
			record = []
			lineno = 0
			if (recno % 1000000) == 0:
				print "Done with Record Number: " + str(recno)
				sys.stdout.flush()
			recno += 1
		lineno += 1
#print alen
#print hist

csvpath = outbase + "_LongestARunCounts.csv"
print "Writing histograms: " + csvpath
sys.stdout.flush()
writer = csv.writer(open(csvpath, 'w'))
writer.writerow(["LongestARunInSequence","Count"])
for key, value in hist.items():
	writer.writerow([key, value])


