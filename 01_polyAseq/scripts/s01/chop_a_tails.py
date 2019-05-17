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
import itertools

parser = argparse.ArgumentParser(description="Remove A tails from polyA-seq reads.", epilog="Creates a new file with polyA trimmed.")
parser.add_argument('inpath', metavar='PATH_TO_GZIPPED_FASTQ', help='Path to .fastq.gz file to use as input.')
parser.add_argument('--outdir', dest='outdir', metavar='PATH_TO_OUTPUT_DIR', action='store', default=".", help='Path to place split .fastq.gz files in (default: current working directory)')
#parser.add_argument('--mm', dest='mm', metavar='ALLOWED_MISMATCHES', action='store', default=0, help='Number of mismatches to allow when finding the A-stretch (default: 0)')

args = parser.parse_args()

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

outbase = re.sub(r'\.fastq\.gz$', '', args.inpath)
outbase = os.path.basename(outbase)
outbase = os.path.join(args.outdir,outbase)

print "Reading FASTQ.GZ file: " + args.inpath
print "Basename for output files: " + outbase
sys.stdout.flush()

a_path = outbase + "_A.fastq.gz"
noa_path = outbase + "_noA.fastq.gz"
alla_path = outbase + "_allA.fastq.gz"

out_a = gzip.open(a_path, 'w+')
out_noa = gzip.open(noa_path, 'w+')
out_alla = gzip.open(alla_path, 'w+')

# Histogram of cut decisions: -1 means no AAAA found, 9 and 10 means AAAA found AND cut, other values means NOT cut but had that number mismatches in the 10 bp window
cut_hist = defaultdict(int)
# Histogram of size of the unique sequence tags IF A-tail was cut
tag_hist = defaultdict(int)

with gzip.open(args.inpath,'r') as f:
	record = []
	lineno = 1
	recno = 1
	for line in f:
		record.append(line)
		# work in blocks of 4
		if lineno == 4:
			seq = record[1].rstrip()
			loc = seq.find("AAAA")
			mypos = -1
			myscore = -1
			all_scores = []
			while loc >= 0:
				score = seq.count("A",loc,loc+10)
				all_scores.append(score)
				if score >= 9:
					mypos = loc
					myscore = score
					break
				else:
					loc = seq.find("AAAA",loc+1,len(seq))
			#seq = record[3].rstrip()
			if (mypos==-1)&(myscore==-1):
				out_noa.write(record[0])
				out_noa.write(record[1])
				out_noa.write(record[2])
				out_noa.write(record[3])
			else:
				#print "Found Site! Pos=" + str(mypos) + " Score=" + str(myscore)
				cut = seq[0:mypos]
				tag_hist[len(cut)] += 1
				#print seq[0:mypos] + " | " + seq[mypos:len(seq)]
				#print record[3][0:mypos] + " | " + record[3][mypos:len(seq)]
				if(len(cut)>0):
					out_a.write(record[0])
					out_a.write(cut + "\n")
					out_a.write(record[2])
					out_a.write(record[3][0:mypos] + "\n")
				else:
					out_alla.write(record[0])
					out_alla.write(record[1])
					out_alla.write(record[2])
					out_alla.write(record[3])

			if len(all_scores)==0:
				topscore = -1
			elif myscore > 0:
				topscore = myscore
			else:
				topscore = max(all_scores)
			#hist[myseq] += 1
			#print topscore
			cut_hist[topscore] += 1

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

#print cut_hist
#print tag_hist

csvpath = outbase + "_TagLenghAfterACut.csv"
print "Writing histogram: " + csvpath
sys.stdout.flush()
writer = csv.writer(open(csvpath, 'w'))
writer.writerow(["Length","Count"])
for key, value in tag_hist.items():
	writer.writerow([key, value])

csvpath = outbase + "_AsInWindow.csv"
print "Writing histogram: " + csvpath
sys.stdout.flush()
writer = csv.writer(open(csvpath, 'w'))
writer.writerow(["AsInWindow","Count"])
for key, value in cut_hist.items():
	writer.writerow([key, value])


out_a.close()
out_noa.close()
out_alla.close()

