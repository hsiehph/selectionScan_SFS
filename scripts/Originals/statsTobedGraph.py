import random, os, re, sys, time, gzip, copy
import numpy as np
import argparse
import tables
import glob, bisect, time
import pandas as pd
from scipy.stats import mannwhitneyu
from tables import open_file


def iter_chrom(fobj):
	currChrom = None
	l_chrom = []
	for line in fobj:
		if line.startswith("#"):
			continue
		tmp = line.strip().split()
		tmp_chrom = tmp[0]
		if currChrom is None:
			currChrom = tmp_chrom
			l_chrom.append(line)
		else:
			if currChrom != tmp_chrom:
				yield l_chrom
				l_chrom = []
				currChrom = tmp_chrom
				l_chrom.append(line)
			else:
				l_chrom.append(line)

	yield l_chrom


def iter_midpointWindows(listLines):
	l_twoWins = []

	for line in listLines:
		if len(l_twoWins) != 2:
			l_twoWins.append(line)
		else:
			yield l_twoWins
			l_twoWins.pop(0)
			l_twoWins.append(line)

	yield l_twoWins


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--columnStat", help="Column index (starts at 0) of the statistic of interest in the input file")
	parser.add_argument("--inputFile", default=True)
	parser.add_argument("--outputFile", default=None)
	args = parser.parse_args()

	lout = []
	idx = int(args.columnStat)

	with gzip.open(args.inputFile) as fin:
		for list_chromlines in iter_chrom(fin):
			newStart = 0
			for twoLines in iter_midpointWindows(list_chromlines):
				if len(twoLines) == 2:
					line1, line2 = twoLines[0].strip().split(), twoLines[1].strip().split()

					# in the case of non-overlapping windows
					if int(line1[2]) < int(line2[1]):
						if newStart == 0:
							lout.append([line1[0], int(line1[1]) - 1 , int(line1[2]), line1[idx]])
						else:
							lout.append([line1[0], newStart , int(line1[2]), line1[idx]])

						newStart = int(line2[1])

					# in the case of overlapping windows
					else:
						mid1 = (float(line1[1]) + float(line1[2])) / 2
						mid2 = (float(line2[1]) + float(line2[2])) / 2
			
						mid_mid12 = int((mid1 + mid2) / 2)

						if newStart == 0:
							lout.append([line1[0], int(line1[1]) - 1 , mid_mid12, line1[idx]])
							newStart = mid_mid12
						else:
							lout.append([line1[0], newStart, mid_mid12, line1[idx]])
							newStart = mid_mid12

#					if line1[0] == "chr2":
#						if int(line1[1]) >= 109029898 and int(line1[1]) <= 109098836:
#							print (lout[-2:])


				# the case of the last interval of a chromosome
				else:
					line1 = twoLines[0].strip().split()
					lout.append([line1[0], newStart , line1[2], line1[idx]])

#	exit()
	if args.outputFile:
		with open(args.outputFile, 'w') as fout:
			for entry in lout:
				fout.write("\t".join([str(x) for x in entry]) + "\n")
	else:
		with sys.stdout as fout:
			for entry in lout:
				fout.write("\t".join([str(x) for x in entry]) + "\n")





