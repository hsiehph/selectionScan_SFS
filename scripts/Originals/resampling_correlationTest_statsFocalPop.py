import pandas as pd
from scipy.stats import pearsonr
import os, sys, re, gzip
from itertools import product

if __name__ == "__main__":

	inputFname = sys.argv[1]
	numberRows = int(sys.argv[2])

	# index starts at 1
	targetCols = sorted([int(x)-1 for x in sys.argv[3].split(",")])
	refCols = sorted([int(x)-1 for x in sys.argv[4].split(",")])
	pairs_columns = list(product(targetCols, refCols))
	
	fileName = sys.argv[5]
	
	if inputFname.startswith("stdin"):
		f_in = sys.stdin
	else:
		f_in = gzip.open(inputFname)

	df = pd.read_table(f_in)
	columnNames = df.columns

	# sample with replacement
	df_sample_regions = df.sample(n=numberRows, replace=True)
	
	header = ["file"]
	out_list = [fileName]

	for (p1, p2) in pairs_columns:
		header.extend(["corr_%s_%s" % (columnNames[p1], columnNames[p2]), "pval_%s_%s" % (columnNames[p1], columnNames[p2])])
		df_sample_regions_clean = df_sample_regions.iloc[:,[p1, p2]].dropna()
		corr, pval = pearsonr(df_sample_regions_clean.iloc[:, 0], df_sample_regions_clean.iloc[:, 1].abs())
		out_list.extend([corr, pval])

	sys.stdout.write("\t".join(header) + "\n")
	sys.stdout.write("\t".join([str(x) for x in out_list]) + "\n")



