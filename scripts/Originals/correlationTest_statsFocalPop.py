import pandas as pd
from scipy.stats import pearsonr
import os, sys, re, gzip
from itertools import product
import numpy as np

if __name__ == "__main__":

	inputFname = sys.argv[1]
	# index starts at 1
	targetCols = sorted([int(x)-1 for x in sys.argv[2].split(",")])
	refCols = sorted([int(x)-1 for x in sys.argv[3].split(",")])
	pairs_columns = list(product(targetCols, refCols))
	
	fileName = sys.argv[4]

	if inputFname.startswith("stdin"):
		f_in = sys.stdin
	else:
		f_in = gzip.open(inputFname)

	df = pd.read_table(f_in, header=None)
	columnNames = df.columns

	header = ["file"]
	out_list = [fileName]

	for (p1, p2) in pairs_columns:
		header.extend(["corr_%s_%s" % (columnNames[p1], columnNames[p2]), "pval_%s_%s" % (columnNames[p1], columnNames[p2])])
		df_clean = df.iloc[:,[p1, p2]].dropna()

		# check if any non-numeric element in df_clean
		# this happens when two segdups are in the same distance from the locus
		tuple_idx_nonNumeric = np.where(np.array(df_clean.applymap(lambda x: isinstance(x, (int,float)))) == False)
		if np.array(tuple_idx_nonNumeric).size != 0:
			for idx in zip(*tuple_idx_nonNumeric):
				currEntry = df_clean.iloc[idx]
				df_clean.iloc[idx] = int(currEntry.split(",")[0])

		corr, pval = pearsonr(df_clean.iloc[:, 0], df_clean.iloc[:, 1].abs())
		out_list.extend([corr, pval])

	sys.stdout.write("\t".join(header) + "\n")
	sys.stdout.write("\t".join([str(x) for x in out_list]) + "\n")


