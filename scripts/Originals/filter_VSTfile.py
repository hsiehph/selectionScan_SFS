import random, os, re, sys, time, gzip, vcf, pysam, copy
import numpy as np
import argparse, pathlib
import _pickle as cPickle
import tables
import glob, bisect, time
import pandas as pd
from scipy.stats import mannwhitneyu
from tables import open_file

def get_masked_base_count(mask_fh, chrom, chromStart, chromEnd):
	region_mask = mask_fh.getNode("/mask/%s" % chrom)[chromStart:chromEnd].sum(1)
	return np.count_nonzero(region_mask)

def vst_intCN(l_pop1CN, l_pop2CN):

	int_pop1CN = np.rint(pop1_genos)
	int_pop2CN = np.rint(pop2_genos)
	int_totalCN = np.append(int_pop1CN, int_pop2CN)

	V_total = np.var(int_totalCN)

	var_int_pop1CN = np.var(int_pop1CN)
	w_pop1 = float(len(int_pop1CN)) / len(int_totalCN)

	var_int_pop2CN = np.var(int_pop2CN)
	w_pop2 = float(len(int_pop2CN)) / len(int_totalCN)

	vst = (V_total - (var_int_pop1CN * w_pop1 + var_int_pop2CN * w_pop2)) / V_total

	if np.isnan(vst):
		vst = 0.0

	return vst


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--fileVST")
	parser.add_argument("--removeINV", default=True)
	parser.add_argument("--removeVariantBySingleCaller", default=True)
	parser.add_argument("--mask_track", default="/net/eichler/vol7/home/psudmant/genomes/mask_tracks/HG19-noWM-pad36",
						help="HDF5 mask track for reference genome")
	args = parser.parse_args()

	mask_fh = open_file(args.mask_track, mode="r")

	with open(args.fileVST) as f_vst:
		for line in f_vst:
			line = line.strip().split()
			
			# remove any inversion variant
			if args.removeINV:
				if "INV" in line[3]:
					continue

			# remove any variant called by a single caller
			if args.removeVariantBySingleCaller:
				l_callers = line[8].split(",")[:-1]
				if len(l_callers) == 1 & "dCGH" not in l_callers :
					continue
		
			# summarizing the VST file
			chrom, chromStart, chromEnd = [line[0], int(line[1]) , int(line[2])]
			len_locus = chromEnd - chromStart

			# get proportion of masked bases
			mask_count = get_masked_base_count(mask_fh, "chr%s" % chrom, chromStart, chromEnd)
			mask_pct = (mask_count / (chromEnd - chromStart))

			pop1_genos = [float(x) for x in line[6].split(",")[:-1]]
			pop2_genos = [float(x) for x in line[7].split(",")[:-1]]
			
			pop1_geno_median = np.median(np.array(pop1_genos))
			pop2_geno_median = np.median(np.array(pop2_genos))
			diff_geno_median = pop1_geno_median - pop2_geno_median
			
			try:
				uStat, uPval = mannwhitneyu(np.rint(pop1_genos), np.rint(pop2_genos))
			except ValueError:
				uPval = 1.0

			vst_integerCN = vst_intCN(pop1_genos, pop2_genos)

			line[0] = "chr" + line[0]
			tmp_output = ['%s' % len_locus , '%s' % mask_count , '%.4f' % mask_pct, '%.4f' % pop1_geno_median, '%.4f' % pop2_geno_median, '%.4f' % diff_geno_median, str(uPval), str(vst_integerCN)]
			line.extend(tmp_output)
			print('\t'.join(line))

