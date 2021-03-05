# This script is used to determine if a set of genes is enriched with highly differential SNPs compared to the rest of the genic SNPs.
# The purpose is to suggest candidates for possible polygenic selection.
# Note that a SNP is genic when it sits within a gene as well as the 1000 bp flanking sequences on the both side.
# In this version, we do not specify the types of mutations regarding their impacts on protein coding level.
#
# Usage:
#		python geneSetEnrich_v3.py  _PBS_pseudoTPEDformat  _fname_geneList  _fname_def_GSEA  _outputFile  >  logfile
#



import argparse
import sys, os, re, gzip, numpy, cPickle, itertools, gc, time, operator
import numpy as np
from scipy import stats
import readline
import rpy2.robjects as robjects
import statsmodels

## Python FDR: http://statsmodels.sourceforge.net/devel/generated/statsmodels.sandbox.stats.multicomp.multipletests.html#statsmodels.sandbox.stats.multicomp.multipletests
def multi_test_correction(pvals, method='bonferroni'):
	corrected_pvals = statsmodels.sandbox.stats.multicomp.multipletests(pvals, method)
	return list(corrected_pvals)


## multiple testing correction
## method can be any of the following:"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
## see http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html

def vcf_chrom_iter(vcf_filename):

	try:
		vcfObj = gzip.open(vcf_filename, 'rb')
		testRead = vcfObj.readline()
		flag_zip = 1
	except IOError:
		vcfObj = open(vcf_filename, 'r')
		flag_zip = 0

	if flag_zip == 1:
		vcfObj.close()
		vcfObj = gzip.open(vcf_filename, 'rb')

	chrom = {}
	curr_chrom = 0

	while True:
		line = vcfObj.readline()
		if line == "":
			yield  (curr_chrom, chrom)
			break
		if line.startswith('#'):
			continue

		l_locus = line.strip().split()
		try:
			chromID = int(l_locus[0])
		except ValueError:
			chromID = int(l_locus[0][3:])

		pos = int(l_locus[1])
		l_pbs = [float(x) for x in l_locus[5:8]] 

		if curr_chrom == chromID:
			value = [pos]
			value.extend(l_pbs)
			chrom[pos] = value
		elif curr_chrom <> chromID:
			if curr_chrom <> 0:
				yield  (curr_chrom, chrom)
			chrom = {}
			curr_chrom = chromID
			value = [pos]
			value.extend(l_pbs)
			chrom[pos] = value
		
		
def r_multi_test_correction(pvals, method='bonferroni'):
	r_p_adj = robjects.r['p.adjust']
	corrected_pvals = r_p_adj(robjects.FloatVector(pvals), method)
	return list(corrected_pvals)


def sliding_window2(iterable, size, step):
	len_iter = len(iterable)
	for i in xrange(0, len_iter - size + 1, step):
		yield iterable[i: i + size]


if __name__ == '__main__':
#	print 'The output file,', sys.argv[4], ', was generated by the following command:\n'
#	print ' '.join(sys.argv)

#	start_time = time.clock()

	parser = argparse.ArgumentParser()
	parser.add_argument("--inputPBS")
	parser.add_argument("--outputFile", default=False)
	args = parser.parse_args()

	# read in all pbs values for the entire genome
	dict_chrom_PBS_genome = {}

	for chrom, dict_pos_rawPBS in vcf_chrom_iter(args.inputPBS):
		l_tuples = []

		# v_list_PosPBS = [pos, PBS_pop1, PBS_pop2, PBS_pop3]
		for k_pos, v_list_PosPBS in dict_pos_rawPBS.iteritems():
			l_tuples.append(v_list_PosPBS)
		# sort pbs based on phy pos
		l_tuples.sort()
		# dict_chrom_PBS_genome[chrom] is a list of four tuples: pos, PBS_pop1, PBS_pop2, PBS_pop3
		dict_chrom_PBS_genome[chrom] = zip(*l_tuples)

	# get the ranks for all PBS values of the entire genome
	l_PBS_pop1_genome = []
	l_PBS_pop2_genome = []
	l_PBS_pop3_genome = []

	for chrom in sorted(dict_chrom_PBS_genome):
		l_PBS_pop1_genome.extend(dict_chrom_PBS_genome[chrom][1])
		l_PBS_pop2_genome.extend(dict_chrom_PBS_genome[chrom][2])
		l_PBS_pop3_genome.extend(dict_chrom_PBS_genome[chrom][3])

	# calc the ranks of PBS values for the entire genome
	l_PBS_pop1_genome = stats.rankdata(l_PBS_pop1_genome)
	l_PBS_pop2_genome = stats.rankdata(l_PBS_pop2_genome)
	l_PBS_pop3_genome = stats.rankdata(l_PBS_pop3_genome)

	# a list indicating the start and the end of each chromosome in the list l_PBS_popX_genome
	list_idx = [0]
	tmp_list = [0]
	for chrom in sorted(dict_chrom_PBS_genome):
		tmp_list.append(len(dict_chrom_PBS_genome[chrom][0]))
		list_idx.append(sum(tmp_list))

	# replace PBS with its rank in dict_chrom_PBS_genome
	for i, chrom in enumerate(sorted(dict_chrom_PBS_genome)):
		try:
			dict_chrom_PBS_genome[chrom][1] = l_PBS_pop1_genome[list_idx[i]: list_idx[i+1]]
			dict_chrom_PBS_genome[chrom][2] = l_PBS_pop2_genome[list_idx[i]: list_idx[i+1]]
			dict_chrom_PBS_genome[chrom][3] = l_PBS_pop3_genome[list_idx[i]: list_idx[i+1]]
		except IndexError:
			print (dict_chrom_PBS_genome.keys())
			print (chrom)
			print (list_idx)
			print ("WTF148!")

	# tie correction for ranks
	T_pop1 = stats.tiecorrect(l_PBS_pop1_genome)
	T_pop2 = stats.tiecorrect(l_PBS_pop2_genome)
	T_pop3 = stats.tiecorrect(l_PBS_pop3_genome)

	if T_pop1 == 0 or T_pop2 == 0 or T_pop3 == 0:
		raise ValueError('All numbers are identical in mannwhitneyu')

	# initialization for MW-U
	n2 = len(l_PBS_pop1_genome)

	with gzip.open(args.outputFile, 'wb') as fout:
		cPickle.dump([dict_chrom_PBS_genome, T_pop1, T_pop2, T_pop3, n2], fout)

