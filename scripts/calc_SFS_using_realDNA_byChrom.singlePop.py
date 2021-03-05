
## This script calculate the G2D and summary (theta and pi) statistics. 
## Note that the input data is DNA SNP data in PLINK TPED format, i.e. each genotype consists of A,T,C, and/or G. 
## The SNPID ends with '_ancestralState'
## Sites that consist of missing data and/or are non-biallelic are excluded (i.e. # segregating alleles must be two).
## A site is excluded if its ancestral allele is not one of those segregating in human samples.
## Also, if the ancestral state is missing, the site is dropped from further analysis.
## The size of windows and the size of step are specified by users.
##
## Usage:
##	python  G2D_using_macsout_TPED.py  _bytemap  _numSNPs_window(numeric)  _numSNPs_step(numeric)  _filename_TPED  _number_pops(numeric)  _number_individuals_pop1  ...  _output_filename_G2D_summaryStat  _output_filename_G2D_summaryStat_windowSFS  _option_importExisting_genomeSFS  >  _logFilename
##
## In this version, a cPickle file '_output_filename_G2D_summaryStat_windowSFS' is generated, in which a dictionary stores all the AFS for all the windows.
## And also, the windowID is a tuple of two numbers (ex. (1,58)), where the first is for chrom, and the second is for the winID 
##
## Note that if the last argument '_option_importExisting_genomeSFS' is given, the script will not calculate the genome SFS; instead, it will load the specified genomeSFS for the G2D calculation.
## 		if the last argument is missing, the script will calculate the genome SFS.




import re, gzip, sys, os, time, math, pickle, numpy, dadi
#from bx.intervals.operations.quicksect import IntervalNode
#from bx.intervals.operations.quicksect import IntervalTree
#from bx.intervals.intersection import Interval, Intersecter


from collections import deque
from itertools import islice



### This func is used to get the joint-freq (used as a key) of a snp for the populations.
## This func is used in both genome-wide and local window calculations.
def get_singleSNP_Freq(SNPrecord, PopSizeArr, PopIndicatorArr):
	# SNP: a row from the tped file.
	# tuple_PopSizeArr: a tuple with two identical popSizeArr list: 
	line = SNPrecord.strip().split()
	chromID = line[0]
	pos = line[1]
	ref_alt = [line[3]]
	alt = line[4]

	alt = alt.split(",")
	ref_alt.extend(alt)

	# ignore INDELs
	for allele in ref_alt:
		if len(allele) != 1:
			return None, None

	snp = []
	for entry in line[9:]:
		if entry != ".":
			entry = entry.split(":")[0]
			entry = re.split("[|/]", entry)
			snp.extend(entry)
		else:
			snp.extend([".", "."])


	# ignore SNPs with missing data or non-biallelic .
	if "." in set(snp):
#		print(snp)
		return None, None
	elif len(set(snp)) != 2:
#		print(snp)
		return None, None

	# extract the ancestral state from the SNPID
	for entry in line[7].split(";"):
		if entry.startswith("AA"):
			ancestral_call = entry.split("=")[1]
			break

	## ignore SNPs with incompatible ancestral allelic information
	if ancestral_call in ['N','-']:
		return None, None

	idx_ancestral_call = ref_alt.index(ancestral_call)
	fs_key = []
	indx_array = numpy.cumsum([0] + PopSizeArr)
	list_start_indx = indx_array[:-1]
	list_end_indx = indx_array[1:]

	# create a key of joint freq. A key is a tuple, in which each entry is a digit representing the #derived allele of a given pop.
	for pop in range(len(indx_array) - 1):
		start_i = list_start_indx[pop]
		end_i = list_end_indx[pop]
		tmp_snp = snp[start_i : end_i]
		num_derived = len(tmp_snp) - tmp_snp.count(str(idx_ancestral_call))
		fs_key.append(num_derived)
	
	# get fs for populations indicated in PopIndicatorArr
	# e.g. in the case of NDL-AFR-EA-OCN, PopIndicatorArr=[1, 1, 1, 1]
	new_fs_key = [ fs_key[i] for i, v in enumerate(PopIndicatorArr) if v == 1]

	# ignore sites fixed in every population of interest
	if new_fs_key == [ PopSizeArr[i] for i, v in enumerate(PopIndicatorArr) if v == 1 ] or new_fs_key == [0] * len(new_fs_key):
		return None, None


	l_currGenos = [chromID, pos]
	l_currGenos.extend(snp)
	return tuple(new_fs_key), ' '.join(l_currGenos)


### This func is used to calculate the genome-wide joint FS
def genome_FS(vcf, popSizeArray, popIndicatorArray, dict_fnames_bytearray_availseq):

	popSizeArray_analysis = [popSizeArray[i] for i, v in enumerate(popIndicatorArray) if v == 1]

	# initialize (num_pops) + 1 dictionaries for constructing the single-pop SFSs
	genome_fs_dict = {}

	total_avail_SNPs = 0
	currChrom = 0

	if vcf == "stdin":
		vcfObj = sys.stdin
	else:
		vcfObj = gzip.open(vcf, 'rb')

	for tmp_locus in vcfObj:
		locus = tmp_locus.decode("utf-8")
		if locus.startswith("#"):
			continue

		try:
			chr_locus = int(locus.strip().split()[0])
		except ValueError:
			chr_locus = int(locus.strip().split()[0][3:])

		pos_locus = int(locus.strip().split()[1])

		# load the bytemap of the current chromosome	
		if currChrom != chr_locus:
			currChrom = chr_locus
			bytearray_availseq = pickle.load(gzip.open(dict_fnames_bytearray_availseq[currChrom]))

		# exclude sites that are not available in the bytemap
		if bytearray_availseq[pos_locus] != 49:
			continue

		# create a key for the genome-wide FS
		FS_key, genos = get_singleSNP_Freq(locus, popSizeArray, popIndicatorArray)

		if FS_key == None:
			continue

		# exclude monomorphic sites
		if FS_key == tuple([0] * len(popSizeArray_analysis)) or FS_key == tuple(popSizeArray_analysis):
			continue

		# update each SFS in the list "list_genome_fs_dict"
		if FS_key[0] not in genome_fs_dict:
			genome_fs_dict[FS_key[0]] = 1
		else:
			genome_fs_dict[FS_key[0]] += 1

		total_avail_SNPs += 1

	# calculate the (genome-wide) maximum likelihood estimates for the prob of observing a dervied allele of each freq class. 
#	for i in range(len(list_genome_fs_dict)):
#		for key in list_genome_fs_dict[i]:
#			tmp_freq = float(list_genome_fs_dict[i][key]) / total_avail_SNPs
#			list_genome_fs_dict[i][key] = tmp_freq

	return (genome_fs_dict, total_avail_SNPs)


if __name__ == '__main__':
	start_time = time.clock()

	path_bytemap_dict_chrom = os.path.abspath(sys.argv[1])
	fnames = os.listdir(path_bytemap_dict_chrom)

	dict_fnames_bytemap = {}
	for name in fnames:
		if name.startswith('chr'):
			chrID = int(re.search('chr([0-9]+).*', name).group(1))
			dict_fnames_bytemap[chrID] = os.path.join(path_bytemap_dict_chrom, name)

	vcf = sys.argv[2]


	# note that this script requires four populations in the order as the following:
	# pop1 = an archaic sample such as the NDL
	# pop2 = a distantly-related pop to pop1 (geographically or phylogenetically; say the AFR)
	# pop3 = a close-related pop to pop1 (geographically or phylogenetically; say the EA as to the OCN)
	# pop4 = the focul population of interset (e.g. OCN)

	popsize_arr = [2 * x for x in eval(sys.argv[3])]
	popIndicator_arr = eval(sys.argv[4])
	out_fname_SFS = sys.argv[5]

#	num_pops = 4

#	popsize_arr = []
#	[popsize_arr.append(int(sys.argv[i+3])*2) for i in range(num_pops)]

#	argv_idx = 2 + len(popsize_arr)
#	out_fname_SFS = sys.argv[argv_idx + 1]

	# calculate the genome-wide SFS using all of the SNPs
	tuple_numSNPs_Genome_FS_dict = genome_FS(vcf, popsize_arr, popIndicator_arr, dict_fnames_bytemap)

	pickle.dump(tuple_numSNPs_Genome_FS_dict, open(out_fname_SFS, 'wb'))

	print ('\nThe output pickled SFS file', out_fname_SFS, 'were generated using the following command:\n')
	print (' '.join(sys.argv))
	print ('running time:', time.clock() - start_time)

