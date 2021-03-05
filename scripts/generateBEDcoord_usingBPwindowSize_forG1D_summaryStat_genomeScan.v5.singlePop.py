
## This script pre-calculate BED coords for windows that have a given number of variants (windowSize) from a VCF
## Sites that consist of missing data and/or are non-biallelic are excluded (i.e. # segregating alleles must be two).
## A site is excluded if 1) its ancestral allele is not one of those segregating in human samples or 2) masked in the given bytemap.
##

import re, gzip, sys, os, time, math, pickle, numpy, dadi


from collections import deque
from itertools import islice

import os

def touch(path):
	with open(path, 'a'):
		os.utime(path, None)
					
def iter_split(a, n):
	n = min(n, len(a))
	k, m = divmod(len(a), n)
	return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

def sliding_window1(iterable, size, step, fillvalue=None):
	if size < 0 or step < 1:
		raise ValueError
	it = iter(iterable)
	q = deque(islice(it, size), maxlen=size)
	if not q:
		return  # empty iterable or size == 0
	q.extend(fillvalue for _ in range(size - len(q)))  # pad to size
	while True:
		yield q  # iter() to avoid accidental outside modifications
		q.append(next(it))
		q.extend(next(it, fillvalue) for _ in range(step - 1))


def sliding_window2(iterable, size, step):
	len_iter = len(iterable)
	if len_iter < size:
		yield iterable
	else:
		for i in range(0, len_iter - size + 1, step):
			yield iterable[i: i + size]
	# for the leftover SNVs at the end of the chromosome if any
	if i+size != len_iter:
		yield iterable[i + size:]


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
			return None, None, None, None

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
		return None, None, None, None
	elif len(set(snp)) != 2:
		return None, None, None, None

	# extract the ancestral state from the SNPID
	for entry in line[7].split(";"):
		if entry.startswith("AA="):
			ancestral_call = entry.split("=")[1]
			break

	## ignore SNPs with incompatible ancestral allelic information
	if ancestral_call in ['N','-']:
		return None, None, None, None

	l_derived_call = [x for x in ref_alt if x != ancestral_call]
	set_derived_call = set(l_derived_call)
	if len(set_derived_call) == 1:
		derived_call = l_derived_call[0]
	else:
		print ("multiallleic sites, skipped")
		return None, None, None, None

	try:
		idx_ancestral_call = ref_alt.index(ancestral_call)
	except ValueError:
		print (ancestral_call, derived_call, ref_alt)
		exit()

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
	# e.g. in the case of simulations with CMP-DNS-NDL-AFR-EA-OCN, PopIndicatorArr=[0, 0, 1, 1, 1, 1]
	new_fs_key = [ fs_key[i] for i, v in enumerate(PopIndicatorArr) if v == 1]

	# ignore sites fixed in every population of interest
	if new_fs_key == [ PopSizeArr[i] for i, v in enumerate(PopIndicatorArr) if v == 1 ] or new_fs_key == [0] * len(new_fs_key):
		return None, None, None, None

	l_currGenos = [chromID, pos]
	l_currGenos.extend(snp)
	return tuple(new_fs_key), ' '.join(l_currGenos), ancestral_call, derived_call


### This func is a chromosome generator. Each time it yields a list with two elements, 
## in which the 1st element(curr_chrom) is the current chromosome ID, and the 2nd element(chrom) is a dictionary.
## In the dictionary, a key is the physical position of a SNP and, and
## the corresponding value of a key is the joint freq of the SNP returned by the func 'get_singleSNP_Freq'
def vcf_chrom_iter(vcf, PopSizeArr, PopIndicatorArr,  dict_fnames_byteMap, WindowSize, StepSize):
	vcfObj = gzip.open(vcf,'rb')
	chromID = None
	PopSizeArr_analysis = [PopSizeArr[i] for i, v in enumerate(PopIndicatorArr) if v == 1]
	dict_chrom = {}
	max_winIDX = numpy.floor(WindowSize / float(StepSize))

	for tmp_locus in vcfObj:
		locus = tmp_locus.decode("utf-8")
		if locus.startswith("#"):
			continue

		if chromID is None:
			try:
				chromID = int(locus.strip().split()[0])
			except ValueError:
				chromID = int(locus.strip().split()[0][3:])

			bytemap_chrom = pickle.load(gzip.open(dict_fnames_byteMap[chromID]))

		# calculate the joint freq of the current SNP
		locus_freq, genos, ancestral, derived = get_singleSNP_Freq(locus, PopSizeArr, PopIndicatorArr)
		if locus_freq == None:
			continue
		
		tmp = locus.strip().split()
		pos = int(tmp[1])

		# exclude SNPs that are 1) monomorphic or 2) not available in the bytemap	
		if locus_freq == tuple([0] * len(PopSizeArr_analysis)) or locus_freq == tuple(PopSizeArr_analysis):
			continue
		try:
			if bytemap_chrom[pos] != 49:
				continue
		except IndexError:
			print (chromID, pos, len(bytemap_chrom))
			exit("WTF!")

#		chrom[pos] = (locus_freq, ancestral, derived)

		winIDX = numpy.floor(pos / float(StepSize)) 
		count = 0
#		print ("pos:%s" % pos)
		while winIDX >= 0 and count < max_winIDX:
			winEND = WindowSize + winIDX * StepSize
#			print (winEND - WindowSize, winEND)
			if winEND not in dict_chrom:
				dict_chrom[winEND] = [chromID, {}]
			dict_chrom[winEND][1][pos] = (locus_freq, ancestral, derived)
			winIDX = winIDX - 1
			count = count + 1

	for endPos in sorted(dict_chrom.keys()):
#		for k in dict_chrom[endPos][1]:
#			print (k, dict_chrom[endPos][1][k])
#		exit()
		yield (endPos, dict_chrom[endPos])

		# yield every XX SNVs (XX as defined by WindowSize)
#		if len(chrom) == WindowSize:
#			yield [chromID, chrom]
			# remove the first YY SNVs in the current window (YY as defiend by StepSize)
#			curr_keys = sorted(chrom.keys())
#			for idx in range(StepSize):
#				del chrom[curr_keys[idx]]
#	if len(chrom) != WindowSize - StepSize:
#		yield [chromID, chrom]
		



if __name__ == '__main__':
	start_time = time.clock()

	path_bytemap_dict_chrom = os.path.abspath(sys.argv[1])
	fnames = os.listdir(path_bytemap_dict_chrom)

	dict_fnames_bytemap = {}
	for name in fnames:
		if name.startswith('chr'):
			chrID = int(re.search('chr([0-9]+).*', name).group(1))
			dict_fnames_bytemap[chrID] = os.path.join(path_bytemap_dict_chrom, name)

	windowSize = int(sys.argv[2])
	stepSize = int(sys.argv[3])
	vcf = os.path.abspath(sys.argv[4])

	popsize_arr = [2 * x for x in eval(sys.argv[5])]
	popIndicator_arr = eval(sys.argv[6])
	out_fname_prefix = sys.argv[7]
	num_batch = int(sys.argv[8])

	l_bedCoords = []
	l_dict_chromFreq = []

	for end_pos, window in vcf_chrom_iter(vcf, popsize_arr, popIndicator_arr, dict_fnames_bytemap, windowSize, stepSize):
		if len(l_bedCoords) == 0:
			chromID = window[0]

		l_dict_chromFreq.append(window[1].copy())
#		keys = sorted(window[1].keys())
#		l_pos, r_pos = (keys[0] - 1), keys[-1]
		l_pos = int(end_pos - windowSize + 1)
		r_pos = int(end_pos)
		l_bedCoords.append([window[0], l_pos, r_pos, len(window[1])])


	for idx, bedCoords in enumerate(iter_split(l_bedCoords, num_batch)):
		out_fnameBED = out_fname_prefix + "chr%s.batch%s.bed" % (chromID, idx)
		with open(out_fnameBED, "w") as fout:
			for bed in bedCoords:
				fout.write('\t'.join([str(x) for x in bed]) + "\n")


	for idx, lsub_dict_chromFreq in enumerate(iter_split(l_dict_chromFreq, num_batch)):
		out_chromFreqDict = out_fname_prefix + "chr%s.batch%s.chromFreqDict_pickle.gz" % (chromID, idx)
		new_dict_chromFreq = {}
		for d in lsub_dict_chromFreq:
			new_dict_chromFreq.update(d)
		with gzip.open(out_chromFreqDict, "wb") as fout:
			pickle.dump([chromID, new_dict_chromFreq], fout)

	if idx < num_batch:
		for i in range(idx+1, num_batch):
			out_fnameBED = out_fname_prefix + "chr%s.batch%s.bed" % (chromID, i)
			with open(out_fnameBED, "a") as fout:
				fout.write("")

			out_chromFreqDict = out_fname_prefix + "chr%s.batch%s.chromFreqDict_pickle.gz" % (chromID, i)
			with gzip.open(out_chromFreqDict, "wb") as fout:
				pickle.dump([], fout)
			

	print ("running time:%s" % (time.clock() - start_time))
