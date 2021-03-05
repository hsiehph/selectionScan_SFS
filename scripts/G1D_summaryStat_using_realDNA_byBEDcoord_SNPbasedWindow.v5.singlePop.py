
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
## In this version, a pickle file '_output_filename_G2D_summaryStat_windowSFS' is generated, in which a dictionary stores all the AFS for all the windows.
## And also, the windowID is a tuple of two numbers (ex. (1,58)), where the first is for chrom, and the second is for the winID 
##
## Note that if the last argument '_option_importExisting_genomeSFS' is given, the script will not calculate the genome SFS; instead, it will load the specified genomeSFS for the G2D calculation.
## 		if the last argument is missing, the script will calculate the genome SFS.




import re, gzip, sys, os, time, math, pickle, numpy, dadi
from collections import deque
from itertools import islice


def calc_SinglePair_Fst(list_popsizes, list_freqs):

	list_popsizes = [float(x) for x in list_popsizes]
	list_freqs = [float(x) for x in list_freqs]

	r = float(len(list_popsizes))
	n = float(sum(list_popsizes))

	p_bar = 0.0
	for i, v in enumerate(list_freqs):
		p_bar += (list_popsizes[i] / n) * list_freqs[i]

	n_bar = sum(list_popsizes) / r

	S_sq_part1 = 1.0 / ((r - 1) * n_bar)

	S_sq_part2 = 0.0
	sum_n_i = 0.0
	sum_n_i_sq = 0.0
	for i, v in enumerate(list_freqs):
		S_sq_part2 += list_popsizes[i] * ((list_freqs[i] - p_bar)**2)
		sum_n_i += list_popsizes[i]
		sum_n_i_sq += (list_popsizes[i])**2


	S_sq = S_sq_part1 * S_sq_part2

	T1 = S_sq - (1.0 / (2*n_bar - 1)) * (p_bar * (1-p_bar) - (r-1)/r * S_sq)

	nc_part1 = 1.0 / (r - 1)
	nc_part2 = (sum_n_i - (sum_n_i_sq / sum_n_i))
	nc = nc_part1 * nc_part2

	T2_part1 = (2 * nc - 1) / (2 * n_bar - 1)
	T2_part2 = p_bar * (1- p_bar)
	T2_part3 = (1 + (2*(r-1)*(n_bar - nc)) / (2 * n_bar - 1))
	T2_part4 = S_sq / r
	T2 = T2_part1 * T2_part2 + T2_part3 * T2_part4

	try:
		locus_fst = T1 / T2
	except ZeroDivisionError:
		locus_fst = 0.0

	return locus_fst


def calc_perSNV_PBS(l_popsizeArray, l_Freq):
	# calculate per SNV pairwise FST and PBS statistics
	l_Freq_p3 = float(l_Freq[3]) / l_popsizeArray[3]	# the focul pop
	l_Freq_p2 = float(l_Freq[2]) / l_popsizeArray[2]	# a close-related pop to the focal pop
	l_Freq_p1 = float(l_Freq[1]) / l_popsizeArray[1]	# a distant-related pop to the focal pop 

	FST_p23 = max(0, calc_SinglePair_Fst([l_popsizeArray[2], l_popsizeArray[3]], [l_Freq_p2, l_Freq_p3]))
	FST_p13 = max(0, calc_SinglePair_Fst([l_popsizeArray[1], l_popsizeArray[3]], [l_Freq_p1, l_Freq_p3]))
	FST_p12 = max(0, calc_SinglePair_Fst([l_popsizeArray[1], l_popsizeArray[2]], [l_Freq_p1, l_Freq_p2]))
	

	if FST_p23 < 1:
		T_FST_p23 = -math.log(1 - FST_p23)
	elif FST_p23 >= 1:
		T_FST_p23 = -math.log(1 - 0.9999999999)

	if FST_p13 < 1:
		T_FST_p13 = -math.log(1 - FST_p13)
	elif FST_p13 >= 1:
		T_FST_p13 = -math.log(1 - 0.9999999999)

	if FST_p12 < 1:
		T_FST_p12 = -math.log(1 - FST_p12)
	elif FST_p12 >= 1:
		T_FST_p12 = -math.log(1 - 0.9999999999)

	pbs_3 = max(0, (T_FST_p23 + T_FST_p13 - T_FST_p12) / 2)
	pbs_2 = max(0, (T_FST_p23 + T_FST_p12 - T_FST_p13) / 2)
	pbs_1 = max(0, (T_FST_p12 + T_FST_p13 - T_FST_p23) / 2)

	return [FST_p23, FST_p13, FST_p12, pbs_1, pbs_2, pbs_3]


def calculate_G2D(GenomeFSDict, chromFreqDict, bytemap_currChrom, fname_windowBED, popsizeArray, popIndicatorArray):
## GenomeFSDict: each key is a tuple presenting the joint freq and its corresponding value is the fraction of SNPs with this joint freq
## chromFreqDict: a list, in which 1st element is the current chromID(numeric), and the 2nd element is a dictionary. 
##                The key of the 2nd is a SNP's physical positions, and its corresponding value is the joint freq of the SNP.
	
	popsizeArray_analysis = [popsizeArray[i] for i, v in enumerate(popIndicatorArray) if v == 1]

	curr_chrom = chromFreqDict[0]
	sorted_SNPkeys = sorted(chromFreqDict[1])
	tmp_outlist_stats_indvPop = []
	tmp_output_dict_snp ={}
	
	set_keys = set(chromFreqDict[1].keys())
	# iteratively find the SNPs in the current window
	fBED = open(fname_windowBED)
	for line in fBED:
			
		currChrom, l_pos, r_pos = line.strip().split()[:3]
		# adjust to 1-based
		l_pos = int(l_pos) + 1
		r_pos = int(r_pos)

		# generate summary statistics
		len_window = r_pos - l_pos + 1
		# skip windows with NO available bases	
		len_avail_seq = bytemap_currChrom[l_pos: r_pos+1].count(b'1')
		if len_avail_seq == 0:
			continue
		
		window_avail_SNPs = 0
		window_FS_dict = {}

		# initialize a list of dictionaries, in which each dictionary is a SFS dictionary for an individual pop
		indpop_window_FS_dict = {}
		indpop_window_fs_forSumStat = numpy.zeros(numpy.array(popsizeArray_analysis[0])+1)

		# construct the FS (window_FS_dict) of the current window.

		for key in set_keys.intersection(range(l_pos, r_pos+1)):
			window_avail_SNPs += 1

			# construct the individual pop SFS for each of the 4 pops
			for i in range(len(popsizeArray_analysis)):
				if chromFreqDict[1][key][0][i] not in indpop_window_FS_dict:
					indpop_window_FS_dict[chromFreqDict[1][key][0][i]] = 1
					indpop_window_fs_forSumStat[chromFreqDict[1][key][0][i]] = 1
				else:
					indpop_window_FS_dict[chromFreqDict[1][key][0][i]] += 1
					indpop_window_fs_forSumStat[chromFreqDict[1][key][0][i]] += 1


		# convert the above SFS to DaDi spectrum objects
		indpop_window_fs_forSumStat = dadi.Spectrum(indpop_window_fs_forSumStat)

		# calculate G1D for individual population		
		CL_genome = 0
		CL_window = 0
		for key in indpop_window_FS_dict:
			num_snps = indpop_window_FS_dict[key]
			CL_genome += num_snps * math.log(GenomeFSDict[key])
			tmp_freq = float(indpop_window_FS_dict[key])/window_avail_SNPs
			CL_window += num_snps * math.log(tmp_freq)
		indpop_G1D = 2 * (CL_window - CL_genome)

		# calculate pi for each of the four pops
		pi = indpop_window_fs_forSumStat.pi() / len_avail_seq

		# calculate Tajima's D for each of the four pops
		tajima_d = indpop_window_fs_forSumStat.Tajima_D()

		# calculate Fay and Wu's H for each of the four pops
		theta_L = indpop_window_fs_forSumStat.theta_L()

		num = pi * len_avail_seq - theta_L 
		n = indpop_window_fs_forSumStat.sample_sizes[0]
		an = numpy.sum(1./numpy.arange(1,n))
		bn = numpy.sum(1./numpy.arange(1,n)**2)
		bn1 = numpy.sum(1./numpy.arange(1,n+1)**2)
		theta  = indpop_window_fs_forSumStat.Watterson_theta()
		s = indpop_window_fs_forSumStat.S()
		theta_sq = s*(s-1.)/(an**2 + bn)
		var = (n-2.)*theta / (6*(n-1.)) + (18*n**2*(3*n+2.)*bn1 - (88*n**3+9*n**2-13*n+6.)) * theta_sq / (9*n*(n-1.)**2)
		h = 2 * num
		h_norm = num / numpy.sqrt(var)

		# output the following lists of statistics
		# list_fD, list_tajima_d, list_win_FST_PBS, list_G2Ds, 
		# make a list of statistics calculated for the current window
		try:
			G1D = '%.6f' % indpop_G1D
			Tajima_d = '%.6f' % tajima_d
			PI = '%.6f' % pi
			H = '%.6f' % h
			H_norm = '%.6f' % h_norm

			curr_out_indvPop = [str(int(x)) for x in [curr_chrom, l_pos-1, r_pos, len_avail_seq]]
			curr_out_indvPop[0] = "chr" + curr_out_indvPop[0]
			curr_out_indvPop.append(str(window_avail_SNPs))
			curr_out_indvPop.append(G1D)
			curr_out_indvPop.append(Tajima_d)
			curr_out_indvPop.append(H)
			curr_out_indvPop.append(H_norm)
			curr_out_indvPop.append(PI)
		except TypeError:
			curr_out_indvPop = [str(int(x)) for x in [curr_chrom, l_pos-1, r_pos, len_avail_seq]]
			curr_out_indvPop[0] = "chr" + curr_out_indvPop[0]
			curr_out_indvPop.append(str(window_avail_SNPs))
			curr_out_indvPop.extend(["NA"] * 5)

		tmp_outlist_stats_indvPop.append(curr_out_indvPop)

	return tmp_outlist_stats_indvPop



if __name__ == '__main__':
	start_time = time.clock()

	path_bytemap_dict_chrom = os.path.abspath(sys.argv[1])
	fnames = os.listdir(path_bytemap_dict_chrom)

	dict_fnames_bytemap = {}
	for name in fnames:
		if name.startswith('chr'):
			chrID = int(re.search('chr([0-9]+).*', name).group(1))
			dict_fnames_bytemap[chrID] = os.path.join(path_bytemap_dict_chrom, name)

	f_windowBED = sys.argv[2]
	f_pickle_chromFreqDict = sys.argv[3]

	# note that this script requires four populations in the order as the following:
	# pop0 = an archaic sample such as the NDL
	# pop1 = a distantly-related pop to pop1 (geographically or phylogenetically; say the AFR)
	# pop2 = a close-related pop to pop1 (geographically or phylogenetically; say the EA as to the OCN)
	# pop3 = the focul population of interset (e.g. OCN)

	popsize_arr = [2 * x for x in eval(sys.argv[4])]
	popIndicator_arr = eval(sys.argv[5])
	popsize_arr_analysis = [popsize_arr[i] for i, v in enumerate(popIndicator_arr) if v == 1]
	out_fname_prefix = sys.argv[6]

	try:
		Genome_FS_dict = pickle.load(gzip.open(sys.argv[7]))
	except IndexError:
		exit("Error! This script requires an input for genome SFS from a dict of pre-calculated SFS")

	out_fname_stat_indvPop = out_fname_prefix + ".statsIndvPop.gz"

	# initialize output column names
	pi_labs = []
	g1d_labs = []
	tajima_labs = []
	H_labs = []
	H_norm_labs = []
	for i in range(len(popsize_arr_analysis)):
		g1d_labs.append('G1D_pop' + str(i))
		tajima_labs.append('TajimaD_pop' + str(i))
		H_labs.append('H_pop' + str(i))
		H_norm_labs.append('H_norm_pop' + str(i))
		pi_labs.append('pi_pop' + str(i))

	output_stats_indvPop_list = [['#chromID', 'lpos_window', 'rpos_window', 'availableSites', 'numSegSites'] + \
								  g1d_labs + tajima_labs + H_labs + pi_labs ]

	# output_dict: a dictionary, whose keys are (chromID, windowID) and the corresponding values are snps in the plink tped format.

	if not os.path.isfile("statsIndvPop.header"):	
		with open("statsIndvPop.header", "w") as fout:
			tmp = ['#chromID', 'lpos_window', 'rpos_window', 'availableSites', 'numSegSites'] + g1d_labs + tajima_labs + H_labs + H_norm_labs + pi_labs
			tmp.extend(["numGenes","genes"])
			fout.write('\t'.join(tmp) + "\n")

	with gzip.open(f_pickle_chromFreqDict) as fpickle:
		chromosome = pickle.load(fpickle)

	try:
		bytemap_chrom = pickle.load(gzip.open(dict_fnames_bytemap[int(chromosome[0])]))
	except IndexError:
		print("Error: possible empty file %s!" % f_windowBED)
		out_fobj = gzip.open(out_fname_stat_indvPop, 'w')
		out_fobj.write("")
		out_fobj.close()
		exit()

	# calculate the G2D and summary Stats for each window
	tmp_sumStat_indvPop = calculate_G2D(Genome_FS_dict, chromosome, bytemap_chrom, f_windowBED, popsize_arr, popIndicator_arr)

	output_stats_indvPop_list.extend(tmp_sumStat_indvPop)


	# output statistics to a file specified by out_fname_g2d
	out_fobj = gzip.open(out_fname_stat_indvPop, 'wb')
	for entry in output_stats_indvPop_list:
		tmp = '\t'.join([str(i) for i in entry]) + '\n'
		out_fobj.write(tmp.encode())
	out_fobj.close()


#	pickle.dump(list_Genome_FS_dict, open(out_fname_prefix + '_genomeFSDict.pickle', 'wb'))
#	pickle.dump(output_dict_windowSFS, gzip.open(out_fname_prefix + '_windowSFS.pickle.gz', 'wb'))

#	print ('\nThe output files', out_fname_stat_focalPop, 'and', out_fname_stat_indvPop, 'were generated using the following command:\n')
#	print (' '.join(sys.argv))
#	print ('running time:', time.clock() - start_time)

