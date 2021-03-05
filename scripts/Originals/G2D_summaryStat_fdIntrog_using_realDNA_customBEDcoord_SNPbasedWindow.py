
## This script calculate the G2D and summary (theta and pi) statistics. 
## Usage:
##	python  G2D_using_macsout_TPED.py  _bytemap  _numSNPs_window(numeric)  _numSNPs_step(numeric)  _filename_TPED  _number_pops(numeric)  _number_individuals_pop1  ...  _output_filename_G2D_summaryStat  _output_filename_G2D_summaryStat_windowSFS  _option_importExisting_genomeSFS  >  _logFilename
##
## Note that if the last argument '_option_importExisting_genomeSFS' is given, the script will not calculate the genome SFS; instead, it will load the specified genomeSFS for the G2D calculation.
## 		if the last argument is missing, the script will calculate the genome SFS.




import pysam, re, gzip, sys, os, time, math, cPickle, numpy, dadi


from collections import deque
from itertools import islice


def get_singleSNP_Freq(VariantRecordSamples, PopSizeArr):
	chromID = VariantRecordSamples.chrom
	pos = VariantRecordSamples.pos
	ref_alt = [VariantRecordSamples.ref]
	ref_alt.extend(list(VariantRecordSamples.alts))
	ancestral_call = VariantRecordSamples.info["AA"][0]

	# ignore INDELs
	for allele in ref_alt:
		if len(allele) != 1:
			return None, None, None, None

	snp = []
	for s in VariantRecordSamples.samples.keys():
		if VariantRecordSamples.samples[s]['GT'] != ".":
			snp.extend([str(x) for x in VariantRecordSamples.samples[s]['GT']])
		else:
			snp.extend([".", "."])
	# ignore SNPs with missing data or non-biallelic .
	if "." in set(snp):
		return None, None, None, None
	elif len(set(snp)) != 2:
		return None, None, None, None

#	snp = [ref_alt.index(x) for x in snp]

#	idx_ancestral_call = ref_alt.index(ancestral_call)
	derived_call = [ x for x in ref_alt if x != ancestral_call][0]
	fs_key = []
	indx_array = numpy.cumsum([0] + PopSizeArr)
	list_start_indx = indx_array[:-1]
	list_end_indx = indx_array[1:]

	# create a key of joint freq. A key is a tuple, in which each entry is a digit representing the #derived allele of a given pop.
	for pop in xrange(len(indx_array) - 1):
		start_i = list_start_indx[pop]
		end_i = list_end_indx[pop]
		tmp_snp = snp[start_i : end_i]
#		num_derived = len(tmp_snp) - tmp_snp.count(str(idx_ancestral_call))
		num_derived = len(tmp_snp) - tmp_snp.count(ancestral_call)
		fs_key.append(num_derived)
	
	l_currGenos = [str(chromID), str(pos)]
	l_currGenos.extend(snp)
	return tuple(fs_key), ' '.join(l_currGenos), ancestral_call, derived_call


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


def calculate_G2D(list_GenomeFSDict, Dict_chrom_chromFreqDict, Dict_chrom_bytemap, list_bedLoci, popsizeArray):
## GenomeFSDict: each key is a tuple presenting the joint freq and its corresponding value is the fraction of SNPs with this joint freq
## chromFreqDict: a dictionary, whose keys are chromID (int) and their correponding values are lists. In each such a list, the 1st element is the current chromID(numeric), and the 2nd element is a dictionary, whose keys are positions and values are corresponding joint freq of the SNVs.
	
	tmp_outlist_stats_focalPop = []
	tmp_outlist_stats_indvPop = []
	tmp_output_dict_snp ={}
	tmp_output_dict_windowSFS ={}
	list_pos_FST_PBS = []
	
	# iteratively find the SNPs in the current window
	for line in list_bedLoci:
		currChrom, l_pos, r_pos = line
		# adjust to 1-based
		curr_chrom = int(currChrom[3:])
		l_pos = int(l_pos) + 1
		r_pos = int(r_pos)

		chromFreqDict = Dict_chrom_chromFreqDict[curr_chrom]
		bytemap_currChrom = cPickle.load(gzip.open(Dict_chrom_bytemap[curr_chrom]))

		# generate summary statistics
		len_window = r_pos - l_pos + 1

		# skip windows with NO available bases
		len_avail_seq = bytemap_currChrom[l_pos: r_pos+1].count(b'1')

		if len_avail_seq == 0:
			continue
			
		set_keys = set(chromFreqDict[1].keys())

		# ignore locus with zero SNV
		if len(set_keys) == 0:
			continue
		if len(set_keys.intersection(xrange(l_pos, r_pos+1))) == 0:
			continue

		window_avail_SNPs = 0
		window_FS_dict = {}

		# initialize a list of dictionaries, in which each dictionary is a SFS dictionary for an individual pop
		list_indpop_window_FS_dict = []
		list_indpop_window_fs_forSumStat = []
		for i in range(len(popsizeArray)):
			list_indpop_window_FS_dict.append({})
			list_indpop_window_fs_forSumStat.append(numpy.zeros(numpy.array(popsizeArray[i])+1))

		# intialize a numpy matrix for the joint frequency spectrum
		joint_fs_forSumStat_p23 = numpy.zeros(numpy.array([popsizeArray[2], popsizeArray[3]])+1)
		joint_fs_forSumStat_p13 = numpy.zeros(numpy.array([popsizeArray[1], popsizeArray[3]])+1)
		joint_fs_forSumStat_p12 = numpy.zeros(numpy.array([popsizeArray[1], popsizeArray[2]])+1)

		# intialization for calculating the fd statistic
		ABBAsum = 0.0
		BABAsum = 0.0
		maxABBAsumD = 0.0
		maxBABAsumD = 0.0

		# construct the FS (window_FS_dict) of the current window.

		for key in set_keys.intersection(xrange(l_pos, r_pos+1)):
			window_avail_SNPs += 1

			# note that this script requires four populations in the order as the following:
			# pop3 = focal pop (e.g. the OCN)
			# pop2 = a close-related pop to pop1 (geographically or phylogenetically; say the EA as to the OCN)
			# pop1 = a distantly-related pop to pop1 (geographically or phylogenetically; say the AFR)
			# pop0 = a archaic pop such as the NDL

			# update ABBA-BABA values for calculating the fD statistic of the current window
			# note since the the allele frequencies were calculated by polarizing w.r.t. chimpanzee (Clint), ABBA-BABA values were calculated using AFR-OCN-NDL only
			# Also note that on the tree, pop1_tree == p1Freq (AFR in VCF), pop2_tree == p3Freq (OCN in VCF), pop3_tree == pop0Freq (NDL in VCF)
			# That is, it is required that in the VCF, the order of populations (p0, p1, p2, and p3) is the NDL, an AFR pop, a pop closely-related to the focal, and the focal (e.g. OCN).

			try:
				p0freq, p1freq, p2freq, p3freq = numpy.array([float(x) for x in chromFreqDict[1][key][0]]) / numpy.array([float(y) for y in popsizeArray])
				pop1Freq_tree = p1freq
				pop2Freq_tree = p3freq
				pop3Freq_tree = p0freq
				pop4Freq_tree = 0  # using only one outgroup sequence to determine ancestral states)

				ABBAsum += (1 - pop1Freq_tree) * pop2Freq_tree * pop3Freq_tree * (1 - pop4Freq_tree)
				BABAsum += pop1Freq_tree * (1 - pop2Freq_tree) * pop3Freq_tree * (1 - pop4Freq_tree)
				if pop3Freq_tree >= pop2Freq_tree:
					maxABBAsumD += (1 - pop1Freq_tree) * pop3Freq_tree * pop3Freq_tree * (1 - pop4Freq_tree)
					maxBABAsumD += pop1Freq_tree * (1 - pop3Freq_tree) * pop3Freq_tree * (1 - pop4Freq_tree) 
				else:
					maxABBAsumD += (1 - pop1Freq_tree) * pop2Freq_tree * pop2Freq_tree * (1 - pop4Freq_tree)
					maxBABAsumD += pop1Freq_tree * (1 - pop2Freq_tree) * pop2Freq_tree * (1 - pop4Freq_tree)
			except:
				continue

			# calculate the pairwise FST and PBS statistics for the current SNV
			curr_perSNV_PBS = ["chr"+str(curr_chrom), key]
			curr_perSNV_PBS.extend(calc_perSNV_PBS(popsizeArray, chromFreqDict[1][key][0]))

			# estimate heterozygosity for pop3 (i.e. the focal population, say the OCN)
			hetZ_pop3 = 2 * float(p3freq) * (1.0 - p3freq) * ( float(popsizeArray[3]) / (float(popsizeArray[3]) - 1) )
			curr_perSNV_PBS.append(hetZ_pop3)

			# append derived allele freq.
			curr_perSNV_PBS.extend(chromFreqDict[1][key][0])

			# add ancestral and derived allele info
			curr_perSNV_PBS.extend(chromFreqDict[1][key][1:])
			list_pos_FST_PBS.append(curr_perSNV_PBS)


			# update window_FS_dict to make a joint SFS for pop2 and pop3 for calculating the G2D statistic
			key_jointFreq_pop23 = chromFreqDict[1][key][0][2:]
			if key_jointFreq_pop23 not in window_FS_dict:
				window_FS_dict[key_jointFreq_pop23] = 1
			else:
				window_FS_dict[key_jointFreq_pop23] += 1

			p0, p1, p2, p3 = chromFreqDict[1][key][0]
			# construct a joint FS of pop-pair  for calculate other summary statistics
			joint_fs_forSumStat_p23[[numpy.array(i) for i in [p2, p3]]] += 1
			joint_fs_forSumStat_p13[[numpy.array(i) for i in [p1, p3]]] += 1
			joint_fs_forSumStat_p12[[numpy.array(i) for i in [p1, p2]]] += 1

			# construct the individual pop SFS for each of the 4 pops
			for i in range(len(popsizeArray)):
				if chromFreqDict[1][key][0][i] not in list_indpop_window_FS_dict[i]:
					list_indpop_window_FS_dict[i][chromFreqDict[1][key][0][i]] = 1
					list_indpop_window_fs_forSumStat[i][chromFreqDict[1][key][0][i]] = 1
				else:
					list_indpop_window_FS_dict[i][chromFreqDict[1][key][0][i]] += 1
					list_indpop_window_fs_forSumStat[i][chromFreqDict[1][key][0][i]] += 1

		try:
		# max PBS of the current window
			currWin_95percPBS = numpy.percentile(zip(*list_pos_FST_PBS)[5], 95)
		except IndexError:
			print (curr_chrom, set_keys)
			print (l_pos, r_pos)
			exit("WTF_283")

		# convert the above SFS to DaDi spectrum objects
		joint_fs_forSumStat_p23 = dadi.Spectrum(joint_fs_forSumStat_p23)
		joint_fs_forSumStat_p13 = dadi.Spectrum(joint_fs_forSumStat_p13)
		joint_fs_forSumStat_p12 = dadi.Spectrum(joint_fs_forSumStat_p12)
		for i in range(len(popsizeArray)):
			list_indpop_window_fs_forSumStat[i] = dadi.Spectrum(list_indpop_window_fs_forSumStat[i])

		# calculate G2D
		list_G2Ds = []
		CL_genome = 0
		CL_window = 0

		for key in window_FS_dict:
			num_snps = window_FS_dict[key]
			try:
				CL_genome += num_snps * math.log(list_GenomeFSDict[0][key])
			except KeyError:
				continue
			tmp_freq = float(window_FS_dict[key])/window_avail_SNPs
			CL_window += num_snps * math.log(tmp_freq)

		joint_G2D = 2 * (CL_window - CL_genome)
		list_G2Ds.append(joint_G2D)

		# calculate G1D for individual population		
		for indx, indpop_GenomeSFS in enumerate(list_GenomeFSDict[1:]):
			CL_genome = 0
			CL_window = 0
			for key in list_indpop_window_FS_dict[indx]:
				num_snps = list_indpop_window_FS_dict[indx][key]
				CL_genome += num_snps * math.log(indpop_GenomeSFS[key])
				tmp_freq = float(list_indpop_window_FS_dict[indx][key])/window_avail_SNPs
				CL_window += num_snps * math.log(tmp_freq)
			indpop_G1D = 2 * (CL_window - CL_genome)
			list_G2Ds.append(indpop_G1D)

		# calculate FST and PBS for the current window
		winFst_p23 = max(0, joint_fs_forSumStat_p23.Fst())
		winFst_p13 = max(0, joint_fs_forSumStat_p13.Fst())
		winFst_p12 = max(0, joint_fs_forSumStat_p12.Fst())

		if winFst_p23 < 1:
			T_winFst_p23 = -math.log(1 - winFst_p23)
		elif winFst_p23 >= 1:
			T_winFst_p23 = -math.log(1 - 0.9999999999)

		if winFst_p13 < 1:
			T_winFst_p13 = -math.log(1 - winFst_p13)
		elif winFst_p13 >= 1:
			T_winFst_p02 = -math.log(1 - 0.9999999999)

		if winFst_p12 < 1:
			T_winFst_p12 = -math.log(1 - winFst_p12)
		elif winFst_p12 >= 1:
			T_winFst_p12 = -math.log(1 - 0.9999999999)

		winPBS_pop1 = max(0, (T_winFst_p12 + T_winFst_p13 - T_winFst_p23) / 2 )
		winPBS_pop2 = max(0, (T_winFst_p12 + T_winFst_p23 - T_winFst_p13) / 2 )
		winPBS_pop3 = max(0, (T_winFst_p23 + T_winFst_p13 - T_winFst_p12) / 2 )
		
		list_win_PBS_FST = [winPBS_pop1, winPBS_pop2, winPBS_pop3, winFst_p12, winFst_p13, winFst_p23]

		# calculate pi for each of the four pops
		list_pi= []
		for i in range(len(popsizeArray)):
			list_pi.append(list_indpop_window_fs_forSumStat[i].pi() / len_avail_seq)

		# calculate Tajima's D for each of the four pops
		list_tajima_d = []
		for i in range(len(popsizeArray)):
			list_tajima_d.append(list_indpop_window_fs_forSumStat[i].Tajima_D())

		# calculate the fD statistic for the current window
		try:
			fD = (ABBAsum - BABAsum) / (maxABBAsumD - maxBABAsumD)
		except ZeroDivisionError:
			fD = "NA"

		list_fD = [fD, ABBAsum, BABAsum]

		# output the following lists of statistics
		# list_fD, list_tajima_d, list_win_FST_PBS, list_G2Ds, 
		# make a list of statistics calculated for the current window
		list_G2Ds = ['%.6f' % elemt for elemt in list_G2Ds]
		list_win_PBS_FST = ['%.6f' % elemt if elemt != "NA" else elemt for elemt in list_win_PBS_FST]
		list_fD = ['%.6f' % elemt if elemt != "NA" else elemt for elemt in list_fD]
		list_tajima_d = ['%.6f' % elemt if elemt != "NA" else elemt for elemt in list_tajima_d]
		list_pi = ['%.6f' % elemt if elemt != "NA" else elemt for elemt in list_pi]

		curr_out_focalPop = [str(int(x)) for x in [curr_chrom, l_pos, r_pos, len_avail_seq]]
		curr_out_focalPop[0] = "chr" + curr_out_focalPop[0]
		curr_out_focalPop.append(str(window_avail_SNPs))
		curr_out_focalPop.append('%.6f' % currWin_95percPBS)
		curr_out_focalPop.extend(list_win_PBS_FST)
		curr_out_focalPop.extend(list_fD)
		curr_out_focalPop.append(list_G2Ds[0])
		tmp_outlist_stats_focalPop.append(curr_out_focalPop)

		curr_out_indvPop = [str(int(x)) for x in [curr_chrom, l_pos, r_pos, len_avail_seq]]
		curr_out_indvPop[0] = "chr" + curr_out_indvPop[0]
		curr_out_indvPop.append(str(window_avail_SNPs))
		curr_out_indvPop.extend(list_G2Ds[1:])
		curr_out_indvPop.extend(list_tajima_d)
		curr_out_indvPop.extend(list_pi)

		tmp_outlist_stats_indvPop.append(curr_out_indvPop)

		# store the SFS of pop2 and pop3 for the current window at *tmp_output_dict_windowSFS*
		tmp_output_dict_windowSFS[(curr_chrom, l_pos, r_pos, len_window)] = joint_fs_forSumStat_p23


	return tmp_outlist_stats_focalPop, tmp_outlist_stats_indvPop, tmp_output_dict_windowSFS, list_pos_FST_PBS



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
	path_vcf = os.path.abspath(sys.argv[3])
	fnames_vcf = os.listdir(path_vcf)

	dict_fnames_vcf = {}
	for name in fnames_vcf:
		if re.search('chr[0-9]+.*\.vcf\.gz$', name):
			chrID = re.search('chr([0-9]+).*\.vcf\.gz$', name)
			dict_fnames_vcf[int(chrID.group(1))] = os.path.join(path_vcf, name)


	# note that this script requires four populations in the order as the following:
	# pop0 = an archaic sample such as the NDL
	# pop1 = a distantly-related pop to pop1 (geographically or phylogenetically; say the AFR)
	# pop2 = a close-related pop to pop1 (geographically or phylogenetically; say the EA as to the OCN)
	# pop3 = the focul population of interset (e.g. OCN)

	num_pops = 4

	popsize_arr = []
	[popsize_arr.append(int(sys.argv[i+4])*2) for i in range(num_pops)]

	argv_idx = 3 + len(popsize_arr)
	out_fname_prefix = sys.argv[argv_idx + 1]

	out_fname_stat_focalPop = out_fname_prefix + ".statsFocalPop.gz"
#	out_fname_stat_indvPop = out_fname_prefix + ".statsIndvPop.gz"
#	out_fname_perSNV_PBS = out_fname_prefix + ".perSNVpbs.gz"

	try:
		list_Genome_FS_dict = cPickle.load(gzip.open(sys.argv[argv_idx+2]))
	except IndexError:
		exit("Error! This script requires an input for genome SFS from a list of pre-calculated SFS")

	# initialize output column names
	pi_labs = []
	g1d_labs = []
	tajima_labs = []
	for i in range(num_pops):
		pi_labs.append('pi_pop' + str(i))
		g1d_labs.append('G1D_pop' + str(i))
		tajima_labs.append('TajimaD_pop' + str(i))


	output_stats_focalPop_list =[['#chromID', 'lpos_window', 'rpos_window', 'availableSites', 'numSegSites', 'win95percPBS', \
								  'winPBS_pop1', 'winPBS_pop2', 'winPBS_pop3', 'winFST_pop12', 'winFST_pop13', 'winFST_pop23', \
								  'fD', 'ABBA', 'BABA', 'G2D_pop23', 'identifySegDup', 'numSegDup' ,'distanceToSegDup']]

#	output_stats_indvPop_list = [['#chromID', 'lpos_window', 'rpos_window', 'availableSites', 'numSegSites'] + \
#								  g1d_labs + tajima_labs + pi_labs ]

#	output_perSNV_FST_PBS_list = [['#chrom', 'pos', 'FST_pop23', 'FST_pop13', 'FST_pop12', 'PBS_pop1', 'PBS_pop2', 'PBS_pop3', 'hetZ_pop3', 'derivedFreq_pop0', 'derivedFreq_pop1', 'derivedFreq_pop2', 'derivedFreq_pop3', 'ancestral', 'derived']]

	# output_dict: a dictionary, whose keys are (chromID, windowID) and the corresponding values are snps in the plink tped format.

	output_dict_windowSFS = {}
	if f_windowBED.startswith("stdin"):
		f_BED = sys.stdin
	else:
		f_BED = open(f_windowBED)

	dict_chrom_BEDloci = {}
	list_BEDloci = []
	for line in f_BED:
		chrom, win_start, win_end = line.strip().split()
		chromID = int(chrom[3:])
		list_BEDloci.append([chrom, win_start, win_end])

		vcf_reader = pysam.VariantFile(dict_fnames_vcf[chromID])
		bytemap_chrom = cPickle.load(gzip.open(dict_fnames_bytemap[chromID]))

		if chromID not in dict_chrom_BEDloci:
			dict_pos_freq = {}
			dict_chrom_BEDloci[chromID] = [chromID, dict_pos_freq]

		if any(True for _ in vcf_reader.fetch(str(chromID), int(win_start), int(win_end))):
			vcf_records = vcf_reader.fetch(str(chromID), int(win_start), int(win_end))
			for record in vcf_records:
				pos = record.pos
				chromID = record.chrom
				locus_freq, genos, ancestral, derived = get_singleSNP_Freq(record, popsize_arr)
				if locus_freq is None:
					continue
				try:
					if bytemap_chrom[pos] != 49:
						continue
				except IndexError:
					print (chromID, pos, len(bytemap_chrom))
					exit("WTF!")
				dict_chrom_BEDloci[int(chromID)][1][int(pos)] = (locus_freq, ancestral, derived)
		else:
			continue

	# calculate the G2D and summary Stats for each window
	tmp_sumStat_focalPop, tmp_sumStat_indvPop, tmp_dict_windowSFS, list_perSNV_FST_PBS = calculate_G2D(list_Genome_FS_dict, dict_chrom_BEDloci, dict_fnames_bytemap, list_BEDloci, popsize_arr)

	output_stats_focalPop_list.extend(tmp_sumStat_focalPop)
#	output_stats_indvPop_list.extend(tmp_sumStat_indvPop)
#	output_perSNV_FST_PBS_list.extend([entry for entry in list_perSNV_FST_PBS if entry not in output_perSNV_FST_PBS_list])
#	output_dict_windowSFS.update(tmp_dict_windowSFS)

	# output statistics to a file specified by out_fname_g2d
	out_fobj = gzip.open(out_fname_stat_focalPop, 'w')
	for entry in output_stats_focalPop_list:
		tmp = '\t'.join([str(i) for i in entry]) + '\n'
		out_fobj.write(tmp)
	out_fobj.close()


"""
	out_fobj = gzip.open(out_fname_stat_indvPop, 'w')
	for entry in output_stats_indvPop_list:
		tmp = '\t'.join([str(i) for i in entry]) + '\n'
		out_fobj.write(tmp)
	out_fobj.close()

	# output a list of perSNV FST and PBS values
	out_fobj = gzip.open(out_fname_perSNV_PBS, 'w')
	for k, entry in enumerate(output_perSNV_FST_PBS_list):
		if k == 0:
			tmp = '\t'.join(entry) + '\n'
		else:
			tmp = '\t'.join(['%.6f' % i if (ii > 1 and ii < 9) else str(i) if (ii==0 or ii > 12) else str(int(i)) for ii,i in enumerate(entry)]) + '\n'
		out_fobj.write(tmp)
	out_fobj.close()
"""

