
## This script calculate the summary "fraction of SNPs that have iHS scores greater than a cutoff". 
## Note that the input data is a standardized iHS score file with eight columns as follows:
## "SNPID, PhyPos, FreqofSNP, EHH0, EHH1, rawScore_iHS, standardized_iHS, isOutlier"

## Usage:
##  python iHS_summaryStat_using_selscanOutput_SNPbasedWindow.v5.py  _path_byteMap_byChrom  _fname_G2D_summary_v5  _path_std_iHS_selscan_byChrom  _iHS_cutoff  _fname_outFile
##
##

import sys, os, time, math, cPickle, numpy, dadi, re, gzip

def calculate_selscanSummary(Dict_pos_currChrom_selscan, bytemap_currChrom, WinInfo, selscanType):
	
	sorted_SNPkeys = numpy.array(sorted(Dict_pos_currChrom_selscan))
	tmp_outlist_selscan_stats = {}

	for win in sorted(WinInfo):
		max_stats = None
		derv_freq = None
		id_maxStats = None

		l_pos = int(win[1]) + 1
		r_pos = int(win[2])

		len_window = r_pos - l_pos + 1

		# generate summary statistics
		len_avail_seq = bytemap_currChrom[l_pos: r_pos+1].count(b'1')

		# skip windows with NO available bases	
		if len_avail_seq == 0:
			continue

		flag_noSNPs = 0
		# find the SNPs in the current window.
		overlap_SNPkeys = sorted_SNPkeys[numpy.where((sorted_SNPkeys >= l_pos) & (sorted_SNPkeys <= r_pos))]
		overlap_SNPkeys = list(overlap_SNPkeys)

		if overlap_SNPkeys != []:
			window_avail_SNPs = len(overlap_SNPkeys)
			count_outliers = 0
			for key in overlap_SNPkeys:
				tmp = Dict_pos_currChrom_selscan[key].strip().split()
				if tmp[-1] == "1":
					count_outliers += 1
				
				# header of xpehh: [id, pos, gpos, p1, ihh1, p2, ihh2, xpehh, normxpehh, crit]
				# header of ihh12: [id, pos, p1, ihh12, normihh12, crit]
				# header of ihs:   [id, pos, p1, ihh1, ihh2, ihs, normihs, crit] 

				if selscanType in ["ihs", "ihh12"]:
					curr_stats, curr_derv_freq, curr_stats_snpid = float(tmp[-2]), float(tmp[2]), tmp[0]
				elif selscanType == "xpehh":
					curr_stats, curr_derv_freq, curr_stats_snpid = float(tmp[-2]), float(tmp[3]), tmp[0]

				if max_stats != None:
					if abs(curr_stats) > abs(max_stats):
						max_stats = curr_stats
						derv_freq = curr_derv_freq
						id_maxStats = curr_stats_snpid
				else:
					max_stats = curr_stats
					derv_freq = curr_derv_freq
					id_maxStats = curr_stats_snpid
		else:
			flag_noSNPs = 1

		if not flag_noSNPs:
			frac_outliers = float(count_outliers)/window_avail_SNPs
		elif flag_noSNPs:
			window_avail_SNPs = 0
			id_maxStats = "NA"
			derv_freq = "NA"
			max_stats = "NA"
			count_outliers = "NA"
			frac_outliers = "NA"

		# make a list of statistics calculated for the current window
		# header = ['#chromID', 'lpos_window', 'rpos_window', 'availableSites', 'numSegSites', 'numOutlierSNVs', 'fracOutlierSNVs']
		
		tmp_outlist_selscan_stats[(win[1], win[2])] = [win[0], win[1], win[2], len_avail_seq, window_avail_SNPs, frac_outliers]

	return tmp_outlist_selscan_stats



if __name__ == '__main__':
	start_time = time.clock()

	path_bytemap_dict_chrom = os.path.abspath(sys.argv[1])
	fnames = os.listdir(path_bytemap_dict_chrom)

	dict_fnames_bytemap = {}
	for name in fnames:
		if re.search("chr[0-9]+", name):
			chromID = re.search("chr([0-9]+)", name)
			dict_fnames_bytemap[int(chromID.group(1))] = os.path.join(path_bytemap_dict_chrom, name)

	f_realDNA_g2d = gzip.open(sys.argv[2])

	dict_realDNA_g2d = {}
	for line in f_realDNA_g2d:
		if line.startswith('#'):
			continue
		tmp = line.split('\t')
		chrom = int(tmp[0][3:])
		if chrom not in dict_realDNA_g2d:
			dict_realDNA_g2d[chrom] = []
		dict_realDNA_g2d[chrom].append((tmp[0], int(tmp[1]), int(tmp[2])))

	path_selscanOutput_bychrom_target = sys.argv[3]
	path_selscanOutput_bychrom_ref = sys.argv[4]
	type_selscanStats = sys.argv[5]

	if type_selscanStats in ["ihs", "xpehh", "ihh12"]:
		header = ['#chromID', 'lpos_window', 'rpos_window', 'availableSites', 'numSegSites_target', 'f_%sOutliers_target' % type_selscanStats, 'numSegSites_ref', 'f_%sOutliers_ref' % type_selscanStats]
	else:
		exit("the 5th argument must be ihs, xpehh, or ihh12!")
	
	sys.stdout.write("\t".join(header) + "\n")
	# header of xpehh: [id, pos, gpos, p1, ihh1, p2, ihh2, xpehh, normxpehh, crit]
	# header of ihh12: [id, pos, p1, ihh12, normihh12, crit]
	# header of ihs: [id, pos, freq, ihh1, ihh2, ihs, normihs, crit] 

	fnames_target = os.listdir(path_selscanOutput_bychrom_target)
	fnames_ref = os.listdir(path_selscanOutput_bychrom_ref)
	dict_fnames_selscanOutput_target = {}
	dict_fnames_selscanOutput_ref = {}
	
	for name in fnames_target:
		if re.search(type_selscanStats, name):
			if name.endswith("norm"):
				ch = re.search('chr([0-9]+)', name)
				dict_fnames_selscanOutput_target[int(ch.group(1))] = os.path.join(path_selscanOutput_bychrom_target, name)
	for name in fnames_ref:
		if re.search(type_selscanStats, name):
			if name.endswith("norm"):
				ch = re.search('chr([0-9]+)', name)
				dict_fnames_selscanOutput_ref[int(ch.group(1))] = os.path.join(path_selscanOutput_bychrom_ref, name)

	# calculate the iHS summary Stats for each window
	for currChrom in sorted(dict_fnames_selscanOutput_target.keys()):
		f_selscan_currChrom_target = open(dict_fnames_selscanOutput_target[currChrom])
		f_selscan_currChrom_ref = open(dict_fnames_selscanOutput_ref[currChrom])
		bytemap_chrom = cPickle.load(gzip.open(dict_fnames_bytemap[currChrom]))
		winInfo_chrom = dict_realDNA_g2d[currChrom]

		# exclude SNPs that are 1) monomorphic or 2) not available in the bytemap	

		dict_pos_curr_selscan_target = {}
		dict_pos_curr_selscan_ref = {}
		for line in f_selscan_currChrom_target:
			if line.startswith("id"):
				print(line)
				continue
			try:
				pos = int(line.split()[1])
			except IndexError:
				print(line)
				exit()
			try:
				pos = int(line.split()[1])
			except ValueError:
				continue
			if bytemap_chrom[pos] == 49:
				dict_pos_curr_selscan_target[pos] = line

		for line in f_selscan_currChrom_ref:
			if line.startswith("id"):
				print(line)
				continue
			try:
				pos = int(line.split()[1])
			except ValueError:
				continue
			if bytemap_chrom[pos] == 49:
				dict_pos_curr_selscan_ref[pos] = line

		# calculate the summary statistics using standardized iHS scores
		tmp_selscan_sumStat_target = calculate_selscanSummary(dict_pos_curr_selscan_target, bytemap_chrom, winInfo_chrom, type_selscanStats)

		tmp_selscan_sumStat_ref = calculate_selscanSummary(dict_pos_curr_selscan_ref, bytemap_chrom, winInfo_chrom, type_selscanStats)
		
#		tmp_outlist_selscan_stats.append([win[0], win[1], win[2], len_avail_seq, window_avail_SNPs, frac_outliers])


		for k in sorted(tmp_selscan_sumStat_target):
			if k in tmp_selscan_sumStat_ref:
				tmp_out = tmp_selscan_sumStat_target[k] + tmp_selscan_sumStat_ref[k][4:]
			else:
				tmp_out = tmp_selscan_sumStat_target[k] + ["NA"] * 2
			sys.stdout.write("\t".join([str(i) for i in tmp_out]) + "\n")



