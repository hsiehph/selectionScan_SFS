# this script is for preparing the input file for iHS/XPEHH. The inputs are a VCF phased file (simulated) and a processed chimpanzee ref seq (Note only works for SINGLE chromosome).
# the genotypes will be replaced based on their states, where 0(macs) -> 1(iHS) for the ancestral state, and 1(macs) -> 0(iHS) for the derived state
# The results are sent to standard output.
#
# The next steps are using unix command 'cut' and java program 'transpose.jar' to group samples from the same population into a single file in the iHS/XPEHH format.
# the SNP map file for the iHS/XPEHH is also generated
#
# Usage:
# 	python convert_BEAGLEphasedFile_to_iHS_haplo_and_MAP_forSimData.py  _VCF_phasedFile  _RecombHotspotFile_inHapMapFormat  _col#_phyPos_recombFile  _col#_geneticPos_recombFile  _prefix_outputFiles

import time, sys, os, re, bisect, gzip

## for a given genetic map file (usually a map of single chrom), build a list of tuples for markers. The structure of the list is [(1_phy_pos,1_gen_pos), (2_phy_pos,2_gen_pos), ...]
def construct_tuple_for_genMap(genMap_filename, col_phy_pos, col_gen_pos):
	f_obj = open(genMap_filename,'r')
	list_of_tuples = []

	# skip the header line if it exists in the genetic map files
	f_obj.readline()

	prev_genPos = float(0)

	while True:
		line = f_obj.readline()
		if line == '':
			break
		if line.startswith("Chrom"):
			continue

		line = line.strip().split()

		## avoid genetic map markers, which have different physical pos, being mapped to the same genetic position.
		## for example, on chrom 1, many markers are mapped to the genPos 174.200926
		tmp_tuple = (int(line[col_phy_pos]), float(line[col_gen_pos]))
		if tmp_tuple[1] != prev_genPos:
			list_of_tuples.append(tmp_tuple)
			prev_genPos = tmp_tuple[1]

	return list_of_tuples

def interpolate_gen_pos(snp, curr_genMap):
	# case where the snp is in between two markers in the genetic map.
	# bisect.bisect_left(): 
	# locate the index position for inserting the current snp (based on phy pos) in curr_genMap[0]. If the snp phy pos exists in curr_genMap[0], the index will be BEFORE any existing entries.
	# the if-statement handles snps within (1st_genMap_marker, the last_genMap_marker] (Noted that this does not handle the case when the snp's position is equal to the 1st_genMap_marker's)
	tmp_snp = snp[:2]

	if curr_genMap[0][0] < int(snp[1]) and curr_genMap[0][-1] >= int(snp[1]):
		snp_ind_pos = bisect.bisect_left(curr_genMap[0], int(snp[1]))
		diff_phy_pos = curr_genMap[0][snp_ind_pos] - curr_genMap[0][snp_ind_pos-1]
		diff_gen_pos = curr_genMap[1][snp_ind_pos] - curr_genMap[1][snp_ind_pos-1]
		gen_dist_perBp = diff_gen_pos / diff_phy_pos
		tmp_snp.append(str(round(curr_genMap[1][snp_ind_pos-1] + gen_dist_perBp * (int(snp[1]) - curr_genMap[0][snp_ind_pos-1]), 12)))

	# case where the snp phy pos is equal to the 1st marker of genetic map, whose gen pos is always zero
	# having this elif-statement is simply because the code in the if-statement above is not able to handle this case.
	elif curr_genMap[0][0] == int(snp[1]):
		tmp_snp.append('0')

	# case where the snp phy pos is after the last marker of genetic map
	elif curr_genMap[0][0] > int(snp[1]) or curr_genMap[0][-1] < int(snp[1]):
		tmp_snp = []

	if tmp_snp != []:
		tmp_snp.extend(snp[2:])

	return tmp_snp

def gziplines(fname):
	from subprocess import Popen, PIPE
	f = Popen(['zcat', fname], stdout=PIPE)
	for line in f.stdout:
		yield line



if __name__ == '__main__':

	fname_vcf = sys.argv[1]
	if fname_vcf.startswith("stdin"):
		fobj_vcf = sys.stdin
	else:
		fobj_vcf = gzip.open(fname_vcf)

	# preparing for interpolating the genetic positions
	fname_genMap = os.path.abspath(sys.argv[2])
	col_phy_pos_genMap, col_gen_pos_genMap = int(sys.argv[3])-1, int(sys.argv[4])-1

	list_snp_tuples = construct_tuple_for_genMap(fname_genMap, col_phy_pos_genMap, col_gen_pos_genMap)
	curr_GenMap = list(zip(*list_snp_tuples))

	list_keptSampleIDs = []

	# subset the output TPED using a list of IDs (list_keptSamples) if the 5th argument is avaiable
	# otherwise, output all samples in the original VCF
	try:
		with open(sys.argv[5]) as f_IDs:
			for i in f_IDs:
				list_keptSampleIDs.append(i.strip())
	except IndexError:
		list_keptSampleIDs = []

	list_idx_keptSampleIDs = []

	out_obj_vcf = sys.stdout

	previousID = None

	for entry in fobj_vcf:
		if not entry.startswith('#'):

			if list_idx_keptSampleIDs == []:
				if list_keptSampleIDs == []:
					list_keptSampleIDs = list_sampleIDs
				for sampleID in list_keptSampleIDs:
					list_idx_keptSampleIDs.append(list_sampleIDs.index(sampleID))
			
			snp = entry.strip().split()

			# get a list of all possible alleles
			list_ref_alt = [snp[3]]
			for x in snp[4].split(","):
				list_ref_alt.append(x)

			# determine ancestral allele coding from the data
			l_info = snp[7].split(";")
			ancestral = "NA"
			for info in l_info:
				if info.startswith("AA"):
					ancestral = info.split("=")[1]
					ancestral_coding = str(list_ref_alt.index(ancestral))
		
			# if AA tag isn't avaliable, use reference allele for "0"
			if ancestral == "NA":
				ancestral = list_ref_alt[0]
				ancestral_coding = "0"

			genos = []
			for g in snp[9:]:
				g = g.split(":")[0]
				g = re.split("[|/]", g)
				genos.extend(g)

			alleles = list(set(genos))
			
			# skip variants > 2 alleles
			if len(alleles) not in [1, 2]:
				continue
			else:
				# recode biallelic variant with non-01 coding to 01-binary coding
				# note that here I force 0 and 1 are ancestral and derived, respectively.
				if len(alleles) == 1:
					if ancestral_coding not in alleles :
						new_genos = ["1" for x in genos]
						derived = list_ref_alt[int(alleles[0])]
					elif ancestral_coding in alleles:
						new_genos = ["0" for x in genos]
						try:
							list_ref_alt.remove(ancestral)
						except ValueError:
							print(list_ref_alt)
							print(ancestral)
							print(genos)
							exit()
						derived = list_ref_alt[0]
				elif len(alleles) == 2:
					if ancestral_coding not in alleles:
						continue
					elif ancestral_coding in alleles:
						new_genos = ["0" if x == ancestral_coding else "1" for x in genos]
						alleles.remove(ancestral_coding)
						derived = list_ref_alt[int(alleles[0])]
				else:
					print("WTF!")
					print(list_ref_alt)
					print(ancestral)
					print(genos)

			# For selscan, the TPED file has the format: [chromID(numeric), snvID(char), genetic_pos, physical_pos] 
			# Note that the snv ID has the format as "chr22_21346962_ref_alt"
			currID = "chr%s_%s_%s_%s" % (snp[0], snp[1], ancestral, derived)
			if currID == previousID:
				continue
			else:
				previousID = currID

			l_curr_snp = [snp[0], currID, None, snp[1]]
			for idx in sorted(list_idx_keptSampleIDs):
				l_curr_snp.extend(new_genos[idx*2: idx*2+2])

			tmp_new_map = interpolate_gen_pos([l_curr_snp[1], l_curr_snp[3], '0', '1'], curr_GenMap)
			if tmp_new_map == []:
#				sys.stderr.write('SNP removed due to not covered by the recomb map:' + entry.strip() + ';chimpState:' + chimp_base + '\n')
				continue
			l_curr_snp[2] = tmp_new_map[2]

			# append the current snp and its map info to the output lists
			out_obj_vcf.write("\t".join(l_curr_snp) + "\n")
		else:
			if entry.startswith("#CHROM"):
				list_sampleIDs = entry.strip().split()[9:]


