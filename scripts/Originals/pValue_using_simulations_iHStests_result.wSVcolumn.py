## this script is for gathering the GSEA results, perSNP Fst + MWU, from the 1000 whole genome simulations
## Usage:
## 	python  GSEA_pValue_using_simulationGSEA_individualFst_result.v3.py  _option_Bonferroni_fdr  _fname_realData_GSEA_UsingIndividualFst  _path_GSEA_simulations  __out_fname  > yri_pygmy_bytemapFiltered_geneSetEnrich_using_realDNA_wcFst_v4.c5_pathways.GSEA.wFracSigOverSims.log
## 
## For "_option_Bonferroni_fdr", one can choose to use Bonferroni-p or fdr-p to run the analysis here.
## Note that the output contains three more columns in addition to the original six columns of the results of real data GSEA results
## There three are : "number_simulation_GSEA", "frac_GSEAsims_less_1e-2", "frac_GSEAsims_less_1e-10" 
## The fraction, for each set, refers to the proportion of the simulations that reach a given significant level (1e-2/1e-10)
## Also, the amount of significant sets, for each simulation, is also printed to standard output.


import sys, os, time, gzip
import scipy.stats as stats
from collections import deque

if __name__ == '__main__':

#	print '#The output file, ' +  sys.argv[3] + ', is generated by:' + ' '.join(sys.argv) + ' >' + sys.argv[3] + '.log'

	realdata_GSEA = gzip.open(sys.argv[1])

	list_path_simulation_GSEAfiles = [os.path.abspath(x) for x in sys.argv[2:-1]]
	list_simulation_GSEAfiles = []
	for path_sim in list_path_simulation_GSEAfiles:
		list_simulation_GSEAfiles.extend([os.path.join(path_sim,x) for x in os.listdir(path_sim)])
	out_fname = sys.argv[-1]

	# start calculating fractions of simulations, whose p values meet the two pvalue cutoffs: 1e-2 and 1e-10 
	dict_GSEA_vals = {}
	dict_sim_sigs = {}	
	dict_out = {}
	out_header1 = []
	idx_pval_bonf_realdata = []

	for currline in realdata_GSEA:
		if currline.startswith('chrom') or currline.startswith('#') or currline.startswith("X.chrom"):
			line = currline.strip().split()
			idx_pval_bonf_realdata.append(line.index("fracOutlierSNVs"))

			out_header1 = currline.strip() + '\tNUMsims\tfracGTobs_fracOutlierSNVs\n'
			continue

		formatted_line = []
		for x in currline.strip().split():
			try:
				formatted_line.append(float(x))
			except ValueError:
				formatted_line.append(x)

		for i in range(1,3):
			formatted_line[i] = int(formatted_line[i])
		key = (int(formatted_line[0][3:]), int(formatted_line[1]), int(formatted_line[2]))

		# the list of each key has the formate [ #sims, #sims >= obs_winPBS_pop3, #sims >= obs_fD, #sims >= G2D_pop23, #sims_BonfP_MWUpbs_pop3 <= 0.05]
		dict_GSEA_vals[key] = [ 0, 0]
		dict_out[key] = formatted_line

	# loop through all GSEA files for simulation data sets.
	for f in list_simulation_GSEAfiles:
		if not f.endswith('ihs_summaryWindows.gz'):
			continue
		dict_win_BonfP = {}
		idx_pval_bonf_sim = []

		f_obj = gzip.open(f)
		for currline in f_obj:
			if currline.startswith('chrom') or currline.startswith('#') or currline.startswith("X.chrom"):
				if idx_pval_bonf_sim == []:
					line = currline.strip().split()
					idx_pval_bonf_sim.append(line.index("f_ihsOutliers_target"))
					continue	

			formatted_line = []
			for x in currline.strip().split():
				try:
					formatted_line.append(float(x))
				except ValueError:
					formatted_line.append(x)
			
			sim_key = (int(formatted_line[0][3:]), int(formatted_line[1]), int(formatted_line[2])-1)
			if sim_key in dict_GSEA_vals:
				# for each set, count number of available simulations
				dict_GSEA_vals[sim_key][0] += 1
				for ii in range(len(idx_pval_bonf_sim)):
					try:
						if formatted_line[idx_pval_bonf_sim[ii]] >= dict_out[sim_key][idx_pval_bonf_realdata[ii]]:
							dict_GSEA_vals[sim_key][1 + ii] += 1
					except IndexError:
						print(f)
						print(ii)
						print(idx_pval_bonf_sim)
						print(idx_pval_bonf_sim[ii])
						print(idx_pval_bonf_realdata[ii])
						print(formatted_line[idx_pval_bonf_sim[ii]])
						print(dict_out[sim_key][idx_pval_bonf_realdata[ii]])
						print(sim_key)
						exit()
	# walk through sets and calculate, for each set, the fractions of significance over the simulations.
	out_list = []
	for k in dict_GSEA_vals:
		try:
			l_frac_sig = [ float(dict_GSEA_vals[k][idx]) / (dict_GSEA_vals[k][0]) for idx in range(1, len(dict_GSEA_vals[k]))]
		except ZeroDivisionError:
			if dict_GSEA_vals[k][0] == 0:
				l_frac_sig = ["NA"] * (len(dict_GSEA_vals[k]) - 1)
			else:
				exit('Error occurred! ', k, ':', dict_GSEA_vals[k])
		l_out = [dict_GSEA_vals[k][0]]
		l_out.extend(l_frac_sig)
		dict_out[k].extend(l_out)
	
	out_fobj = gzip.open(out_fname, "wb")

	out_fobj.write(out_header1)
	for key in sorted(dict_out):
		out_fobj.write('\t'.join([str(x) for x in dict_out[key]]) + '\n')

