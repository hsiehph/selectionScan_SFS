import re, gzip, sys, os, time, math, pickle, numpy, dadi

popsize_arr = [2 * x for x in eval(sys.argv[1])]
popIndicator_arr = eval(sys.argv[2])
popsize_arr_analysis = [popsize_arr[i] for i, v in enumerate(popIndicator_arr) if v == 1]

DIR_genomeSFS = os.path.abspath(sys.argv[3])
out_pickle = sys.argv[4]


Genome_FS_dict = {}

total_avail_SNPs = 0

for f in os.listdir(DIR_genomeSFS):
	if f.endswith("pickle"):
		genome_fs_dict, avail_SNPs = pickle.load(open(os.path.join(DIR_genomeSFS,f),"rb"))
		total_avail_SNPs += avail_SNPs
		
		for key in genome_fs_dict:
			if key not in Genome_FS_dict:
				Genome_FS_dict[key] = genome_fs_dict[key]
			else:
				Genome_FS_dict[key] += genome_fs_dict[key]

for key in Genome_FS_dict:
	tmp_freq = float(Genome_FS_dict[key]) / total_avail_SNPs
	Genome_FS_dict[key] = tmp_freq

with gzip.open(out_pickle , "wb") as fout:
	pickle.dump(Genome_FS_dict, fout)


