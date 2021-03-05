import re, gzip, sys, os, time, math, cPickle, numpy, dadi

popsize_arr = [2 * x for x in eval(sys.argv[1])]
popIndicator_arr = eval(sys.argv[2])
popsize_arr_analysis = [popsize_arr[i] for i, v in enumerate(popIndicator_arr) if v == 1]

DIR_genomeSFS = os.path.abspath(sys.argv[3])
out_pickle = sys.argv[4]


list_Genome_FS_dict = []
for i in range(len(popsize_arr_analysis)+1):
	list_Genome_FS_dict.append({})

total_avail_SNPs = 0

for f in os.listdir(DIR_genomeSFS):
	if f.endswith("pickle"):
		list_genome_fs_dict, avail_SNPs = cPickle.load(open(os.path.join(DIR_genomeSFS,f)))
		total_avail_SNPs += avail_SNPs
		
		for i in range(len(list_genome_fs_dict)):
			for key in list_genome_fs_dict[i]:
				if key not in list_Genome_FS_dict[i]:
					list_Genome_FS_dict[i][key] = list_genome_fs_dict[i][key]
				else:
					list_Genome_FS_dict[i][key] += list_genome_fs_dict[i][key]

for i in range(len(list_Genome_FS_dict)):
	for key in list_Genome_FS_dict[i]:
		tmp_freq = float(list_Genome_FS_dict[i][key]) / total_avail_SNPs
		list_Genome_FS_dict[i][key] = tmp_freq

with gzip.open(out_pickle , "wb") as fout:
	cPickle.dump(list_Genome_FS_dict, fout)


