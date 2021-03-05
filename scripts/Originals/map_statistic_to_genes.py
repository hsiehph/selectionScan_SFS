
### This script finds the overlaps between genes (in the refSeq_Genes.txt from UCSC) and windows of selection scan.
### Usage:  python  map_statistic_to_genes.py  refSeq_Genes.txt  __File_selectionScan_statistics  __output_File

### refSeq_Genes.txt: 6 columns: #name	chrom	strand	txStart	txEnd	name2 
### __File_selectionScanP_Statistics can be either G2D or XPCLR results

### G2D output with p_value
### 7 columans: chorm,  windowID,  window_lpos,  window_rpos,  number_simulation,  G2D,  p-value


### XPCLR output
### chr#	grid#	#ofSNPs_in_window	physical_pos	genetic_pos		XPCLR_score		max_s


import sys, csv, time, interval

if __name__ == '__main__':
	start_time = time.clock()
	
	# using refSeq_Genes.txt downloaded from the UCSC table RefSeq Genes track
	fobj_gene = open(sys.argv[1], 'rb')
	fobj_stats = open(sys.argv[2], 'rb')

	# build a lookup table of gene names along with their related information.
	# A key is a chromID(numeric) and its corresponding value is a list of elements: [Start_pos, End_pos, nc_accession, gene, strand, CDS]; the corresponding values are the CCDS information of the genes.
	genes_dict = {}
	for gene in fobj_gene:
		if gene.startswith('#'):
			continue
		gene = gene.strip().split()
		try:
			key = int(gene[1][3:])
		except ValueError:
			continue

		tmp_interval = interval.Interval(int(gene[3])-100000+1, int(gene[4])+100000+1)
		if key not in genes_dict:
			genes_dict[key] = []
#			genes_dict[key].append((gene, tmp_interval))
			genes_dict[key].append((gene[5], tmp_interval))
		else:
#			genes_dict[key].append((gene, tmp_interval))
			genes_dict[key].append((gene[5], tmp_interval))
			

	# loop through the windows of selection scan to find possible overlaps with genes.
	list_output_windows = []
	for window in fobj_stats:
		if window.startswith('#'):
			if window.startswith('#chr'):
				list_output_windows.append([window.strip() + '\tGene\n'])
			else:
				list_output_windows.append([window])
			continue
	
		window = window.strip().split('\t')
		key = int(window[0])
		try:
			window_range = interval.Interval(int(float(window[2])), int(float(window[3])))
		except ValueError:
			print window
			exit()
		current_chrom_genes = genes_dict[key]
		tmp_gene = []
		for gene in current_chrom_genes:
			if window_range.overlaps(gene[1]):
#				window.extend(gene[0])
				if gene[0] not in tmp_gene:
					tmp_gene.append(gene[0])
		window.append(','.join(tmp_gene))
		list_output_windows.append(window)

	fobj_output_windows = open(sys.argv[3],'wb')
	for window in list_output_windows:
		if len(window) == 1:
			fobj_output_windows.write(window)
		else:
			fobj_output_windows.write('\t'.join(window) + '\n')

	print 'running time:', time.clock() - start_time
