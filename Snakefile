import sys,os,numpy,re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

DIRS_TO_MAKE = ["log", "genomeFS_byChrom/"]
for folder in DIRS_TO_MAKE:
	if not os.path.exists(folder):
		os.makedirs(folder)

configfile: "config.yaml"

DIR_byteMap = config["DIR_byteMap"]
#outFname_prefix = config["outFname_prefix"]
windowSize = config["windowSize"]
stepSize = config["stepSize"]
numBatch = config["numBatch"]
l_suffix = config["list_suffix"]
l_batchID = list(range(int(numBatch)))

list_numSamples_indicator = "\"" + str(config["list_numSamples_indicator"]) + "\""

list_chroms = ["chr%s" % i for i in range(1,23)]
#list_chroms = ["chr7"]

rule dummy:
	input:	
			["selectionStats/%s/%s/%s.statsIndvPop.gz.tbi" % (group, pop, pop) for group in config["VCF"].keys() for pop in config["superPop"][group].keys() ] ,
			["genomeFS_byChrom/%s/%s/%s.genomeFS.pickle.gz" % (group, pop, pop) for group in config["VCF"].keys() for pop in config["superPop"][group].keys() ] ,
#			expand("selectionStats/{group}/{pop}/{pop}.statsIndvPop.gz.tbi", group=config["VCF"].keys(), pop=config["superPop"][wildcards.group].keys()),
#			expand("genomeFS_byChrom/{group}/{pop}/{pop}.genomeFS.pickle.gz", group=config["VCF"].keys(), pop=config["superPop"][wildcards.group].keys())

rule tabix:
	input:	"selectionStats/{group}/{pop}/{pop}.statsIndvPop.gz" ,
	output: "selectionStats/{group}/{pop}/{pop}.statsIndvPop.gz.tbi"
	params: sge_opts = "-l mfree=2G -l h_rt=24:00:00"
	priority: 19
	run:
		shell("""  tabix -p vcf {input} """)


rule merge_selectionStats:
	input: expand("selectionStats_byBED/{{group}}/{{pop}}/{{pop}}.{CHROM}.batch{batchID}.statsIndvPop.gz", CHROM=list_chroms, batchID=l_batchID)
	output:	"selectionStats/{group}/{pop}/{pop}.statsIndvPop.gz",
	params: sge_opts = "-l mfree=10G -l h_rt=24:00:00"
	priority: 25
	run:
		shell(""" zcat {input} | grep -v "#" | sort -k1,1 -k2,2n | bedtools intersect -a - -b /net/eichler/vol27/projects/human_diversity/nobackups/hsiehph/genomicData/AnnoHumanGenes/refSeq.hg38.uniq.bed -loj | awk '{{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14}}' | bedtools groupby -full -g 1-3 -c 11,11 -o count,collapse | cut --complement -f11 | cat ./statsIndvPop.header - | bgzip -c > {output} """)



rule selectionStats_byBED:
	input:	"winBEDs_DictFreq/{group}/{pop}/{pop}.{CHROM}.batch{batchID}.bed",
			"winBEDs_DictFreq/{group}/{pop}/{pop}.{CHROM}.batch{batchID}.chromFreqDict_pickle.gz",
			"genomeFS_byChrom/{group}/{pop}/{pop}.genomeFS.pickle.gz"
	output: 
			temp("selectionStats_byBED/{group}/{pop}/{pop}.{CHROM}.batch{batchID}.statsIndvPop.gz")
	params: sge_opts = "-l mfree=2G -l h_rt=24:00:00"
	priority: 40
	run:
		file_sampleID = config["superPop"][wildcards.group][wildcards.pop]
		with open(file_sampleID) as fin:
			num_lines = sum(1 for line in fin if line.rstrip()) 
		list_numSamples_sim = "\"[" + str(num_lines) + "]\""
		shell(
		""" python scripts/G1D_summaryStat_using_realDNA_byBEDcoord_SNPbasedWindow.v5.singlePop.py {DIR_byteMap} {input[0]}  {input[1]} {list_numSamples_sim} {list_numSamples_indicator}  $TMPDIR/{wildcards.pop}.{wildcards.CHROM}.batch{wildcards.batchID} {input[2]} ; """
		""" rsync --exclude="$TMPDIR/.[!.]*"  --bwlimit=100000 $TMPDIR/{wildcards.pop}.{wildcards.CHROM}.batch{wildcards.batchID}.statsIndvPop.gz  selectionStats_byBED/{wildcards.group}/{wildcards.pop} ; """)


rule winBED_DictFreq:
	input: "superPopVCF/{group}/{pop}/{pop}.{CHROM}.vcf.gz"
	output: 
			temp(["winBEDs_DictFreq/{group}/{pop}/{pop}.{CHROM}.batch%s.bed" % batchID for batchID in l_batchID]),
			temp(["winBEDs_DictFreq/{group}/{pop}/{pop}.{CHROM}.batch%s.chromFreqDict_pickle.gz" % batchID for batchID in l_batchID])
	params: sge_opts = "-l mfree=4G -l h_rt=24:00:00"
	priority: 45
	run:
		file_sampleID = config["superPop"][wildcards.group][wildcards.pop]
		with open(file_sampleID) as fin:
			num_lines = sum(1 for line in fin if line.rstrip()) 
		list_numSamples_sim = "\"[" + str(num_lines) + "]\""
		shell(
		""" python scripts/generateBEDcoord_usingBPwindowSize_forG1D_summaryStat_genomeScan.v5.singlePop.py  {DIR_byteMap}  {windowSize}  {stepSize}  {input}  {list_numSamples_sim} {list_numSamples_indicator}  $TMPDIR/{wildcards.pop}.  {numBatch} ; """
		""" rsync --exclude="$TMPDIR/.[!.]*" --bwlimit=100000 $TMPDIR/{wildcards.pop}.{wildcards.CHROM}.batch*.bed  winBEDs_DictFreq/{wildcards.group}/{wildcards.pop}/ ; """
		""" rsync --exclude="$TMPDIR/.[!.]*" --bwlimit=100000 $TMPDIR/{wildcards.pop}.{wildcards.CHROM}.batch*.chromFreqDict_pickle.gz  winBEDs_DictFreq/{wildcards.group}/{wildcards.pop}/ """)


rule sum_GenomeFS:
	input: expand("genomeFS_byChrom/{{group}}/{{pop}}/{{pop}}.{CHROM}.genomeFS.pickle", CHROM=list_chroms)
	output: "genomeFS_byChrom/{group}/{pop}/{pop}.genomeFS.pickle.gz"
	params: sge_opts= "-l mfree=4G -l h_rt=12:00:00", DIR_genomeFS_byChrom="genomeFS_byChrom/{group}/{pop}/"
	priority: 45
	run:
		file_sampleID = config["superPop"][wildcards.group][wildcards.pop]
		with open(file_sampleID) as fin:
			num_lines = sum(1 for line in fin if line.rstrip()) 
		list_numSamples_sim = "\"[" + str(num_lines) + "]\""
		shell(""" python scripts/sumALLpickle_genomeFS.singlePop.py {list_numSamples_sim}  {list_numSamples_indicator}  {params.DIR_genomeFS_byChrom}  {output} """)


rule calc_GenomeFS:
	input: "superPopVCF/{group}/{pop}/{pop}.{CHROM}.vcf.gz"
	output: "genomeFS_byChrom/{group}/{pop}/{pop}.{CHROM}.genomeFS.pickle"
	params: sge_opts = "-l mfree=4G -l h_rt=24:00:00"
	priority: 50
	run:
		file_sampleID = config["superPop"][wildcards.group][wildcards.pop]
		with open(file_sampleID) as fin:
			num_lines = sum(1 for line in fin if line.rstrip()) 
		list_numSamples_sim = "\"[" + str(num_lines) + "]\""
		shell(""" python scripts/calc_SFS_using_realDNA_byChrom.singlePop.py {DIR_byteMap} {input}  {list_numSamples_sim} {list_numSamples_indicator}  {output} """)


rule subsetPop:
	input: lambda wc: config["VCF"][wc.group][wc.CHROM]
	output: "superPopVCF/{group}/{pop}/{pop}.{CHROM}.vcf.gz"
	params: sge_opts = "-l mfree=4G -l h_rt=48:00:00 -pe serial 4"
	priority: 50
	run:
		file_sampleID = config["superPop"][wildcards.group][wildcards.pop]
		shell(""" bcftools view --threads 4 -S %s --force-samples -f .,PASS -M2 -v snps -Ov {input} | bgzip -c > {output} ; tabix -p vcf {output} """ % file_sampleID)




