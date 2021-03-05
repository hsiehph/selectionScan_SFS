# selectionScan_SFS
A simple SNAKEMAKE pipeline for whole-genome selection scan using site-frequency-spectrum (SFS) based methods

This pipeline implements several classic single-population, SFS-based selection scan methods, including Tajima's D, Fay and Wu's H, composite-likelihood-ratio test, and pi.
The main input data is a VCF and a file that lists sample IDs for a particular population. It can also take multiple populations by modifying the config.yaml file.

References:
1) Nielsen R, Williamson S, Kim Y, Hubisz MJ, Clark AG, Bustamante C. Genomic scans for selective sweeps using SNP data. Genome Res. 2005 Nov;15(11):1566-75. doi: 10.1101/gr.4252305. PMID: 16251466; PMCID: PMC1310644.
2) Tajima F. The effect of change in population size on DNA polymorphism. Genetics. 1989 Nov;123(3):597-601. PMID: 2599369; PMCID: PMC1203832.
3) Fay JC, Wu CI. Hitchhiking under positive Darwinian selection. Genetics. 2000 Jul;155(3):1405-13. PMID: 10880498; PMCID: PMC1461156.

