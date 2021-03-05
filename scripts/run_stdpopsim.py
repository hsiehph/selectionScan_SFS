
import stdpopsim, sys, os


if __name__ == "__main__":

	contigID = sys.argv[1]
	geneMap = sys.argv[2]
	numDiploidSamples = int(sys.argv[3])
	popID = sys.argv[4]

	species = stdpopsim.get_species("HomSap")
	model = species.get_demographic_model("Melanesian_6H19")
	samples = model.get_samples(0,0,0,0,0,numDiploidSamples*2)
#	model = stdpopsim.PiecewiseConstantSize(species.population_size)
	contig = species.get_contig(contigID, genetic_map = geneMap)  # default is a flat genetic map
#	new_contig = stdpopsim.Contig(
#		mutation_rate = 2.35e-8,
#		recombination_map = contig.recombination_map,
#		genetic_map = contig.genetic_map,
#	)
	engine = stdpopsim.get_engine("msprime")
	ts = engine.simulate(model, contig, samples)

	n_dip_indv = int(ts.num_samples / 2)
	indv_names = [f"{popID}_sample{str(i)}" for i in range(n_dip_indv)]

	with sys.stdout as vcf_file:
		ts.write_vcf(vcf_file, contig_id=contigID, ploidy=2, individual_names=indv_names)


