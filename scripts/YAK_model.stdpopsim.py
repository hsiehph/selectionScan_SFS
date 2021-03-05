import math
import msprime, stdpopsim
import numpy as np

def _siberian_3():
    id = "SiberiaYAK_NGA_3G19"
    description = "Siberian demographic model"
    long_description = """
        The Siberian model was drawn from the Model-A in Hsieh et al. 2019.
        It describes the ancestral human population in Eurasia,
        and the subsequent split between two Siberian populations - 
        the Yakut and Nganasan.
        Model parameters are the maximum likelihood values of the
        various parameters given in Table 1 of Hsieh et al.
        Note that parameters were based on a mutation rate of 2.35e-8 per
        base per generation and a generation time of 25 years.
    """
    populations = [
            stdpopsim.Population(
                id="CHB", 
                description="1000 Genomes CHB (Han Chinese in Beijing, China)",
                sampling_time=0),
            stdpopsim.Population(
                id="YAK",
                description="Yakuts from Siberia",
                sampling_time=0),
            stdpopsim.Population(
                id="NGA",
                description="Nganasans from Siberia",
                sampling_time=0)
    ]


    citations = [
        stdpopsim.Citation(
            author="Hsieh et al.",
            year=2019,
            doi="https://doi.org/10.1093/molbev/msx226",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]

    generation_time = 25

    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_a = 7396
    Nb_CHB_YAK_NGA = 939
    N_YAK = 16707
    N_NGA = 2840
    # Times are provided in years, so we convert into generations.

    T_CHB_YAK_NGA = 14052 / generation_time
    T_YAK_NGA = 9614 / generation_time

    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_growth = 0.008
    Ng_CHB = Nb_CHB_YAK_NGA / math.exp(-r_growth * T_CHB_YAK_NGA)
    Nb_YAK = Nb_CHB_YAK_NGA / math.exp(-r_growth * (T_CHB_YAK_NGA-T_YAK_NGA))
    # Migration rates during the various epochs.
    m_CHB_YAK_NGA = 3.02e-5

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=CHB, 1=YAK and 2=NGA
        # initially.
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=Ng_CHB, growth_rate=r_growth, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_YAK, metadata=populations[1].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_NGA, metadata=populations[2].asdict()
            )
        ],
        migration_matrix=[
            [0, m_CHB_YAK_NGA, m_CHB_YAK_NGA],  # noqa
            [m_CHB_YAK_NGA, 0, 0],  # noqa
            [m_CHB_YAK_NGA, 0, 0],  # noqa
        ],
        demographic_events=[
            # YAK and NGA merge into YAK_NGA with rate changes at T_YAK_NGA
            msprime.MassMigration(
                time=T_YAK_NGA, source=2, destination=1, proportion=1.0
            ),
            msprime.MigrationRateChange(time=T_YAK_NGA, rate=0),
            msprime.MigrationRateChange(time=T_YAK_NGA, rate=m_CHB_YAK_NGA, matrix_index=(0, 1)),
            msprime.MigrationRateChange(time=T_YAK_NGA, rate=m_CHB_YAK_NGA, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_YAK_NGA, initial_size=Nb_YAK, growth_rate=r_growth, population_id=1),
            # Population YAK_NGA merges into CHB at T_CHB_YAK_NGA
            msprime.MassMigration(time=T_CHB_YAK_NGA, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=T_CHB_YAK_NGA, rate=0),
            # Size changes to N_A at T_AF
            msprime.PopulationParametersChange(
                time=T_CHB_YAK_NGA, initial_size=N_a, growth_rate=0, population_id=0)
        ],
    )
#    dd = msprime.DemographyDebugger(
#        population_configurations=population_configurations,
#        migration_matrix=migration_matrix,
#        demographic_events=demographic_events)
#    dd.print_history()
_species.add_demographic_model(_siberian_3())



