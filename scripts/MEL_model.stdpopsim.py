import math
import msprime, stdpopsim
import numpy as np

def _melanesian_6():
    id = "Melanesian_6H19"
    description = "Melanesian demographic model"
    long_description = """
        This Melanesian model was from Hsieh et al. 2019.
        Parameters are based on Table S7 and Figure S8 (hominin only).
        It describes the hominin evolution without archaic admixture events.
        Models for the three modern lienages, including Africans, East Asians,
        and Melanesians, were inferred using Dadi (Gutenkunst et al. 2009),
        while parameters for the three archaic lineages, including 
        the Denisovan and the two Neanderthals from Altai and Vindija,
        were drwan from other studies (Malaspinas et al. 2016, Prufer 
        et al. 2017).Note that parameters were based on a mutation rate of 
        1.5e-8 per base per generation and a generation time of 29 years.
    """
    citations = [
        stdpopsim.Citation(
            author="Hsieh et al.",
            year=2019,
            doi="https://doi.org/10.1126/science.aax2083",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]

    populations = [
            stdpopsim.Population(
                id="DNS", 
                description="Denisovan, Altai, Russia",
                sampling_time=1775),
            stdpopsim.Population(
                id="NDL_Altai",
                description="Neanderthal, Altai, Russia",
                sampling_time=2252),
            stdpopsim.Population(
                id="NDL_Vindija",
                description="Neanderthal, Vindija, Croatia",
                sampling_time=1715),
            stdpopsim.Population(
                id="AFR", 
                description="African populations (SGDP)",
                sampling_time=0),
            stdpopsim.Population(
                id="EA",
                description="East Asian population (SGDP)",
                sampling_time=0),
            stdpopsim.Population(
                id="MEL",
                description="Melanesian population (SGDP)",
                sampling_time=0)
    ]

    generation_time = 29

    # First we set out the maximum likelihood values of the various parameters
    # given in Table S7.
    N_hominin = 32671
    N_AMH = 13601
    N_NDL_Altai = 826
    N_NDL_Vindija = 13249
    N_DNS = 5083
    N_AMH_exp = 27847
    N_AFR = 56022
    N_Eurasian = 679
    N_MEL = 4103
    # Times are provided in years, so we convert into generations.

    Ts_ARC_AMH = 586525 / generation_time
    Ts_DNS_NDL = 437610 / generation_time
    Ts_AMH_exp = 360940 / generation_time
    Ts_Altai_Vindija = 97875 / generation_time
    Ts_AFR_Eurasian = 74260 / generation_time
    Ts_EA_MEL = 52021 / generation_time

    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_growth = 0.0016
    N_EA = N_Eurasian / math.exp(-r_growth * Ts_AFR_Eurasian)
    # Migration rates during the various epochs.
    m_AFR_Eurasian = 0.00015
    m_AFR_EA = 1.33e-5
    m_AFR_MEL = 1.8e-6
    m_EA_MEL = 3.05e-5
    m_MEL_EA = 8.23e-5

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,

        # Population IDs correspond to their indexes in the population configuration array. 
        # Therefore, we initially have 0=DNS, 1=NDL_Altai, 2=DNL_Vindija, 3=AFR, 4=EA, 5=MEL.
        # Sampling times are assuming 29 years per generation.

        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_DNS, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_NDL_Altai, metadata=populations[1].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_NDL_Vindija, metadata=populations[2].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_AFR, metadata=populations[3].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EA, growth_rate=r_growth, metadata=populations[4].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_MEL, metadata=populations[5].asdict()
            )
        ],
        migration_matrix=[
            [0, 0, 0, 0, 0, 0],  # noqa
            [0, 0, 0, 0, 0, 0],  # noqa
            [0, 0, 0, 0, 0, 0],  # noqa
            [0, 0, 0, 0, m_AFR_EA, m_AFR_MEL],  # noqa
            [0, 0, 0, m_AFR_EA, 0, m_EA_MEL],  # noqa
            [0, 0, 0, m_AFR_MEL, m_MEL_EA, 0],  # noqa
        ],
        demographic_events=[
            # MEL and EA merge into Eurasian with rate changes at Ts_EA_MEL
            msprime.MassMigration(
                time=Ts_EA_MEL, source=5, destination=4, proportion=1.0
            ),
            msprime.MigrationRateChange(time=Ts_EA_MEL, rate=0),
            msprime.MigrationRateChange(time=Ts_EA_MEL, rate=m_AFR_Eurasian, matrix_index=(4, 3)),
            msprime.MigrationRateChange(time=Ts_EA_MEL, rate=m_AFR_Eurasian, matrix_index=(3, 4)),

            msprime.PopulationParametersChange(
                time=Ts_AFR_Eurasian, initial_size=N_Eurasian, growth_rate=r_growth, population_id=4),

            # Population Eurasian merges into AFR at Ts_AFR_Eurasian
            msprime.PopulationParametersChange(
                time=Ts_AFR_Eurasian, initial_size=N_Eurasian, growth_rate=0, population_id=4),
            msprime.MassMigration(time=Ts_AFR_Eurasian, source=4, destination=3, proportion=1.0),
            msprime.MigrationRateChange(time=Ts_AFR_Eurasian, rate=0),

            # Size changes N_AFR to N_AMH_exp  at Ts_AFR_Eurasian
            msprime.PopulationParametersChange(
                time=Ts_AFR_Eurasian, initial_size=N_AMH_exp, growth_rate=0, population_id=3),

            # Altai Neanderthal merges to Vindija at Ts_AFR_Eurasian
            msprime.MassMigration(time=Ts_Altai_Vindija, source=1, destination=2, proportion=1.0),

            # Size changes N_AMH_exp to N_AMH at Ts_AMH_exp
            msprime.PopulationParametersChange(
                time=Ts_AMH_exp, initial_size=N_AMH, growth_rate=0, population_id=3),

            # Neanderthal and Denisovan merge at Ts_DNS_NDL
            msprime.MassMigration(time=Ts_DNS_NDL, source=0, destination=2, proportion=1.0),

            # AMH and archaic merge at Ts_ARC_AMH
            msprime.MassMigration(time=Ts_ARC_AMH, source=3, destination=2, proportion=1.0),
            # Size changes N_Vindija to N_hominin at Ts_ARC_AMH
            msprime.PopulationParametersChange(
                time=Ts_ARC_AMH, initial_size=N_hominin, growth_rate=0, population_id=2)

        ],
    )
#    dd = msprime.DemographyDebugger(
#        population_configurations=population_configurations,
#        migration_matrix=migration_matrix,
#        demographic_events=demographic_events)
#    dd.print_history()

_species.add_demographic_model(_melanesian_6())



