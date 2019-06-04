# AT1_AcetylationStoich
Data analysis to calculate acetylation stoichiometry

AT1 Stoichiometry Data Processing Workflow:

Main Goal: To calculate acetylation stoichiometry of the AT-1S113R/+ (KI) and AT-1 sTg and their wild type litter-mates. 

Raw Data Access:
The input datafiles for these scripts are too large to be uploaded to GitHub. They are available through the proteomics data repository, ProteomExchange (ID PXD014013). The input for the spectral_library_generation script can be found under ./SPECTRUM_LIBRARY/Puglielli_complete_library_20180927.csv. The inputs for the stoich_processing_AT1 scripts can be found under ./RESULT/20181015_DenuPuglielli_AT1_cyto_Spectronaut10_Report.csv and ./RESULT/20181023_DenuPuglielli_mito_Spectronaut10_stoich_Report.csv


Order of scripts:

1_spectral_library_generation:
	•	Adds a H vs L label to each fragment ion based on their m/z values to be later merged onto the data

2_stoich_processing_AT1:
	•	Takes the raw output from Spectronaut to calculate the stoichiometry of a given acetylation site (k_site) for 			a given peptide
	•	Uses an additional script (Isotopic_distribution_BRAIN.R) to correct for the natural isotopic distribution of 			the light peak
	•	Cytosolic fractions and mitochondrial fractions were processes separately, so there are two scripts that can 			be run in any order

3_stoichsummary_stats_AT1:
	•	Calculates average stoichiometry for each of the four conditions (AT-1S113R/+ (KI) and AT-1 sTg (TG) and their 			wild type litter-mates)
	•	Calculates statistical significance using one-way analysis of variance (ANOVA)

Cite as:
alawton2. (2019, June 4). alawton2/AT1_AcetylationStoich: First release AcStoich (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.3238525
