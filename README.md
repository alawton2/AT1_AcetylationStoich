# AT1_AcetylationStoich
Data analysis to calculate acetylation stoichiometry

AT1 Stoichiometry Data Processing Workflow:

Main Goal: To calculate acetylation stoichiometry of the AT-1S113R/+ (KI) and AT-1 sTg and their wild type litter-mates. 

Order of scripts:

1_spectral_library_generation:
	•	Adds a H vs L label to each fragment ion based on their m/z values to be later merged onto the data

2_stoich_processing_AT1:
	•	Takes the raw output from Spectronaut to calculate the stoichiometry of a given acetylation site (k_site) for a given peptide
	•	Uses an additional script (Isotopic_distribution_BRAIN.R) to correct for the natural isotopic distribution of the light peak
	•	Cytosolic fractions and mitochondrial fractions were processes separately, so there are two scripts that can be run in any order

3_stoichsummary_stats_AT1:
	•	Calculates average stoichiometry for each of the four conditions (AT-1S113R/+ (KI) and AT-1 sTg (TG) and their wild type litter-mates)
	•	Calculates statistical significance using one-way analysis of variance (ANOVA)
