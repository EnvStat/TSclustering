# TSclustering

This repository includes the necessary files to reproduce the results from the following paper: Shafik, N., Rodrigues, A., Sihvonen, P., and Vanhatalo, J.. Model-based clustering of ecological time series. 

### Instructions for use:

'cos_effort_covariates.stan': The file includes the model appropriate for full data. The model includes sampling effort and covariates as well.

Run the code 'moth_clustering_species.R": To get the data and run the model (moth species study case)

### For plotting all the figures in the manuscript, 
use plotting_freq_prob_all.R

### For data: 
# For species clustering: 
abundance_species_trait_30.RData has species abundance, and all needed traits.
# For site: 
Site_abundance_covariates.RData: contains the sparse abundance, effort and climatic covariates for 49 sampling sites and siteID.
habitat_covariate_for_site.Rdata: contains the habitat covariates of each site. 



### For simulations:

1- "simulations_replicate.R": Simulations for 'No covariates' case scenario.

2- "simulations_cov_replicate.R": Simulations for 'with covariates' case scenario.

3- "simulations_missing_replicate.R": Simulations for 'sparse' case scenario. 

For simulations accuracy and coverage calculations:

Run code "simulations_accuracy_coverage.R"




