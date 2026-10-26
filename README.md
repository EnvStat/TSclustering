# TSclustering

This repository includes the necessary files to reproduce the results from the following paper: Shafik, N., Rodrigues, A., Sihvonen, P., and Vanhatalo, J.. Model-based clustering of ecological time series. 

Instructions for use:

For simulations:

1- "simulations_replicate.R": Simulations for 'No covariates' case scenario.

2- "simulations_cov_replicate.R": Simulations for 'with covariates' case scenario.

3- "simulations_missing_replicate.R": Simulations for 'sparse' case scenario. 

For accuracy and coverage calculations:

Run code "simulations_accuracy_coverage.R"

'cos_effort_covariates.stan': The file includes the model appropriate for full data. The model includes sampling effort and covariates as well.

Run the code 'moth_clustering_species.R": To get the data and run the model.


