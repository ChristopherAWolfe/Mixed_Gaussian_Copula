# Summary
The code and files associated with this repository may be used to recreate the results in *A Bayesian framework to infer ontogenetic relationships and predict associated parameters using human growth and development traits* by Christopher A. Wolfe and Kyra E. Stull. All data in the present analyses derives from the open-access repository associated with the **Subadult Virtual Anthropology Database**. The data file can be found at [Stull & Corron, 2022](https://zenodo.org/records/5193208). More information about the repository can be found at the [homepage](https://www.unr.edu/anthropology/research-and-facilities/subadult-database).
> Stull, K., & Corron, L. (2021). SVAD_US (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5193208.

> Stull, K., & Corron, L. (2022). Subadult Virtual Anthropology Database (SVAD) Data Collection Protocol: Epiphyseal Fusion, Diaphyseal Dimensions, Dental Development Stages, Vertebral Neural Canal Dimensions. Zenodo. https://doi.org/10.5281/zenodo.7293977.

The following steps allow a user to complete the analyses from the publication:
1. Execute the commands in 'data_prep.R' to prepare the data for analyses
   1. This step requires the user to import the data from Stull & Corron, 2022
2. Using the data adapted from above, sample from the log posterior density of a Gaussian copula in 'model_fit.R'
   1. Requires the stan code 'MixGaussCop_Growth.stan'
3. The remainder of the analyses can be completed in any order and utilize the posterior samples initially saved in 'model_fit.R'
4. To make the correlation plots use "correlation_plots.R".
5. To complete eigendecomposition use "eigen_analyses.R".
6. To complete the imputation task use "imputation.R".
7. To complete the age estimation task use "age_estim.R".

The code here is designed to evaluate posterior samples of human growth and development data. That said, the basic mode fitting steps can be completed for any series of data types.
