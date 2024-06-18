# Wolfe_Stull_Copula_PNAS
The code and files associated with this repository may be used to recreate the results in *A Bayesian framework to infer ontogenetic relationships and predict associated parameters using human growth and development traits* by Christopher A. Wolfe and Kyra E. Stull. All data in the present analyses derives from the open-access repository associated with the **Subadult Virtual Anthropology Database**. The data file can be found at [Stull & Corron, 2022](https://zenodo.org/records/5193208). More information about the repository can be found at the [homepage](https://www.unr.edu/anthropology/research-and-facilities/subadult-database).
> Stull, K., & Corron, L. (2021). SVAD_US (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5193208.

> Stull, K., & Corron, L. (2022). Subadult Virtual Anthropology Database (SVAD) Data Collection Protocol: Epiphyseal Fusion, Diaphyseal Dimensions, Dental Development Stages, Vertebral Neural Canal Dimensions. Zenodo. https://doi.org/10.5281/zenodo.7293977.

The following steps allow a user to complete the analyses from the publication:
1. Execute the commands in 'data_prep.R' to prepare the data for analyses
   1. This step requires the user to import the data from Stull & Corron, 2022
2. Using the data adapted from above, sample from the log posterior density of a Gaussian copula in 'model_fit.R'
   1. Requires the stan code 'MixGaussCop_Growth.stan'
3. The remainder of the analyses can be completed in any order and utilize the posterior samples initially saved in 'model_fit.R'
   1. For Figure 3 ('imputation_Fig3.R'), a user needs the stan code 'impute1.stan' and 'impute2.stan'
   2. For Figure 4 ('age_estim_Fig4.R') a user needs the stan code 'cop_age_3var.stan'
      - To complete this analysis a user must first perform age estimation using the MCP algorithm by following the [vignette](https://rpubs.com/elainechu/mcp_vignette/) (Stull et al., 2023)
> Stull, K.E., Chu, E.Y., Corron, L.K., & Price, M.H. (2023). The Mixed Cumulative Probit: A multivariate generalization of transition analysis that accommodates variation in the shape, spread, and structure of data. Royal Society Open Science.

The code here is designed to evaluate posterior samples of human growth and development data. That said, the basic mode fitting steps can be completed for any series of data types. 

To cite usage of this code please use:
Christopher Wolfe. (2024). ChristopherAWolfe/Wolfe_Stull_Copula_PNAS: Pre-submission code (v1.0.0). Zenodo.
[![DOI](https://zenodo.org/badge/816498904.svg)](https://zenodo.org/doi/10.5281/zenodo.12037912)
