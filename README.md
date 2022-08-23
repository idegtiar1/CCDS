# Conditional Cross-Design Synthesis
This repository contains code and simulated data for "Conditional Cross-Design Synthesis Estimators for Generalizability in Medicaid" by Irina Degtiar, Tim Layton, Jacob Wallace, and Sherri Rose (2021), [arXiv:2109.13288](https://arxiv.org/abs/2109.13288). The conditional cross-design synthesis (CCDS) estimators presented in this paper combine randomized and observational data, while addressing their respective biases---lack of overlap in the randomized data and unmeasured confounding in the observational data. They thus allow for causal inference on a target population represented by a (possibly reweighted) union of the data sources. The estimators include outcome regression, propensity weighting, and double robust approaches. All use the covariate overlap between the randomized and observational data to remove potential unmeasured confounding bias

# Instructions for Use
## CCDS estimators
The R function “get_estimates_with_bootstrap” in the Real_data/Helper_files/helper_functions.R file allows users to fit the CCDS and comparison estimators presented in the paper.

## Reproducibility
The code contains all functions that were used to generate simulation results (found in Section 5 and Supplemental Materials) and Medicaid study results (found in Section 6 and Supplemental Materials). All tables and figures for the paper’s simulation results can be reproduced by running the provided code; for the real data analysis, all code to reproduce tables & figures is provided but the necessary data cannot be provided due to the data use agreement with the Centers for Medicare and Medicaid Services.

To reproduce results, run the scripts in the “Scripts” folder in numerical order:
- Simulations/Scripts
  *	0_Generate_data.R (was run locally)
  *	1a_Run_simulations.R (was run on server)
  *	1b_Run_simulations_full_ensemble.R (was run on server)
  *	2_Summarize_simulation_results.Rmd (was run locally)
  *	3_Generate_plots_for_presentation.Rmd (was run locally)
- Real_data/Scripts (requires data so cannot be run)
  *	1_Data_cleaning.R
  *	2_EDA.Rmd
  *	3_Estimates.R
  *	4_Summarize_Medicaid_results.Rmd
  *	Appendix_Accounting_for_CMY/Fixed_effects/1_Data_cleaning.R
  *	Appendix_Accounting_for_CMY/Fixed_effects/3_Estimates.R
  *	Appendix_Accounting_for_CMY/Random_effects/1_Data_cleaning.R
  *	Appendix_Accounting_for_CMY/Random_effects/3_Estimates.R

These scripts call the corresponding Helper_files to execute the code. Tables/figures corresponding to manuscript tables/figures are noted with “*PAPER TABLE/FIGURE Number* (table/figure:reference_label)” (e.g., *PAPER FIGURE 2* (figure:results_base_case)).

## Code Description
The provided R code has no licensing information. The code is hosted on the corresponding author’s github. At the time of submission, the commit reference was f26266b.

While there are no hardware nor operating system requirements per se (besides that the implemented parallelization does not work in Windows), due to the computational time needed for the bootstrap and simulation replications, we recommend running code on a cluster. For the simulation, the full set of 84 scenarios examined took 5 days to run point estimates with the ensemble consisting of a “fast” subset containing SL.glm, SL.glm.interaction, SL.earth, and SL.nnet, 5 additional days to re-run the full ensemble result point estimates, 2 days to run the two linear regression bootstraps, approximately 16 days to run the fast ensemble bootstrap, and approximately 40 days to run the full ensemble bootstrap. Results with the fast and full ensembles were similar. These runs were performed in parallel on 100 cores with an Intel(R) Xeon(R) CPU E7-8895 v2 @ 2.80GHz processor. For the Medicaid study, the analysis including bootstrap took approximately 12 days to run in parallel on 20 cores with an Intel(R) Xeon(R) CPU E5-2630 v4 @ 2.20GHz processor.

Simulations were run using R version 3.5.1 (2018-07-02). The following R packages versions were used:

| kernlab_0.9-27     |  ranger_0.10.1      |  qs_0.21.2          |  fst_0.9.4          |  magrittr_2.0.1     
| gridExtra_2.3      |  ggthemes_4.2.0     |  kableExtra_0.9.0   |  knitr_1.31         |  SuperLearner_2.0-24
| nnls_1.4           |  Hmisc_4.1-1        |  Formula_1.2-4      |  survival_2.42-3    |  lattice_0.20-35    
| glmnet_2.0-16      |  foreach_1.5.1      |  Matrix_1.2-14      |  nnet_7.3-12        |  forcats_0.5.1      
| stringr_1.4.0      |  dplyr_1.0.4        |  purrr_0.3.4        |  readr_1.4.0        |  tidyr_1.1.2        
| tibble_3.0.6       |  ggplot2_3.3.3      |  tidyverse_1.3.0    | | 

Simulated data generation and results summaries were run using R version 4.0.2 (2020-06-22). The following R package versions were used:

kernlab_0.9-29      mefa4_0.3-7         qs_0.24.1           fst_0.9.4          
magrittr_2.0.1      gridExtra_2.3       ggthemes_4.2.4      kableExtra_1.3.4   
knitr_1.31          SuperLearner_2.0-26 nnls_1.4            Hmisc_4.5-0        
Formula_1.2-4       survival_3.1-12     lattice_0.20-41     glmnet_4.1-1       
Matrix_1.2-18       nnet_7.3-14         forcats_0.5.1       stringr_1.4.0      
dplyr_1.0.5         purrr_0.3.4         readr_1.4.0         tidyr_1.1.3        
tibble_3.1.0        ggplot2_3.3.3       tidyverse_1.3.0

The code also used the pw_overlap function available from https://github.com/rachelnethery/overlap (commit reference c54975a).
