#################################################################################################################################################
### Name
### 3/17/21
### Medicaid data analysis with conditional cross-design synthesis (CCDS) estimators
### Obtain study and target sample estimates of treatment-specific mean log total spending for Medicaid patients on different health insurance plans 
### Study: use TMLE and AIPW to obtain study population estimates
### Target: use propensity-model based, outcome model based, and double robust CCDS estimators of target population estimates
#################################################################################################################################################
print(Sys.time())

### Load data & libraries
rm(list=ls())
path = getwd()

library(qs) # read data in qs format
library(survey)  
library(SuperLearner) 
library(tmle) 
library(Hmisc)
library(pbapply)
library(parallel)
library(data.table) # fast data manipulation
library(dplyr) # data manipulation
library(madr) # model-averaged double-robust estimator
library(profvis) # memory usage
library(tictoc) # runtime
library(tidyverse) # for data manipulation
library(magrittr) # for data manipulation (set_colnames)  
library(kernlab) # for ksvm modeling
library(ranger) # for random forest modeling
library(tmle) # for TMLE sample estimates
library(BART) # for multinomial BART modeling
library(dbarts) # for BART modeling

data_to_use = "real" # or data_to_use = "simulated"

if(data_to_use == "simulated"){
  setwd("/Code/Real data/New_Medicaid_data")
  source("../Helper_files/pw_overlap.R") # for determining regions of overlap and non-overlap, from Nethery et al. 2019: https://github.com/rachelnethery/overlap
  source("../Helper_files/helper_functions.R") # for fitting models and simulating data
} else {
  source("/Helper_files/pw_overlap.R") # for determining regions of overlap and non-overlap, from Nethery et al. 2019: https://github.com/rachelnethery/overlap
  source("/Helper_files/helper_functions.R") # for fitting models and simulating data
}
tic("Overall")

#################################################################################################################################################
## Set hyperparameters
#################################################################################################################################################
### Create library of algorithms
# Fast library 
SL.glmnet.alt = create.Learner("SL.glmnet", tune = list(alpha = 0.5)) # Default alpha = 1, nfolds = 10, nlambda = 100
my_SL_library = c("SL.glm", SL.glmnet.alt$names, "SL.gam", "SL.nnet")
file_suffix = "_fast"

# Set up cluster
no_cores = 20

# Function parameters across modeling approaches
SL_library=my_SL_library # specifies SuperLearner library for complex_fit_models="ensemble"
overlap_function=F # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
X_interaction = NULL # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression 
propensity_trim_threshold = 0.001 # threshold at which to trim propensity scores and their products used in weight denominators
m=500 # number of bootstrap iterations, should be multiple of n_save
n_save = 100 # Save every n_save iterations
n_obs_select = "All" # label in saved file names indicating proportion of obs data used
random_seed=123 # seed

#################################################################################################################################################
### Read in & process data
#################################################################################################################################################
############################################
### Load data 
if(data_to_use == "simulated"){
  # Testing with simulated data
  target_sample = qread("/Data - simulated/target_sample_relevant_sim.qs")
  
  # Make smaller for testing
  target_sample = target_sample[1:5000,] %>% filter(plan_baseline_clean %in% c("A","B","C"))
  
  # Fix log spending not actually being logged
  target_sample$lpay0tot_0to5 = log(target_sample$lpay0tot_0to5 + 1)
  
  # Separate rand and obs data
  rand_sample = target_sample %>% filter(aa_sample == 2)
  obs_sample = target_sample %>% filter(aa_sample == 1)
} else {
  target_sample = qread("/data/NYS_Medicaid/Name/target_sample_relevant.qs") # Restricted to those with 6 mo of observations (can measure outcome) -- will use for analysis

  # Separate rand and obs data
  rand_sample = target_sample %>% filter(aa_sample == 1)
  obs_sample = target_sample %>% filter(aa_sample == 0)
}


# Calculate proportion of target sample made up of obs data: will need to weight obs component of estimate to account for using random sample of obs observations
# Target sample is a 10% random sample of the obs data; these numbers are from the full dataset; they were provided by Julia
# n_rand = 67619 
# n_obs = 894884 
w_obs = 894884/(67619+894884) # n_obs/(n_rand + n_obs)

# Relevant covariates
X_vars = c("age_clean120", #"age_scaled",
           "female", #"race", 
           "neighborhood",
           "aid_group_clean120", "ssi",
           "baseline_decile_clean", "baseline_decile_missing", "percent_poverty") # age, sex, neighborhood, aid group, SSI eligibility, baseline spending decile, indicator of missing baseline spending, neighborhood povery rate


### Randomized data
# Create dummy variables for plan names
A_rand = subset(rand_sample,select="plan_baseline_clean")
colnames(A_rand) = "A"

# Extract relevant variables
Y_rand = rand_sample$lpay0tot_0to5 # log total annual spending over 6 months post-randomization
X_rand = as.data.frame(model.matrix(~.,subset(rand_sample,select=X_vars)))[,-1]

## Remove spaces from variable names
colnames(X_rand) = gsub(" ", "_", colnames(X_rand))
colnames(X_rand) = gsub("-", "_", colnames(X_rand))


### Observational data
# Create dummy variables for plan names
A_obs = subset(obs_sample,select="plan_baseline_clean")
colnames(A_obs) = "A"

# Extract relevant variables
Y_obs = obs_sample$lpay0tot_0to5 # log total annual spending 1 month post-obsomization
X_obs = as.data.frame(model.matrix(~.,subset(obs_sample,select=X_vars)))[,-1]

## Remove spaces from variable names
colnames(X_obs) = gsub(" ", "_", colnames(X_obs))
colnames(X_obs) = gsub("-", "_", colnames(X_obs))

# Clean up environment
rm(rand_sample, 
   obs_sample)
gc()

#################################################################################################################################################
## Point estimates
#################################################################################################################################################
tic("Estimates - linear")
estimates_linear = get_estimates(Y_rand,
                                 A_rand, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                 X_rand,
                                 Y_obs,
                                 A_obs,
                                 X_obs,
                                 X_interaction, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression
                                 weight_obs = w_obs, # survey weight for obs data to reflect that we're using random subset of observational data but all randomized data
                                 complex_fit_models = F, # whether to use ensemble ("ensemble"), random forest ("ranger") or ksvm ("ksvm") to fit models for Y, A, and S. Note: ksvm will give error if after subsetting by treatment, some X_matrix columns just have 1 value
                                 SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                 S_ridge = F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                 overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                 propensity_trim_threshold,
                                 include_plot = T)

saveRDS(estimates_linear,
        paste0("/Results/Generalizability_results/estimates_linear_nobs",n_obs_select,".rds"))

toc()

tic("Estimates - ensemble")
estimates_ensemble = get_estimates(Y_rand,
                                 A_rand, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                 X_rand,
                                 Y_obs,
                                 A_obs,
                                 X_obs,
                                 X_interaction, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression
                                 weight_obs = w_obs, # survey weight for obs data to reflect that we're using random subset of observational data but all randomized data
                                 complex_fit_models = "ensemble", # whether to use ensemble ("ensemble"), random forest ("ranger") or ksvm ("ksvm") to fit models for Y, A, and S. Note: ksvm will give error if after subsetting by treatment, some X_matrix columns just have 1 value
                                 SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                 S_ridge = F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                 overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                 propensity_trim_threshold,
                                 include_plot = T)

saveRDS(estimates_ensemble,
        paste0("/Results/Generalizability_results/estimates_ensemble",file_suffix,"_nobs",n_obs_select,".rds"))


toc()

tic("Estimates - linear overlap2")
estimates_linear_overlap2 = get_estimates(Y_rand,
                                 A_rand, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                 X_rand,
                                 Y_obs,
                                 A_obs,
                                 X_obs,
                                 X_interaction, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression
                                 weight_obs = w_obs, # survey weight for obs data to reflect that we're using random subset of observational data but all randomized data
                                 complex_fit_models = F, # whether to use ensemble ("ensemble"), random forest ("ranger") or ksvm ("ksvm") to fit models for Y, A, and S. Note: ksvm will give error if after subsetting by treatment, some X_matrix columns just have 1 value
                                 SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                 S_ridge = F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                 overlap_function = "overlap2", # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                 propensity_trim_threshold,
                                 include_plot = T)

saveRDS(estimates_linear_overlap2,
        paste0("/Results/Generalizability_results/estimates_linear_overlap2_nobs",n_obs_select,".rds"))

toc()


#################################################################################################################################################
## Bootstrap CI
#################################################################################################################################################
# Run bootstrap, saving every 100 iterations
n_loops = ceiling(m/n_save) # Number of loops to run

# Open "cluster"
cl=makeCluster(no_cores, outfile="",type="FORK")

#################################################################################################################################################
tic("Bootstrap - linear")

# Run m bootstrap iterations, saving every n_save iterations
for(l in 1:n_loops){
  tic(paste0("Bootstrap - linear ",l))
  clusterSetRNGStream(cl, iseed = (random_seed+l))
  boot_linear_l = get_estimates_with_bootstrap(Y_rand,
                                               A_rand, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                               X_rand,
                                               Y_obs,
                                               A_obs,
                                               X_obs,
                                               X_interaction, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression
                                               weight_obs = w_obs, # survey weight for obs data to reflect that we're using random subset of observational data but all randomized data
                                               complex_fit_models = F, # whether to use linear ("ensemble"), random forest ("ranger") or ksvm ("ksvm") to fit models for Y, A, and S. Note: ksvm will give error if after subsetting by treatment, some X_matrix columns just have 1 value
                                               SL_library, # specifies SuperLearner library for complex_fit_models="linear"
                                               S_ridge = F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                               overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                               propensity_trim_threshold,
                                               M=n_save)
  saveRDS(boot_linear_l,
          paste0("/Results/Generalizability_results/boot_linear_",l,"_nobs",n_obs_select,"_M",m,".rds"))
  
  
  toc()
  rm(boot_linear_l);gc()
}

# Load data and combine bootstrap iterations
boot_linear = NULL
for(l in 1:n_loops){
  boot_linear = rbind(boot_linear, readRDS(paste0("/Results/Generalizability_results/boot_linear_",
                                                  l,file_suffix,"_nobs",n_obs_select,"_M",m,".rds")))
}

saveRDS(boot_linear,
        paste0("/Results/Generalizability_results/boot_linear","_nobs",n_obs_select,"_M",m,".rds"))
rm(boot_linear);gc()

toc()


#################################################################################################################################################
tic("Bootstrap - ensemble")

# Run m bootstrap iterations, saving every n_save iterations
for(l in 1:n_loops){
  tic(paste0("Bootstrap - ensemble ",l))
  clusterSetRNGStream(cl, iseed = (random_seed+l))
  boot_ensemble_l = get_estimates_with_bootstrap(Y_rand,
                                                 A_rand, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                                 X_rand,
                                                 Y_obs,
                                                 A_obs,
                                                 X_obs,
                                                 X_interaction, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression
                                                 weight_obs = w_obs, # survey weight for obs data to reflect that we're using random subset of observational data but all randomized data
                                                 complex_fit_models = "ensemble", # whether to use ensemble ("ensemble"), random forest ("ranger") or ksvm ("ksvm") to fit models for Y, A, and S. Note: ksvm will give error if after subsetting by treatment, some X_matrix columns just have 1 value
                                                 SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                                 S_ridge = F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                                 overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                                 propensity_trim_threshold,
                                                 M=n_save)
  saveRDS(boot_ensemble_l,
          paste0("/Results/Generalizability_results/boot_ensemble_",l,file_suffix,"_nobs",n_obs_select,"_M",m,".rds"))
  
  
  toc()
  rm(boot_ensemble_l);gc()
}

# Load data and combine bootstrap iterations
boot_ensemble = NULL
for(l in 1:n_loops){
  boot_ensemble = rbind(boot_ensemble, readRDS(paste0("/Results/Generalizability_results/boot_ensemble_",
                                                      l,file_suffix,"_nobs",n_obs_select,"_M",m,".rds")))
}

saveRDS(boot_ensemble,
        paste0("/Results/Generalizability_results/boot_ensemble",file_suffix,"_nobs",n_obs_select,"_M",m,".rds"))
rm(boot_ensemble);gc()

toc()



#################################################################################################################################################
tic("Bootstrap - linear overlap2")

# Run m bootstrap iterations, saving every n_save iterations
for(l in 1:n_loops){
  tic(paste0("Bootstrap - linear overlap2",l))
  clusterSetRNGStream(cl, iseed = (random_seed+l))
  boot_linear_overlap2_l = get_estimates_with_bootstrap(Y_rand,
                                                        A_rand, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                                        X_rand,
                                                        Y_obs,
                                                        A_obs,
                                                        X_obs,
                                                        X_interaction, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression
                                                        weight_obs = w_obs, # survey weight for obs data to reflect that we're using random subset of observational data but all randomized data
                                                        complex_fit_models = F, # whether to use linear ("ensemble"), random forest ("ranger") or ksvm ("ksvm") to fit models for Y, A, and S. Note: ksvm will give error if after subsetting by treatment, some X_matrix columns just have 1 value
                                                        SL_library, # specifies SuperLearner library for complex_fit_models="linear"
                                                        S_ridge = F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                                        overlap_function = "overlap2", # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                                        propensity_trim_threshold,
                                                        M=n_save)
  saveRDS(boot_linear_overlap2_l,
          paste0("/Results/Generalizability_results/boot_linear_overlap2_",l,"_nobs",n_obs_select,"_M",m,".rds"))
  
  
  toc()
  rm(boot_linear_overlap2_l);gc()
}

# Load data and combine bootstrap iterations
boot_linear_overlap2 = NULL
for(l in 1:n_loops){
  boot_linear_overlap2 = rbind(boot_linear_overlap2, readRDS(paste0("/Results/Generalizability_results/boot_linear_overlap2_",
                                                                    l,file_suffix,"_nobs",n_obs_select,"_M",m,".rds")))
}

saveRDS(boot_linear_overlap2,
        paste0("/Results/Generalizability_results/boot_linear_overlap2","_nobs",n_obs_select,"_M",m,".rds"))
rm(boot_linear_overlap2);gc()

toc()


# Stop "cluster"
stopCluster(cl)

toc()
print(Sys.time())
print("DONE")
