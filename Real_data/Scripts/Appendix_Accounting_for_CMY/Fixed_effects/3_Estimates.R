#################################################################################################################################################
### Name
### 8/18/21
### Fixed effects for country-month-year
#################################################################################################################################################
print(Sys.time())

### Load data & libraries
rm(list=ls())
path = getwd()

library(qs) # read data in qs format
library(Hmisc)
library(pbapply)
library(parallel)
library(data.table) # fast data manipulation
library(dplyr) # data manipulation
library(tictoc) # runtime
library(tidyverse) # for data manipulation
library(magrittr) # for data manipulation (set_colnames)  
library(lme4)

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
} else {
  target_sample = qread("/data/NYS_Medicaid/Name/target_sample_relevant_cmy.qs") # Restricted to those with 6 mo of observations (can measure outcome) -- will use for analysis

  # Separate rand and obs data
  rand_sample = target_sample %>% filter(aa_sample == 1)
}


# Calculate proportion of target sample made up of obs data: will need to weight obs component of estimate to account for using random sample of obs observations
# Target sample is a 10% random sample of the obs data; these numbers are from the full dataset; they were provided by Julia Yates
# n_rand = 67619 
# n_obs = 894884 
w_obs = 894884/(67619+894884)

# Relevant covariates
X_vars = c("age_clean120", #"age_scaled",
           "female", #"race", 
           "neighborhood",
           "aid_group_clean120", "ssi",
           "baseline_decile_clean", "baseline_decile_missing", "percent_poverty",
           "countyyearmonthfe") # age, sex, race, neighborhood, aid group, SSI eligibility, baseline spending decile, indicator of missing baseline spending, neighborhood povery rate, county-month-year

## Change fixed effects that exist in target_sample and don't exist in rand_sample to closests ones in rand_sample
new_countyyearmonthfe = unique(target_sample$countyyearmonthfe)[!(unique(target_sample$countyyearmonthfe) %in% unique(rand_sample$countyyearmonthfe))]
new_countyyearmonthfe # two new combos
sum(target_sample$countyyearmonthfe %in% new_countyyearmonthfe) # 312 observations out of 982427

target_sample$countyyearmonthfe = ifelse(target_sample$countyyearmonthfe==246, 245,
       ifelse(target_sample$countyyearmonthfe==247, 248, target_sample$countyyearmonthfe))

#################################################################################################################################################
## Point estimates, county-month-year fixed effects
#################################################################################################################################################
# Fit model
rand_fit = lm(lpay0tot_0to5 ~ ., data = rand_sample %>% 
                select("lpay0tot_0to5", "plan_baseline_clean", X_vars))

# Rand estimates
estimates_low = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "I", rand_sample[,X_vars]),
                        interval="confidence")
estimates_high = predict(rand_fit, cbind.data.frame("plan_baseline_clean" = "D", 
                                                    rand_sample[,X_vars]),
                         interval="confidence")

mean_low = apply(estimates_low, 2, function(x) mean(x, na.rm=T))
mean_high = apply(estimates_high, 2, function(x) mean(x, na.rm=T))

mean_low
mean_high
(mean_high - mean_low)/mean_low # 12%


# Target estimates
estimates_low_target = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "I", target_sample[,X_vars]),
                        interval="confidence")
estimates_high_target = predict(rand_fit, cbind.data.frame("plan_baseline_clean" = "D", 
                                                    target_sample[,X_vars]),
                         interval="confidence")

mean_low_target = apply(estimates_low_target, 2, function(x) mean(x, na.rm=T))
mean_high_target = apply(estimates_high_target, 2, function(x) mean(x, na.rm=T))

mean_low_target
mean_high_target
(mean_high_target - mean_low_target)/mean_low_target # 15%



print(Sys.time())
print("DONE")
