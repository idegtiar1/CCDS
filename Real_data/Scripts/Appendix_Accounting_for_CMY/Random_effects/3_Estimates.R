#################################################################################################################################################
### Name
### 8/6/21
### Analysis with random effects (longitudinal analysis)
#################################################################################################################################################
print(Sys.time())

### Load data & libraries
rm(list=ls())
# start.time1 <- proc.time()
path = getwd()
#path = "/home/degtiar/Research with Sherri - NYC_Medicaid/Code - Randomized data"

library(qs) # read data in qs format
#library(ck37r) # parallel TMLE
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
  target_sample_long = qread("/Data - simulated/target_sample_relevant_sim.qs")
  
  # Make smaller for testing
  target_sample_long = target_sample_long[1:5000,] %>% filter(plan_baseline_clean %in% c("A","B","C"))
  
  # Fix log spending not actually being logged
  target_sample_long$lpay0tot_0to5 = log(target_sample_long$lpay0tot_0to5 + 1)
  
  # Separate rand and obs data
  rand_sample_long = target_sample_long %>% filter(aa_sample == 2)
} else {
  target_sample_long = qread("/data/NYS_Medicaid/Name/target_sample_relevant_cmy_long.qs") # Restricted to those with 6 mo of observations (can measure outcome) -- will use for analysis

  # Separate rand and obs data
  rand_sample_long = target_sample_long %>% filter(aa_sample == 1)
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
           "baseline_decile_clean", "baseline_decile_missing", "percent_poverty") # age, sex, race, neighborhood, aid group, SSI eligibility, baseline spending decile, indicator of missing baseline spending, neighborhood povery rate

## Change fixed effects that exist in target_sample_long and don't exist in rand_sample_long to closest ones in rand_sample_long
new_countyyearmonthfe = unique(target_sample_long$countyyearmonthfe)[!(unique(target_sample_long$countyyearmonthfe) %in% unique(rand_sample_long$countyyearmonthfe))]
new_countyyearmonthfe # two new combos
sum(target_sample_long$countyyearmonthfe %in% new_countyyearmonthfe) # 312 observations out of 982427

target_sample_long$countyyearmonthfe = ifelse(target_sample_long$countyyearmonthfe==246, 245,
       ifelse(target_sample_long$countyyearmonthfe==247, 248, target_sample_long$countyyearmonthfe))

#################################################################################################################################################
## Point estimates, county-month-year random effects
#################################################################################################################################################
# Fit model
my_formula = as.formula(paste("lpay0tot", paste(paste(c("plan_baseline_clean", X_vars), collapse=" + "),
                                   " + (1|countyyearmonthfe)"), sep=" ~ "))
rand_fit = lmer(my_formula, data = rand_sample_long %>% 
                  select("lpay0tot", "plan_baseline_clean", X_vars, "countyyearmonthfe"))

# Rand estimates
estimates_low = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "I", 
                                         rand_sample_long[,c(X_vars, "countyyearmonthfe")]))
estimates_high = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "D", 
                                         rand_sample_long[,c(X_vars, "countyyearmonthfe")]))


mean_low = mean(estimates_low) 
mean_high = mean(estimates_high) 
mean_low
mean_high
(mean_high - mean_low)/mean_low # 32%

# Target estimates
estimates_low_target = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "I", 
                                         target_sample_long[,c(X_vars, "countyyearmonthfe")]))
estimates_high_target = predict(rand_fit, 
                         cbind.data.frame("plan_baseline_clean" = "D", 
                                          target_sample_long[,c(X_vars, "countyyearmonthfe")]))


mean_low_target = mean(estimates_low_target) 
mean_high_target = mean(estimates_high_target) 
mean_low_target # 1.34
mean_high_target # 1.97
(mean_high_target - mean_low_target)/mean_low_target # 47%

#################################################################################################################################################
## Point estimates, county-month-year and bene random effects
#################################################################################################################################################
# Fit model
my_formula = as.formula(paste("lpay0tot", paste(paste(c("plan_baseline_clean", X_vars), 
                                                      collapse=" + "),
                                                " + (1|countyyearmonthfe)  + (1|recip_id)"), 
                              sep=" ~ "))
rand_fit = lmer(my_formula, data = rand_sample_long %>% 
                  select("lpay0tot", "plan_baseline_clean", X_vars, "countyyearmonthfe", "recip_id"))

# Rand estimates
estimates_low = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "I", 
                                         rand_sample_long[,c(X_vars, "countyyearmonthfe", "recip_id")]))
estimates_high = predict(rand_fit, 
                         cbind.data.frame("plan_baseline_clean" = "D", 
                                          rand_sample_long[,c(X_vars, "countyyearmonthfe", "recip_id")]))


mean_low = mean(estimates_low) 
mean_high = mean(estimates_high) 
mean_low # 1.98
mean_high # 2.45
(mean_high - mean_low)/mean_low # 23%


#################################################################################################################################################
## Point estimates, county-month-year random effects in 6-mo aggregate spending analysis
#################################################################################################################################################
# Read in aggregated data
target_sample = qread("/data/NYS_Medicaid/Name/target_sample_relevant_cmy.qs") # Restricted to those with 6 mo of observations (can measure outcome) -- will use for analysis

# Separate rand and obs data
rand_sample = target_sample %>% filter(aa_sample == 1)

## Change fixed effects that exist in target_sample and don't exist in rand_sample to closests ones in rand_sample
new_countyyearmonthfe = unique(target_sample$countyyearmonthfe)[!(unique(target_sample$countyyearmonthfe) %in% unique(rand_sample$countyyearmonthfe))]
new_countyyearmonthfe # two new combos
sum(target_sample$countyyearmonthfe %in% new_countyyearmonthfe) # 52 benes out of 163822

target_sample$countyyearmonthfe = ifelse(target_sample$countyyearmonthfe==246, 245,
                                         ifelse(target_sample$countyyearmonthfe==247, 248, target_sample$countyyearmonthfe))

# Fit model
my_formula = as.formula(paste("lpay0tot_0to5", paste(paste(c("plan_baseline_clean", X_vars), collapse=" + "),
                                                " + (1|countyyearmonthfe)"), sep=" ~ "))
rand_fit = lmer(my_formula, data = rand_sample %>% 
                  select("lpay0tot_0to5", "plan_baseline_clean", X_vars, "countyyearmonthfe"))

# Rand estimates
estimates_low = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "I", 
                                         rand_sample[,c(X_vars, "countyyearmonthfe")]))
estimates_high = predict(rand_fit, 
                         cbind.data.frame("plan_baseline_clean" = "D", 
                                          rand_sample[,c(X_vars, "countyyearmonthfe")]))


mean_low = mean(estimates_low) 
mean_high = mean(estimates_high) 
mean_low # 4.2
mean_high # 4.7
(mean_high - mean_low)/mean_low # 12%

# Target estimates
estimates_low_target = predict(rand_fit, 
                        cbind.data.frame("plan_baseline_clean" = "I", 
                                         target_sample[,c(X_vars, "countyyearmonthfe")]))
estimates_high_target = predict(rand_fit, 
                         cbind.data.frame("plan_baseline_clean" = "D", 
                                          target_sample[,c(X_vars, "countyyearmonthfe")]))


mean_low_target = mean(estimates_low_target) 
mean_high_target = mean(estimates_high_target) 
mean_low_target # 3.4
mean_high_target # 3.9
(mean_high_target - mean_low_target)/mean_low_target # 15%



print(Sys.time())
print("DONE")
