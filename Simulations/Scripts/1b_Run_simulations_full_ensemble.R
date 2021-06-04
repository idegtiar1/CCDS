####################################################################################################
### Title: 1_Run_simulations.R
### Author: 
### Description: Run paper 2 simulations across all data generating mechanisms
####################################################################################################

## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
Sys.time()

## ----load_libraries, output=F, message=FALSE------------------------------------------------------------------------------------------------------------------------
### Load libraries
rm(list=ls())  
path = getwd()     

library(tidyverse) # for data manipulation
library(nnet) # for nnet object output (multinomial regression fits)
library(glmnet) # for lasso regression output 
library(Hmisc) # for efficiently sampling from multinomial distribution
library(parallel) # for parallel processing
library(SuperLearner) # for ensembling

library(knitr) # pretty tables
library(kableExtra) # prettier tables
library(ggthemes) # pretty colors for ggplots 
library(gridExtra) # for arranging ggplots
library(magrittr) # for data manipulation (set_colnames) 
library(fst) # for storing large data frames/tables
library(qs) # for storing large data
library(ranger) # for random forest modeling
library(kernlab) # for ksvm modeling

source("/Helper_files/pw_overlap.R") # for determining regions of overlap and non-overlap, from Nethery et al. 2019: https://github.com/rachelnethery/overlap
source("/Helper_files/sim_functions.R") # for fitting models and simulating data   
source("/Helper_files/Lu2019_estimators.R") # for comparison estimators from Lu et al. 2019 paper

# Hyperparameters
m = 2000 # number of replications 
my_N=1000000 # target population size
my_n_target = 10000 # target sample size 
p_A_S1 = c(0.4,0.6) # P(A=1|S=1),P(A=2|S=1)
my_beta_AU = 0.625 # strength of relationship between U and A
my_beta_YU = c(10,0) # strength of relationship between U and Y: main effect and interaction with X1
my_beta_AX_linear = c(0.125,0.1,0.075,0.05) # coefficients for relationship between X and A
my_beta_AX = c(my_beta_AX_linear,0.1) # coefficients for relationship between X and A
my_beta_YX_linear = c(4,4,3,2) # coefficients for relationship between X and Y
my_beta_YX = c(my_beta_YX_linear,0.4) # coefficients for relationship between X and Y
my_beta_YAX_linear = c(4) # coefficient for interaction between X and A in Y regression: X1*A, (X1+1.0)^3*A
my_beta_YAX = c(my_beta_YAX_linear,0) # coefficient for interaction between X and A in Y regression: X1*A, 
random_seed=123 # random seed
no_cores = 100 # number of cores to use for simulations

# Set ggplot theme 
theme_set(theme_bw())

### Create library of algorithms
SL.glmnet.alt = create.Learner("SL.glmnet", tune = list(alpha = 0.5)) # Default alpha = 1, nfolds = 10, nlambda = 100
SL.ranger.alt = create.Learner("SL.ranger", tune = list(num.trees = 300, min.node.size=floor(my_n_target*0.05)), name_prefix = "SL.ranger_trees") # Default num.trees = 100, mtry = floor(sqrt(ncol(X))), min.node.size = ifelse(family$family == "gaussian", 5, 1)
SL.nnet.alt = create.Learner("SL.nnet", tune = list(size = 2)) 
my_SL_library = c("SL.glm", "SL.glm.interaction", SL.glmnet.alt$names,
               SL.ranger.alt$names, 
               SL.nnet.alt$names,
               "SL.earth", "SL.gam", "SL.kernelKnn")
file_suffix = "_full"


# Open "cluster"
cl=makeCluster(no_cores, outfile="",type="FORK")

# ## ----gen_data_complex-----------------------------------------------------------------------------------------------------------------------------------------------
# Read in each object from sim_complex
sim_data_complex = qread("/Results/Simulation_results/sim_data_complex.qs")
clusterExport(cl=cl,"sim_data_complex")

# Remove X1_cube term to fit linear models in simulation
sim_data_complex_noX1cube = subset(sim_data_complex,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_complex_noX1cube")

# Add X1_sq term to fit misspecified models in simulation
sim_data_complex_X1sq = sim_data_complex_noX1cube %>% mutate(X1_sq = X1^2)
clusterExport(cl=cl,"sim_data_complex_X1sq")

## ----run_sim_ensemble-----------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size 
                                  target_pop=sim_data_complex_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed 
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions: whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                                  SL_library=my_SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_ensemble,
        paste0("/Results/Simulation_results/sim_results_ensemble",file_suffix,".rds"))


## ----run_sim_ensemble_t_o-----------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_ensemble_t_o = parSapply(cl,1:m,
                                     function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size 
                                                                         target_pop=sim_data_complex_noX1cube, # target population from which to draw
                                                                         observed_U=F, # whether U is observed 
                                                                         complex_fit_models="ensemble", # algorithm to use to fit regressions # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                                                                         SL_library=my_SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                                                         overlap_function="ifelse(X_target[,1]>qnorm(0.9),0,
                                                                                         ifelse(X_target[,1]<qnorm(0.5),0,1))")}) # overlap region specification



# Save results
saveRDS(sim_results_ensemble_t_o,
        paste0("/Results/Simulation_results/sim_results_ensemble_t_o",file_suffix,".rds"))


# ----run_sim_n200_ensemble-----------------------------------------------------------------------------------------------------------------------------------------
# Change random forest hyperparameters corresponding to sample size
SL.ranger.alt = create.Learner("SL.ranger", tune = list(num.trees = 30, min.node.size=floor(200*0.05)), name_prefix = "SL.ranger_trees") # Default num.trees = 100, mtry = floor(sqrt(ncol(X))), min.node.size = ifelse(family$family == "gaussian", 5, 1)
clusterExport(cl=cl,"SL.ranger.alt")



clusterSetRNGStream(cl, iseed = random_seed)
sim_results_n200_ensemble = parSapply(cl,1:m,
                                      function(i,...){run_simulation_remove_errors(n_target =200, # target sample size
                                                                                   target_pop=sim_data_complex_noX1cube, # target population from which to draw
                                                                                   observed_U=F, # whether U is observed
                                                                                   complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                                   overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_n200_ensemble,
        paste0("/Results/Simulation_results/sim_results_n200_ensemble",file_suffix,".rds"))


# ----run_sim_n2000_ensemble-----------------------------------------------------------------------------------------------------------------------------------------
# Change random forest hyperparameters corresponding to sample size
SL.ranger.alt = create.Learner("SL.ranger", tune = list(num.trees = 300, min.node.size=floor(2000*0.05)), name_prefix = "SL.ranger_trees") # Default num.trees = 100, mtry = floor(sqrt(ncol(X))), min.node.size = ifelse(family$family == "gaussian", 5, 1)
clusterExport(cl=cl,"SL.ranger.alt")



clusterSetRNGStream(cl, iseed = random_seed)
sim_results_n2000_ensemble = parSapply(cl,1:m,
                                       function(i,...){run_simulation_remove_errors(n_target =2000, # target sample size
                                                                                    target_pop=sim_data_complex_noX1cube, # target population from which to draw
                                                                                    observed_U=F, # whether U is observed 
                                                                                    complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                                    overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_n2000_ensemble,
        paste0("/Results/Simulation_results/sim_results_n2000_ensemble",file_suffix,".rds"))


# ----run_sim_n50000_ensemble----------------------------------------------------------------------------------------------------------------------------------------
# Change random forest hyperparameters corresponding to sample size
SL.ranger.alt = create.Learner("SL.ranger", tune = list(num.trees = 300, min.node.size=floor(50000*0.05)), name_prefix = "SL.ranger_trees") # Default num.trees = 100, mtry = floor(sqrt(ncol(X))), min.node.size = ifelse(family$family == "gaussian", 5, 1)
clusterExport(cl=cl,"SL.ranger.alt")



clusterSetRNGStream(cl, iseed = random_seed)
sim_results_n50000_ensemble = parSapply(cl,1:m,
                                        function(i,...){run_simulation_remove_errors(n_target =50000, # target sample size
                                                                                     target_pop=sim_data_complex_noX1cube, # target population from which to draw
                                                                                     observed_U=F, # whether U is observed
                                                                                     complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                                     overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_n50000_ensemble,
        paste0("/Results/Simulation_results/sim_results_n50000_ensemble",file_suffix,".rds"))

# Change random forest hyperparameters corresponding to sample size
SL.ranger.alt = create.Learner("SL.ranger", tune = list(num.trees = 300, min.node.size=floor(my_n_target*0.05)), name_prefix = "SL.ranger_trees") # Default num.trees = 100, mtry = floor(sqrt(ncol(X))), min.node.size = ifelse(family$family == "gaussian", 5, 1)
clusterExport(cl=cl,"SL.ranger.alt")

# # # ## ----gen_data_complex2----------------------------------------------------------------------------------------------------------------------------------------------
# Read in each object from sim_complex2
sim_data_complex2 = qread("/Results/Simulation_results/sim_data_complex2.qs")
clusterExport(cl=cl,"sim_data_complex2")

# Remove X1_cube term to fit linear models in simulation
sim_data_complex2_noX1cube = subset(sim_data_complex2,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_complex2_noX1cube")

## ----run_sim_complex2_ensemble_t_o--------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_complex2_ensemble_t_o = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size 
                                  target_pop=sim_data_complex2_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                                  overlap_function="ifelse(X_target[,1]>qnorm(0.9),0,
                                  ifelse(X_target[,1]<qnorm(0.5),0,1))")}) # whether U is observed #



# Save results
saveRDS(sim_results_complex2_ensemble_t_o,
        paste0("/Results/Simulation_results/sim_results_complex2_ensemble_t_o",file_suffix,".rds"))



## ----run_sim_complex2_ensemble--------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_complex2_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size 
                                  target_pop=sim_data_complex2_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                                  overlap_function=F)}) # overlap region specification


#
# Save results
saveRDS(sim_results_complex2_ensemble,
        paste0("/Results/Simulation_results/sim_results_complex2_ensemble",file_suffix,".rds"))



## ----gen_data_linear------------------------------------------------------------------------------------------------------------------------------------------------
####################################################################################################
# # Read in each object from sim_linear
sim_data_linear = qread("/Results/Simulation_results/sim_data_linear.qs")
clusterExport(cl=cl,"sim_data_linear")

# ----run_sim_complex--------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_linear_ensemble = parSapply(cl,1:m,
                                        function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                                     target_pop=sim_data_linear, # target population from which to draw
                                                                                     observed_U=F, # whether U is observed
                                                                                     complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                                     overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_linear_ensemble,
        paste0("/Results/Simulation_results/sim_results_linear_ensemble",file_suffix,".rds"))


# ## ----gen_data_knot--------------------------------------------------------------------------------------------------------------------------------------------------
# #####################################################################################################
# # Read in each object from sim_knot
sim_data_knot = qread("/Results/Simulation_results/sim_data_knot.qs")
clusterExport(cl=cl,"sim_data_knot")

# Add knot term
sim_data_knot_withknot = sim_data_knot %>% mutate(X1_knot = (X1 < -1)*(X1 + 1))
clusterExport(cl=cl,"sim_data_knot_withknot")

sim_data_knot_noX1cube = subset(sim_data_knot,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_knot_noX1cube")

## ----run_sim_knot_ensemble------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_knot_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                  target_pop=sim_data_knot_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions
                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_knot_ensemble,
        paste0("/Results/Simulation_results/sim_results_knot_ensemble",file_suffix,".rds"))


## ----gen_data_knot2-------------------------------------------------------------------------------------------------------------------------------------------------
#####################################################################################################
# # Read in each object from sim_knot2
sim_data_knot2 = qread("/Results/Simulation_results/sim_data_knot2.qs")
clusterExport(cl=cl,"sim_data_knot2")

sim_data_knot2_noX1cube = subset(sim_data_knot2,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_knot2_noX1cube")

# Add knot term
sim_data_knot2_withknot = sim_data_knot2 %>% mutate(X1_knot = (X1 < -0.5)*(X1 + 0.5))
clusterExport(cl=cl,"sim_data_knot2_withknot")

# ----run_sim_knot2_ensemble-------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_knot2_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                  target_pop=sim_data_knot2_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions
                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_knot2_ensemble,
        paste0("/Results/Simulation_results/sim_results_knot2_ensemble",file_suffix,".rds"))

## ----gen_data_UX1---------------------------------------------------------------------------------------------------------------------------------------------------
#####################################################################################################
# Read in each object from sim_UX1
sim_data_UX1 = qread("/Results/Simulation_results/sim_data_UX1.qs")
clusterExport(cl=cl,"sim_data_UX1")

sim_data_UX1_noX1cube = subset(sim_data_UX1,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_UX1_noX1cube")

## ----run_sim_UX1_ensemble-------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_UX1_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                  target_pop=sim_data_UX1_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions
                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_UX1_ensemble,
        paste0("/Results/Simulation_results/sim_results_UX1_ensemble",file_suffix,".rds"))


## ----gen_data_less_conf---------------------------------------------------------------------------------------------------------------------------------------------
####################################################################################################
# # Read in each object from sim_less_conf
sim_data_less_conf = qread("/Results/Simulation_results/sim_data_less_conf.qs")
clusterExport(cl=cl,"sim_data_less_conf")

sim_data_less_conf_noX1cube = subset(sim_data_less_conf,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_less_conf_noX1cube")


## ----gen_data_more_conf---------------------------------------------------------------------------------------------------------------------------------------------
####################################################################################################
# # Read in each object from sim_more_conf
sim_data_more_conf = qread("/Results/Simulation_results/sim_data_more_conf.qs")
clusterExport(cl=cl,"sim_data_more_conf")

sim_data_more_conf_noX1cube = subset(sim_data_more_conf,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_more_conf_noX1cube")

# ## ----run_sim_U_measured_ensemble------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_U_measured_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                  target_pop=sim_data_complex_noX1cube, # target population from which to draw
                                  observed_U=T, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions
                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_U_measured_ensemble,
        paste0("/Results/Simulation_results/sim_results_U_measured_ensemble",file_suffix,".rds"))
# #
## ----run_sim_less_conf_ensemble--------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_less_conf_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                  target_pop=sim_data_less_conf_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions
                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_less_conf_ensemble,
        paste0("/Results/Simulation_results/sim_results_less_conf_ensemble",file_suffix,".rds"))
# #
## ----run_sim_more_conf_ensemble--------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_more_conf_ensemble = parSapply(cl,1:m,
               function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                  target_pop=sim_data_more_conf_noX1cube, # target population from which to draw
                                  observed_U=F, # whether U is observed
                                  complex_fit_models="ensemble", # algorithm to use to fit regressions
                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_more_conf_ensemble,
        paste0("/Results/Simulation_results/sim_results_more_conf_ensemble",file_suffix,".rds"))
 

## ----gen_data_less_overlap------------------------------------------------------------------------------------------------------------------------------------------
####################################################################################################
# # Read in each object from sim_less_overlap
sim_data_less_overlap = qread("/Results/Simulation_results/sim_data_less_overlap.qs")
clusterExport(cl=cl,"sim_data_less_overlap")

sim_data_less_overlap_noX1cube = subset(sim_data_less_overlap,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_less_overlap_noX1cube")

# ----run_sim_less_overlap_ensemble-----------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_less_overlap_ensemble = parSapply(cl,1:m,
                                         function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                             target_pop=sim_data_less_overlap_noX1cube, # target population from which to draw
                                                                             observed_U=F, # whether U is observed
                                                                             complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                             overlap_function=F)}) # overlap region specification


# Save results
saveRDS(sim_results_less_overlap_ensemble,
        paste0("/Results/Simulation_results/sim_results_less_overlap_ensemble",file_suffix,".rds"))


# ## ----gen_data_more_overlap------------------------------------------------------------------------------------------------------------------------------------------
# ####################################################################################################
# # Read in each object from sim_more_overlap
sim_data_more_overlap = qread("/Results/Simulation_results/sim_data_more_overlap.qs")
clusterExport(cl=cl,"sim_data_more_overlap")

sim_data_more_overlap_noX1cube = subset(sim_data_more_overlap,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_more_overlap_noX1cube")

# ----run_sim_more_overlap_ensemble-----------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_more_overlap_ensemble = parSapply(cl,1:m,
                                                      function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                                     target_pop=sim_data_more_overlap_noX1cube, # target population from which to draw
                                                                                     observed_U=F, # whether U is observed
                                                                                     complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                                     overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_more_overlap_ensemble,
        paste0("/Results/Simulation_results/sim_results_more_overlap_ensemble",file_suffix,".rds"))


# # ## ----gen_data_cbv---------------------------------------------------------------------------------------------------------------------------------------------------
# # #####################################################################################################
# # Read in each object from sim_cbv
sim_data_cbv = qread("/Results/Simulation_results/sim_data_cbv.qs")
clusterExport(cl=cl,"sim_data_cbv")

# Add knot term
sim_data_cbv_withknot = sim_data_cbv %>% mutate(X1_knot = (X1 < -0.5)*(X1 + 0.5))
clusterExport(cl=cl,"sim_data_cbv_withknot")

sim_data_cbv_noX1cube = subset(sim_data_cbv,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_cbv_noX1cube")


## ----run_sim_cbv_ensemble-------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_cbv_ensemble = parSapply(cl,1:m,
                                      function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                     target_pop=sim_data_cbv_noX1cube, # target population from which to draw
                                                                     observed_U=F, # whether U is observed
                                                                     complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                     overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_cbv_ensemble,
        paste0("/Results/Simulation_results/sim_results_cbv_ensemble",file_suffix,".rds"))


# # ## ----gen_data_U_f_X1------------------------------------------------------------------------------------------------------------------------------------------------
# # #####################################################################################################
# # Read in each object from sim_U_f_X1
sim_data_U_f_X1 = qread("/Results/Simulation_results/sim_data_U_f_X1.qs")
clusterExport(cl=cl,"sim_data_U_f_X1")

sim_data_U_f_X1_noX1cube = subset(sim_data_U_f_X1,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_U_f_X1_noX1cube")


## ----run_sim_U_f_X1_ensemble----------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_U_f_X1_ensemble = parSapply(cl,1:m,
                                     function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                    target_pop=sim_data_U_f_X1_noX1cube, # target population from which to draw
                                                                    observed_U=F, # whether U is observed
                                                                    complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                    overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_U_f_X1_ensemble,
        paste0("/Results/Simulation_results/sim_results_U_f_X1_ensemble",file_suffix,".rds"))



# ## ----gen_data_SU----------------------------------------------------------------------------------------------------------------------------------------------------
# #####################################################################################################
# # Read in each object from sim_SU
sim_data_SU = qread("/Results/Simulation_results/sim_data_SU.qs")
clusterExport(cl=cl,"sim_data_SU")

sim_data_SU_noX1cube = subset(sim_data_SU,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_SU_noX1cube")


## ----run_sim_SU_ensemble--------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_SU_ensemble = parSapply(cl,1:m,
                                   function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                  target_pop=sim_data_SU_noX1cube, # target population from which to draw
                                                                  observed_U=F, # whether U is observed
                                                                  complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                  overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_SU_ensemble,
        paste0("/Results/Simulation_results/sim_results_SU_ensemble",file_suffix,".rds"))


# # # ## ----gen_data_1to1--------------------------------------------------------------------------------------------------------------------------------------------------
# # # #####################################################################################################
# # Read in each object from sim_1to1
sim_data_1to1 = qread("/Results/Simulation_results/sim_data_1to1.qs")
clusterExport(cl=cl,"sim_data_1to1")

sim_data_1to1_noX1cube = subset(sim_data_1to1,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_1to1_noX1cube")

## ----run_sim_1to1_ensemble------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_1to1_ensemble =parSapply(cl,1:m,
                                     function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                    target_pop=sim_data_1to1_noX1cube, # target population from which to draw
                                                                    observed_U=F, # whether U is observed
                                                                    complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                    overlap_function=F)}) # overlap region specification


# Save results
saveRDS(sim_results_1to1_ensemble,
        paste0("/Results/Simulation_results/sim_results_1to1_ensemble",file_suffix,".rds"))


# # ## ----gen_data_1to4--------------------------------------------------------------------------------------------------------------------------------------------------
# # #####################################################################################################
# Read in each object from sim_1to4
sim_data_1to4 = qread("/Results/Simulation_results/sim_data_1to4.qs")
clusterExport(cl=cl,"sim_data_1to4")

sim_data_1to4_noX1cube = subset(sim_data_1to4,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_1to4_noX1cube")

## ----run_sim_1to4_ensemble------------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_1to4_ensemble = parSapply(cl,1:m,
                                       function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                      target_pop=sim_data_1to4_noX1cube, # target population from which to draw
                                                                      observed_U=F, # whether U is observed
                                                                      complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                      overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_1to4_ensemble,
        paste0("/Results/Simulation_results/sim_results_1to4_ensemble",file_suffix,".rds"))


# # ## ----gen_data_1to30-------------------------------------------------------------------------------------------------------------------------------------------------
# # #####################################################################################################
# Read in each object from sim_1to30
sim_data_1to30 = qread("/Results/Simulation_results/sim_data_1to30.qs")
clusterExport(cl=cl,"sim_data_1to30")

sim_data_1to30_noX1cube = subset(sim_data_1to30,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_1to30_noX1cube")

## ----run_sim_1to30_ensemble-----------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_1to30_ensemble = parSapply(cl,1:m,
                                        function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                       target_pop=sim_data_1to30_noX1cube, # target population from which to draw
                                                                       observed_U=F, # whether U is observed
                                                                       complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                       overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_1to30_ensemble,
        paste0("/Results/Simulation_results/sim_results_1to30_ensemble",file_suffix,".rds"))

## ----run_sim_overlap1_ensemble----------------------------------------------------------------------------------------------------------
sim_results_overlap1_ensemble = parSapply(cl,1:m,
function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                   target_pop = sim_data_complex_noX1cube, # target population from which to draw
                   observed_U = F, # whether U is observed
                   complex_fit_models = "ensemble", # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                   overlap_function = "overlap1")}) # # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3)



# Save results
saveRDS(sim_results_overlap1_ensemble,
        paste0("/Results/Simulation_results/sim_results_overlap1_ensemble",file_suffix,".rds"))



# # ## ----run_sim_overlap2_ensemble----------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_overlap2_ensemble = parSapply(cl,1:m,
function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                   target_pop = sim_data_complex_noX1cube, # target population from which to draw
                   observed_U = F, # whether U is observed
                   complex_fit_models = "ensemble", # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                   overlap_function = "overlap2")}) # # true overlap region or way to determine overlap (F=default, overlap2, overlap2, overlap3)



# Save results
saveRDS(sim_results_overlap2_ensemble,
        paste0("/Results/Simulation_results/sim_results_overlap2_ensemble",file_suffix,".rds"))

# # ## ----run_sim_overlap3_ensemble----------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_overlap3_ensemble = parSapply(cl,1:m,
function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                   target_pop = sim_data_complex_noX1cube, # target population from which to draw
                   observed_U = F, # whether U is observed
                   complex_fit_models = "ensemble", # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                   overlap_function = "overlap3")}) # # true overlap region or way to determine overlap (F=default, overlap3, overlap3, overlap3)



# Save results
saveRDS(sim_results_overlap3_ensemble,
        paste0("/Results/Simulation_results/sim_results_overlap3_ensemble",file_suffix,".rds"))



# # ## ----run_sim_overlap4_ensemble----------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_overlap4_ensemble = parSapply(cl,1:m,
function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                   target_pop = sim_data_complex_noX1cube, # target population from which to draw
                   observed_U = F, # whether U is observed
                   complex_fit_models = "ensemble", # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                   overlap_function = "overlap4")}) # # true overlap region or way to determine overlap (F=default, overlap4, overlap4, overlap4)



# Save results
saveRDS(sim_results_overlap4_ensemble,
        paste0("/Results/Simulation_results/sim_results_overlap4_ensemble",file_suffix,".rds"))



#####################################################################################################################
#####################################################################################################################
# # ## ----gen_data_s_probabilistic-------------------------------------------------------------------------------------------------------------------------------------------------
# # #####################################################################################################
# Read in each object from sim_s_probabilistic
sim_data_s_probabilistic = qread("/Results/Simulation_results/sim_data_s_probabilistic.qs")
clusterExport(cl=cl,"sim_data_s_probabilistic")

sim_data_s_probabilistic_noX1cube = subset(sim_data_s_probabilistic,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_s_probabilistic_noX1cube")

## ----run_sim_s_probabilistic_ensemble-----------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_s_probabilistic_ensemble = parSapply(cl,1:m,
                                                 function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                                              target_pop=sim_data_s_probabilistic_noX1cube, # target population from which to draw
                                                                                              observed_U=F, # whether U is observed
                                                                                              complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                                              overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_s_probabilistic_ensemble,
        paste0("/Results/Simulation_results/sim_results_s_probabilistic_ensemble",file_suffix,".rds"))


#####################################################################################################################
# # ## ----gen_data_s_probabilistic_c-------------------------------------------------------------------------------------------------------------------------------------------------
# # #####################################################################################################
# Read in each object from sim_s_probabilistic_c
sim_data_s_probabilistic_c = qread("/Results/Simulation_results/sim_data_s_probabilistic_c.qs")
clusterExport(cl=cl,"sim_data_s_probabilistic_c")

sim_data_s_probabilistic_c_noX1cube = subset(sim_data_s_probabilistic_c,select=c(-X1_cube))
clusterExport(cl=cl,"sim_data_s_probabilistic_c_noX1cube")

## ----run_sim_s_probabilistic_c_ensemble-----------------------------------------------------------------------------------------------------------------------------------------
clusterSetRNGStream(cl, iseed = random_seed)
sim_results_s_probabilistic_c_ensemble = parSapply(cl,1:m,
                                                   function(i,...){run_simulation_remove_errors(n_target = my_n_target, # target sample size
                                                                                                target_pop=sim_data_s_probabilistic_c_noX1cube, # target population from which to draw
                                                                                                observed_U=F, # whether U is observed
                                                                                                complex_fit_models="ensemble", # algorithm to use to fit regressions
                                                                                                overlap_function=F)}) # overlap region specification



# Save results
saveRDS(sim_results_s_probabilistic_c_ensemble,
        paste0("/Results/Simulation_results/sim_results_s_probabilistic_c_ensemble",file_suffix,".rds"))

## ----run_sim_boot_ensemble-----------------------------------------------------------------------------------------------------------------------------------------------
# clusterSetRNGStream(cl, iseed = random_seed)
# sim_results_boot_ensemble = parLapply(cl,1:m,
#                                      function(i,...){run_simulation_with_bootstrap(n_target = my_n_target, # target sample size 
#                                                         target_pop=sim_data_complex_noX1cube, # target population from which to draw
#                                                         observed_U=F, # whether U is observed
#                                                         complex_fit_models="ensemble", # algorithm to use to fit regressions # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
#                                                         SL_library=my_SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
#                                                         overlap_function=F, # overlap region specification
#                                                         M=500)}) # Number of bootstrap resamples
# 
# 
# 
# # Save results
# saveRDS(sim_results_boot_ensemble,
#         paste0("/Results/Simulation_results/sim_results_boot_ensemble",file_suffix,".rds"))

# Run m simulation iterations, each with M=1000 bootstrap replications, saving every iteration
for(l in 1:m){
  paste0("Bootstrap - ensemble ",l)
  clusterSetRNGStream(cl, iseed = (random_seed+l))
  sim_results_boot_ensemble_l = run_simulation_with_bootstrap(n_target = my_n_target, # target sample size 
                                                              target_pop=sim_data_complex_noX1cube, # target population from which to draw
                                                              observed_U=F, # whether U is observed
                                                              complex_fit_models="ensemble", # algorithm to use to fit regressions # whether to use ensemble or random forest (ranger) to fit models for Y, A, and S
                                                              SL_library=my_SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                                              overlap_function=F, # overlap region specification
                                                              M=1000, # Number of bootstrap resamples
                                                              parallelize = T) # whether to parallelize bootstrap
  saveRDS(sim_results_boot_ensemble_l,
          paste0("/Results/Simulation_results/boot_ensemble/sim_results_boot_ensemble_",l,file_suffix,".rds"))
  
  rm(sim_results_boot_ensemble_l);gc()
}

# Load data and combine bootstrap iterations
sim_results_boot_ensemble = NULL
for(l in 1:m){
  sim_results_boot_ensemble = rbind(boot_linear_overlap2, readRDS(paste0("/Results/Simulation_results/boot_ensemble/sim_results_boot_ensemble_",
                                                                         l,file_suffix,".rds")))
}

# Save results
saveRDS(sim_results_boot_ensemble,
        paste0("/Results/Simulation_results/sim_results_boot_ensemble",file_suffix,".rds"))


# Stop "cluster"
stopCluster(cl)

Sys.time()
print("DONE")