####################################################################################################
### Title: 0_Generate_data.R
### Author: 
### Description: Generate data from paper 2 simulations across all data generating mechanisms
###              See 2_summarize_simulation_results.Rmd for further details
####################################################################################################

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
library(mefa4) # for %notin% function 
#library(ranger) # for random forest modeling
library(kernlab) # for ksvm modeling

source('../Helper_files/pw_overlap.R') # for determining regions of overlap and non-overlap, from Nethery et al. 2019: https://github.com/rachelnethery/overlap
source('../Helper_files/sim_functions.R') # for fitting models and simulating data    

# Hyperparameters 
m = 1000 # number of replications 
my_N=1000000 # target population size
my_n_target = 10000 # target sample size 
my_p_A_S1 = c(0.4,0.6) # P(A=1|S=1),P(A=2|S=1)
my_beta_AU = 0.625 # strength of relationship between U and A 
my_beta_YU = c(10,0) # strength of relationship between U and Y: main effect and interaction with X1
my_beta_AX_linear = c(0.125,0.1,0.075,0.05) # coefficients for relationship between X and A
my_beta_AX = c(my_beta_AX_linear,0.1) # coefficients for relationship between X and A
my_beta_YX_linear = c(0.2,4,3,2) # coefficients for relationship between X and Y
my_beta_YX = c(my_beta_YX_linear,1) # coefficients for relationship between X and Y
my_beta_YAX_linear = c(2) # coefficient for interaction between X and A in Y regression: X1*A, (X1-0.5)^3*A
my_beta_YAX = c(my_beta_YAX_linear,0) # coefficient for interaction between X and A in Y regression: X1*A, 
random_seed=123 # random seed

####################################################################################################
# Base data generating mechanism (DGM)
sim_complex = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_complex = sim_complex$sim_data_conf
p_A_complex = sim_complex$p_A
p_S_adj_complex = sim_complex$p_S_adj
plot_p_A_complex = sim_complex$plot_p_A
plot_Y_complex = sim_complex$plot_Y
ave_bias_complex = sim_complex$ave_bias

# Save objects
qsave(sim_data_complex,
          "/Results & output/Simulations/sim_data_complex.qs")
qsave(p_A_complex,
          "/Results & output/Simulations/p_A_complex.qs")
qsave(p_S_adj_complex,
          "/Results & output/Simulations/p_S_adj_complex.qs")
ggsave(plot=plot_p_A_complex,filename="/Results & output/Simulations/sim_complex_plot_p_A.png") # *PAPER APPENDIX FIGURE 2A* (first panel of figure:propensity_scores)
ggsave(plot=plot_Y_complex,filename="/Results & output/Simulations/sim_complex_plot_Y.png")
qsave(ave_bias_complex,"/Results & output/Simulations/sim_complex_ave_bias.qs")


#####################################################################################################
# DGM with A*X1^3 interaction in outcome regression
sim_complex2 = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = c(2,2), # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_complex2 = sim_complex2$sim_data_conf
p_A_complex2 = sim_complex2$p_A
p_S_adj_complex2 = sim_complex2$p_S_adj
plot_p_A_complex2 = sim_complex2$plot_p_A
plot_Y_complex2 = sim_complex2$plot_Y
ave_bias_complex2 = sim_complex2$ave_bias

# Save objects
qsave(sim_data_complex2,
          "/Results & output/Simulations/sim_data_complex2.qs")
qsave(p_A_complex2,
          "/Results & output/Simulations/p_A_complex2.qs")
qsave(p_S_adj_complex2,
          "/Results & output/Simulations/p_S_adj_complex2.qs")
ggsave(plot=plot_p_A_complex2,filename="/Results & output/Simulations/sim_complex2_plot_p_A.png")
ggsave(plot=plot_Y_complex2,filename="/Results & output/Simulations/sim_complex2_plot_Y.png")
qsave(ave_bias_complex2,"/Results & output/Simulations/sim_complex2_ave_bias.qs")

#####################################################################################################
# DGM with knot
sim_knot = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "knot", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_knot = sim_knot$sim_data_conf
p_A_knot = sim_knot$p_A
p_S_adj_knot = sim_knot$p_S_adj
plot_p_A_knot = sim_knot$plot_p_A
plot_Y_knot = sim_knot$plot_Y
ave_bias_knot = sim_knot$ave_bias

# Save objects
qsave(sim_data_knot,
          "/Results & output/Simulations/sim_data_knot.qs")
qsave(p_A_knot,
          "/Results & output/Simulations/p_A_knot.qs")
qsave(p_S_adj_knot,
          "/Results & output/Simulations/p_S_adj_knot.qs")
ggsave(plot=plot_p_A_knot,filename="/Results & output/Simulations/sim_knot_plot_p_A.png")
ggsave(plot=plot_Y_knot,filename="/Results & output/Simulations/sim_knot_plot_Y.png")
qsave(ave_bias_knot,"/Results & output/Simulations/sim_knot_ave_bias.qs")


#####################################################################################################
# DGM with more severe knot
sim_knot2 = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "knot2", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_knot2 = sim_knot2$sim_data_conf
p_A_knot2 = sim_knot2$p_A
p_S_adj_knot2 = sim_knot2$p_S_adj
plot_p_A_knot2 = sim_knot2$plot_p_A
plot_Y_knot2 = sim_knot2$plot_Y
ave_bias_knot2 = sim_knot2$ave_bias

# Save objects
qsave(sim_data_knot2,
          "/Results & output/Simulations/sim_data_knot2.qs")
qsave(p_A_knot2,
          "/Results & output/Simulations/p_A_knot2.qs")
qsave(p_S_adj_knot2,
          "/Results & output/Simulations/p_S_adj_knot2.qs")
ggsave(plot=plot_p_A_knot2,filename="/Results & output/Simulations/sim_knot2_plot_p_A.png")
ggsave(plot=plot_Y_knot2,filename="/Results & output/Simulations/sim_knot2_plot_Y.png")
qsave(ave_bias_knot2,"/Results & output/Simulations/sim_knot2_ave_bias.qs")

#####################################################################################################
# DGM with knot inside overlap region
sim_knot3 = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "knot3", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_knot3 = sim_knot3$sim_data_conf
p_A_knot3 = sim_knot3$p_A
p_S_adj_knot3 = sim_knot3$p_S_adj
plot_p_A_knot3 = sim_knot3$plot_p_A
plot_Y_knot3 = sim_knot3$plot_Y
ave_bias_knot3 = sim_knot3$ave_bias

# Save objects
qsave(sim_data_knot3,
          "/Results & output/Simulations/sim_data_knot3.qs")
qsave(p_A_knot3,
          "/Results & output/Simulations/p_A_knot3.qs")
qsave(p_S_adj_knot3,
          "/Results & output/Simulations/p_S_adj_knot3.qs")
ggsave(plot=plot_p_A_knot3,filename="/Results & output/Simulations/sim_knot3_plot_p_A.png")
ggsave(plot=plot_Y_knot3,filename="/Results & output/Simulations/sim_knot3_plot_Y.png")
qsave(ave_bias_knot3,"/Results & output/Simulations/sim_knot3_ave_bias.qs")

#####################################################################################################
# DGM in which U is a function of X1
sim_UX1 = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "UX1", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_UX1 = sim_UX1$sim_data_conf
p_A_UX1 = sim_UX1$p_A
p_S_adj_UX1 = sim_UX1$p_S_adj
plot_p_A_UX1 = sim_UX1$plot_p_A
plot_Y_UX1 = sim_UX1$plot_Y
ave_bias_UX1 = sim_UX1$ave_bias

# Save objects
qsave(sim_data_UX1,
          "/Results & output/Simulations/sim_data_UX1.qs")
qsave(p_A_UX1,
          "/Results & output/Simulations/p_A_UX1.qs")
qsave(p_S_adj_UX1,
          "/Results & output/Simulations/p_S_adj_UX1.qs")
ggsave(plot=plot_p_A_UX1,filename="/Results & output/Simulations/sim_UX1_plot_p_A.png")
ggsave(plot=plot_Y_UX1,filename="/Results & output/Simulations/sim_UX1_plot_Y.png")
qsave(ave_bias_UX1,"/Results & output/Simulations/sim_UX1_ave_bias.qs")


####################################################################################################
# DGM in which there is less unmeasured confounding
sim_less_conf = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = 0.1, # strength of relationship between U and A
                     beta_YU = c(5,0), # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = F, # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_less_conf = sim_less_conf$sim_data_conf
p_A_less_conf = sim_less_conf$p_A
p_S_adj_less_conf = sim_less_conf$p_S_adj
plot_p_A_less_conf = sim_less_conf$plot_p_A
plot_Y_less_conf = sim_less_conf$plot_Y
ave_bias_less_conf = sim_less_conf$ave_bias

# Save objects
qsave(sim_data_less_conf,
          "/Results & output/Simulations/sim_data_less_conf.qs")
qsave(p_A_less_conf,
          "/Results & output/Simulations/p_A_less_conf.qs")
qsave(p_S_adj_less_conf,
          "/Results & output/Simulations/p_S_adj_less_conf.qs")
ggsave(plot=plot_p_A_less_conf,filename="/Results & output/Simulations/sim_less_conf_plot_p_A.png")
ggsave(plot=plot_Y_less_conf,filename="/Results & output/Simulations/sim_less_conf_plot_Y.png")
qsave(ave_bias_less_conf,"/Results & output/Simulations/sim_less_conf_ave_bias.qs")

####################################################################################################
# DGM in which there is more confounding
sim_more_conf = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = 1.5, # strength of relationship between U and A
                     beta_YU = c(20,0), # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = F, # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_more_conf = sim_more_conf$sim_data_conf
p_A_more_conf = sim_more_conf$p_A
p_S_adj_more_conf = sim_more_conf$p_S_adj
plot_p_A_more_conf = sim_more_conf$plot_p_A
plot_Y_more_conf = sim_more_conf$plot_Y
ave_bias_more_conf = sim_more_conf$ave_bias

# Save objects
qsave(sim_data_more_conf,
          "/Results & output/Simulations/sim_data_more_conf.qs")
qsave(p_A_more_conf,
          "/Results & output/Simulations/p_A_more_conf.qs")
qsave(p_S_adj_more_conf,
          "/Results & output/Simulations/p_S_adj_more_conf.qs")
ggsave(plot=plot_p_A_more_conf,filename="/Results & output/Simulations/sim_more_conf_plot_p_A.png")
ggsave(plot=plot_Y_more_conf,filename="/Results & output/Simulations/sim_more_conf_plot_Y.png")
qsave(ave_bias_more_conf,"/Results & output/Simulations/sim_more_conf_ave_bias.qs")


####################################################################################################
# DGM with linear outcome model
sim_linear = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y
                     beta_AX = my_beta_AX_linear, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX_linear, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX_linear, # coefficient for interaction between X and A in Y regression
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_linear = sim_linear$sim_data_conf
p_A_linear = sim_linear$p_A
p_S_adj_linear = sim_linear$p_S_adj
plot_p_A_linear = sim_linear$plot_p_A
plot_Y_linear = sim_linear$plot_Y
ave_bias_linear = sim_linear$ave_bias

# Save objects
qsave(sim_data_linear,
          "/Results & output/Simulations/sim_data_linear.qs")
qsave(p_A_linear,
          "/Results & output/Simulations/p_A_linear.qs")
qsave(p_S_adj_linear,
          "/Results & output/Simulations/p_S_adj_linear.qs")
ggsave(plot=plot_p_A_linear,filename="/Results & output/Simulations/sim_linear_plot_p_A.png")
ggsave(plot=plot_Y_linear,filename="/Results & output/Simulations/sim_linear_plot_Y.png")
qsave(ave_bias_linear,"/Results & output/Simulations/sim_linear_ave_bias.qs")


####################################################################################################
# DGM with less overlap (smaller overlap region)
sim_less_overlap = gen_conf_sim_data_simple(N=my_N, # target population size
                     overlap_factor=1, # controls amount of overlap
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = F, # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_less_overlap = sim_less_overlap$sim_data_conf
p_A_less_overlap = sim_less_overlap$p_A
p_S_adj_less_overlap = sim_less_overlap$p_S_adj
plot_p_A_less_overlap = sim_less_overlap$plot_p_A
plot_Y_less_overlap = sim_less_overlap$plot_Y
ave_bias_less_overlap = sim_less_overlap$ave_bias

# Save objects
qsave(sim_data_less_overlap,
          "/Results & output/Simulations/sim_data_less_overlap.qs")
qsave(p_A_less_overlap,
          "/Results & output/Simulations/p_A_less_overlap.qs")
qsave(p_S_adj_less_overlap,
          "/Results & output/Simulations/p_S_adj_less_overlap.qs")
ggsave(plot=plot_p_A_less_overlap,filename="/Results & output/Simulations/sim_less_overlap_plot_p_A.png")
ggsave(plot=plot_Y_less_overlap,filename="/Results & output/Simulations/sim_less_overlap_plot_Y.png")
qsave(ave_bias_less_overlap,"/Results & output/Simulations/sim_less_overlap_ave_bias.qs")


####################################################################################################
# DGM with more overlap (larger overlap region)
sim_more_overlap = gen_conf_sim_data_simple(N=my_N, # target population size
                     overlap_factor=4, # controls amount of overlap
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = F, # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_more_overlap = sim_more_overlap$sim_data_conf
p_A_more_overlap = sim_more_overlap$p_A
p_S_adj_more_overlap = sim_more_overlap$p_S_adj
plot_p_A_more_overlap = sim_more_overlap$plot_p_A
plot_Y_more_overlap = sim_more_overlap$plot_Y
ave_bias_more_overlap = sim_more_overlap$ave_bias

# Save objects
qsave(sim_data_more_overlap,
          "/Results & output/Simulations/sim_data_more_overlap.qs")
qsave(p_A_more_overlap,
          "/Results & output/Simulations/p_A_more_overlap.qs")
qsave(p_S_adj_more_overlap,
          "/Results & output/Simulations/p_S_adj_more_overlap.qs")
ggsave(plot=plot_p_A_more_overlap,filename="/Results & output/Simulations/sim_more_overlap_plot_p_A.png")
ggsave(plot=plot_Y_more_overlap,filename="/Results & output/Simulations/sim_more_overlap_plot_Y.png")
qsave(ave_bias_more_overlap,"/Results & output/Simulations/sim_more_overlap_ave_bias.qs")


#####################################################################################################
# DGM with constant bias violation
sim_cbv = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "cbv", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_cbv = sim_cbv$sim_data_conf
p_A_cbv = sim_cbv$p_A
p_S_adj_cbv = sim_cbv$p_S_adj
plot_p_A_cbv = sim_cbv$plot_p_A
plot_Y_cbv = sim_cbv$plot_Y
ave_bias_cbv = sim_cbv$ave_bias

# Save objects
qsave(sim_data_cbv,
          "/Results & output/Simulations/sim_data_cbv.qs")
qsave(p_A_cbv,
          "/Results & output/Simulations/p_A_cbv.qs")
qsave(p_S_adj_cbv,
          "/Results & output/Simulations/p_S_adj_cbv.qs")
ggsave(plot=plot_p_A_cbv,filename="/Results & output/Simulations/sim_cbv_plot_p_A.png")
ggsave(plot=plot_Y_cbv,filename="/Results & output/Simulations/sim_cbv_plot_Y.png")
qsave(ave_bias_cbv,"/Results & output/Simulations/sim_cbv_ave_bias.qs")

#####################################################################################################
# DGM with assumption violation: U is a function of X1
sim_U_f_X1 = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "U_f_X1", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_U_f_X1 = sim_U_f_X1$sim_data_conf
p_A_U_f_X1 = sim_U_f_X1$p_A
p_S_adj_U_f_X1 = sim_U_f_X1$p_S_adj
plot_p_A_U_f_X1 = sim_U_f_X1$plot_p_A
plot_Y_U_f_X1 = sim_U_f_X1$plot_Y
ave_bias_U_f_X1 = sim_U_f_X1$ave_bias

# Save objects
qsave(sim_data_U_f_X1,
          "/Results & output/Simulations/sim_data_U_f_X1.qs")
qsave(p_A_U_f_X1,
          "/Results & output/Simulations/p_A_U_f_X1.qs")
qsave(p_S_adj_U_f_X1,
          "/Results & output/Simulations/p_S_adj_U_f_X1.qs")
ggsave(plot=plot_p_A_U_f_X1,filename="/Results & output/Simulations/sim_U_f_X1_plot_p_A.png")
ggsave(plot=plot_Y_U_f_X1,filename="/Results & output/Simulations/sim_U_f_X1_plot_Y.png")
qsave(ave_bias_U_f_X1,"/Results & output/Simulations/sim_U_f_X1_ave_bias.qs")

#####################################################################################################
# DGM with exchangeability of sample selection assumption violation: S is a function of U
sim_SU = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "S_f_U", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_SU = sim_SU$sim_data_conf
p_A_SU = sim_SU$p_A
p_S_adj_SU = sim_SU$p_S_adj
plot_p_A_SU = sim_SU$plot_p_A
plot_Y_SU = sim_SU$plot_Y
ave_bias_SU = sim_SU$ave_bias

# Save objects
qsave(sim_data_SU,
          "/Results & output/Simulations/sim_data_SU.qs")
qsave(p_A_SU,
          "/Results & output/Simulations/p_A_SU.qs")
qsave(p_S_adj_SU,
          "/Results & output/Simulations/p_S_adj_SU.qs")
ggsave(plot=plot_p_A_SU,filename="/Results & output/Simulations/sim_SU_plot_p_A.png")
ggsave(plot=plot_Y_SU,filename="/Results & output/Simulations/sim_SU_plot_Y.png")
qsave(ave_bias_SU,"/Results & output/Simulations/sim_SU_ave_bias.qs")

#####################################################################################################
# DGM with 1:1 ratio of RCT to obs observations
sim_1to1 = gen_conf_sim_data_simple(N=my_N, # target population size
                     desired_p_S = 1/2, # P(S=1)
                     overlap_factor = 4, # controls amount of overlap
                     X1_upper_bound_quantile=2, # percent of the normal distribution in which we only observe rand data
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_1to1 = sim_1to1$sim_data_conf
p_A_1to1 = sim_1to1$p_A
p_S_adj_1to1 = sim_1to1$p_S_adj
plot_p_A_1to1 = sim_1to1$plot_p_A
plot_Y_1to1 = sim_1to1$plot_Y
ave_bias_1to1 = sim_1to1$ave_bias

# Save objects
qsave(sim_data_1to1,
          "/Results & output/Simulations/sim_data_1to1.qs")
qsave(p_A_1to1,
          "/Results & output/Simulations/p_A_1to1.qs")
qsave(p_S_adj_1to1,
          "/Results & output/Simulations/p_S_adj_1to1.qs")
ggsave(plot=plot_p_A_1to1,filename="/Results & output/Simulations/sim_1to1_plot_p_A.png")
ggsave(plot=plot_Y_1to1,filename="/Results & output/Simulations/sim_1to1_plot_Y.png")
qsave(ave_bias_1to1,"/Results & output/Simulations/sim_1to1_ave_bias.qs")


#####################################################################################################
# DGM with 1:4 ratio of RCT to obs observations
sim_1to4 = gen_conf_sim_data_simple(N=my_N, # target population size
                     desired_p_S = 1/5, # P(S=1)
                     overlap_factor = 4, # controls amount of overlap
                     X1_upper_bound_quantile=2, # percent of the normal distribution in which we only observe rand data
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_1to4 = sim_1to4$sim_data_conf
p_A_1to4 = sim_1to4$p_A
p_S_adj_1to4 = sim_1to4$p_S_adj
plot_p_A_1to4 = sim_1to4$plot_p_A
plot_Y_1to4 = sim_1to4$plot_Y
ave_bias_1to4 = sim_1to4$ave_bias

# Save objects
qsave(sim_data_1to4,
          "/Results & output/Simulations/sim_data_1to4.qs")
qsave(p_A_1to4,
          "/Results & output/Simulations/p_A_1to4.qs")
qsave(p_S_adj_1to4,
          "/Results & output/Simulations/p_S_adj_1to4.qs")
ggsave(plot=plot_p_A_1to4,filename="/Results & output/Simulations/sim_1to4_plot_p_A.png")
ggsave(plot=plot_Y_1to4,filename="/Results & output/Simulations/sim_1to4_plot_Y.png")
qsave(ave_bias_1to4,"/Results & output/Simulations/sim_1to4_ave_bias.qs")


#####################################################################################################
# DGM with 1:30 ratio of RCT to obs observations
sim_1to30 = gen_conf_sim_data_simple(N=10*my_N,  # target population size; increased to ensure sufficient individuals to draw from for rand data
                     desired_p_S = 1/31, # P(S=1)
                     overlap_factor = 4, # controls amount of overlap
                     X1_upper_bound_quantile=2, # percent of the normal distribution in which we only observe rand data
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_1to30 = sim_1to30$sim_data_conf
p_A_1to30 = sim_1to30$p_A
p_S_adj_1to30 = sim_1to30$p_S_adj
plot_p_A_1to30 = sim_1to30$plot_p_A
plot_Y_1to30 = sim_1to30$plot_Y
ave_bias_1to30 = sim_1to30$ave_bias

# Save objects
qsave(sim_data_1to30,
          "/Results & output/Simulations/sim_data_1to30.qs")
qsave(p_A_1to30,
          "/Results & output/Simulations/p_A_1to30.qs")
qsave(p_S_adj_1to30,
          "/Results & output/Simulations/p_S_adj_1to30.qs")
ggsave(plot=plot_p_A_1to30,filename="/Results & output/Simulations/sim_1to30_plot_p_A.png")
ggsave(plot=plot_Y_1to30,filename="/Results & output/Simulations/sim_1to30_plot_Y.png")
qsave(ave_bias_1to30,"/Results & output/Simulations/sim_1to30_ave_bias.qs")


#####################################################################################################
# DGM where S probabilistic not deterministic based on X1
sim_s_probabilistic = gen_conf_sim_data_simple(N=my_N, # target population size
                     p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                     beta_AU = my_beta_AU, # strength of relationship between U and A
                     beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                     beta_AX = my_beta_AX, # coefficients for relationship between X and A
                     beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                     beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                     complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                     data_gen_scenario = "s_probabilistic", # which alternative data generating process to apply
                     my_seed = random_seed) # random seed

# Extract objects
sim_data_s_probabilistic = sim_s_probabilistic$sim_data_conf
p_A_s_probabilistic = sim_s_probabilistic$p_A
p_S_adj_s_probabilistic = sim_s_probabilistic$p_S_adj
plot_p_A_s_probabilistic = sim_s_probabilistic$plot_p_A
plot_Y_s_probabilistic = sim_s_probabilistic$plot_Y
ave_bias_s_probabilistic = sim_s_probabilistic$ave_bias

# Save objects
qsave(sim_data_s_probabilistic,
          "/Results & output/Simulations/sim_data_s_probabilistic.qs")
qsave(p_A_s_probabilistic,
          "/Results & output/Simulations/p_A_s_probabilistic.qs")
qsave(p_S_adj_s_probabilistic,
          "/Results & output/Simulations/p_S_adj_s_probabilistic.qs")
ggsave(plot=plot_p_A_s_probabilistic,filename="/Results & output/Simulations/sim_s_probabilistic_plot_p_A.png")
ggsave(plot=plot_Y_s_probabilistic,filename="/Results & output/Simulations/sim_s_probabilistic_plot_Y.png")
qsave(ave_bias_s_probabilistic,"/Results & output/Simulations/sim_s_probabilistic_ave_bias.qs")


#####################################################################################################
# DGM where S probabilistic not deterministic based on X1 -- diff DGM
sim_s_probabilistic_b = gen_conf_sim_data_simple(N=my_N, # target population size
                                               p_A_S1 = my_p_A_S1, # P(A=1|S=1),P(A=2|S=1)
                                               beta_AU = my_beta_AU, # strength of relationship between U and A
                                               beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                                               beta_AX = my_beta_AX, # coefficients for relationship between X and A
                                               beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                                               beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                                               complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                                               data_gen_scenario = "s_probabilistic_b", # which alternative data generating process to apply
                                               my_seed = random_seed) # random seed

# Extract objects
sim_data_s_probabilistic_b = sim_s_probabilistic_b$sim_data_conf
p_A_s_probabilistic_b = sim_s_probabilistic_b$p_A
p_S_adj_s_probabilistic_b = sim_s_probabilistic_b$p_S_adj
plot_p_A_s_probabilistic_b = sim_s_probabilistic_b$plot_p_A
plot_Y_s_probabilistic_b = sim_s_probabilistic_b$plot_Y
ave_bias_s_probabilistic_b = sim_s_probabilistic_b$ave_bias

# Save objects
qsave(sim_data_s_probabilistic_b,
      "/Results & output/Simulations/sim_data_s_probabilistic_b.qs")
qsave(p_A_s_probabilistic_b,
      "/Results & output/Simulations/p_A_s_probabilistic_b.qs")
qsave(p_S_adj_s_probabilistic_b,
      "/Results & output/Simulations/p_S_adj_s_probabilistic_b.qs")
ggsave(plot=plot_p_A_s_probabilistic_b,filename="/Results & output/Simulations/sim_s_probabilistic_b_plot_p_A.png")
ggsave(plot=plot_Y_s_probabilistic_b,filename="/Results & output/Simulations/sim_s_probabilistic_b_plot_Y.png")
qsave(ave_bias_s_probabilistic_b,"/Results & output/Simulations/sim_s_probabilistic_b_ave_bias.qs")


#####################################################################################################
# DGM where S probabilistic not deterministic based on X1 -- diff DGM
sim_s_probabilistic_c = gen_conf_sim_data_simple(N=my_N, # target population size
                                                 p_A_S1 = c(0.1, 0.9), # P(A=1|S=1),P(A=2|S=1)
                                                 beta_AU = my_beta_AU, # strength of relationship between U and A
                                                 beta_YU = my_beta_YU, # strength of relationship between U and Y: main effect and interaction with X1
                                                 beta_AX = my_beta_AX, # coefficients for relationship between X and A
                                                 beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                                                 beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                                                 complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                                                 data_gen_scenario = "s_probabilistic_c", # which alternative data generating process to apply
                                                 my_seed = random_seed) # random seed

# Extract objects
sim_data_s_probabilistic_c = sim_s_probabilistic_c$sim_data_conf
p_A_s_probabilistic_c = sim_s_probabilistic_c$p_A
p_S_adj_s_probabilistic_c = sim_s_probabilistic_c$p_S_adj
plot_p_A_s_probabilistic_c = sim_s_probabilistic_c$plot_p_A
plot_Y_s_probabilistic_c = sim_s_probabilistic_c$plot_Y
ave_bias_s_probabilistic_c = sim_s_probabilistic_c$ave_bias

# Save objects
qsave(sim_data_s_probabilistic_c,
      "/Results & output/Simulations/sim_data_s_probabilistic_c.qs")
qsave(p_A_s_probabilistic_c,
      "/Results & output/Simulations/p_A_s_probabilistic_c.qs")
qsave(p_S_adj_s_probabilistic_c,
      "/Results & output/Simulations/p_S_adj_s_probabilistic_c.qs")
ggsave(plot=plot_p_A_s_probabilistic_c,filename="/Results & output/Simulations/sim_s_probabilistic_c_plot_p_A.png")
ggsave(plot=plot_Y_s_probabilistic_c,filename="/Results & output/Simulations/sim_s_probabilistic_c_plot_Y.png")
qsave(ave_bias_s_probabilistic_c,"/Results & output/Simulations/sim_s_probabilistic_c_ave_bias.qs")

