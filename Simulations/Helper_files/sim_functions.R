#####################################################################################################
#####################################################################################################
## Function to generate data according to simple data generating mechanism involving 5 covariates
#####################################################################################################
#####################################################################################################
gen_conf_sim_data_simple = function(N=1000000, # target population size
                                    desired_p_S = 1/5, # P(S=1)
                                    overlap_factor = 2, # factor by which overlap region differs from 50/50 split - determines size of overlap region. Bigger factor = larger overlap region
                                    X1_upper_bound_quantile=10, # percent of the normal distribution in which we only observe RCT data
                                    p_A_S1 = c(0.4,0.6), # P(A=1|S=1),P(A=2|S=1)
                                    beta_AU = 0.6, # strength of relationship between U and A
                                    beta_YA = -3, # strength of relationship between A and Y
                                    beta_YU = c(2,1), # strength of relationship between U and Y: main effect and interaction with X1
                                    beta_AX = c(0.5,0.4,0.3,0.2), # coefficients for relationship between X and A
                                    beta_YX = c(10,8,6,4), # coefficients for relationship between Y and c(X1-X4,X1^3) (X1^3 included if complex_gen_model == T)
                                    beta_YAX = c(2), # coefficient for interaction between A and X1 in Y regression
                                    # Ya = -1.5 - beta_YA*sim_A + sim_X %*% beta_YX + beta_YU[1]*sim_U + beta_YU[2]*sim_U*sim_X[,1] + 
                                    #      beta_YAX[1]*sim_X[,1]*sim_A + beta_YAX[2]*X1_cube*sim_A + rnorm(N,0,1)
                                    complex_gen_model = F, # whether to include complex model terms (X1^3) in the function for S, A, and Y: F - don't include, T - include (X1-0.5)^3, Diff - include (X1+1)^3
                                    data_gen_scenario = F, # which alternative data generating process to apply
                                    my_seed=123){ # random seed
  # Define functions for ease
  logit = qlogis
  expit = plogis
  
  # Generate data
  if(complex_gen_model==F){ # without X1^3 term
    sim_X = cbind(rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1))
    colnames(sim_X) = c("X1","X2","X3","X4")
  } else if (complex_gen_model==T){ # with X1^3 term
    sim_X = cbind(rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1))
    X1_cube = (sim_X[,1]+1.0)^3
    sim_X = cbind(sim_X,X1_cube)
    colnames(sim_X) = c("X1","X2","X3","X4","X1_cube")
    
  } else {
    sim_X = cbind(rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1))
    X1_cube = (sim_X[,1]-0.5)^3
    sim_X = cbind(sim_X,X1_cube)
    colnames(sim_X) = c("X1","X2","X3","X4","X1_cube")
  }
  
  if(data_gen_scenario == "U_f_X1"){
    p_U = expit(30*sim_X[,1])
    sim_U = rbinom(length(p_U),size=1,p_U)
  } else {
    sim_U = rbinom(N,1,0.5)
  }
  
  
  # Create dataset with measured and unmeasured covariates
  sim_XU = cbind.data.frame(sim_X,sim_U)
  sim_XU_matrix = model.matrix(~.,sim_XU)[,-1]
  sim_X_matrix = sim_X
  
  #####################################################################################################
  ## Update S to ensure approximate n_rand:n_obs ratio desired
  #####################################################################################################
  print("Generating S to approximately ensure desired n_rand:n_obs ratio")
  X1_upper_bound = 1-X1_upper_bound_quantile/100 # everything above bound will only have support in the RCT data: half the RCT data will be above this threshold, half below
  X1_lower_bound = X1_upper_bound - overlap_factor/5 # everything below bound will only have support in the obs data
  
  set.seed(my_seed)
  if(data_gen_scenario == "S_f_U"){
    sim_S = ifelse(sim_X[,1]>qnorm(X1_upper_bound),1,
                   ifelse(sim_X[,1]<qnorm(X1_lower_bound),0,
                          rbinom(N,1,(0.25+sim_U*0.5)/overlap_factor))) # Assumes desired_p_S of 1/5 and X1_upper_bound_quantile=10!!!
  } else if (data_gen_scenario == "s_probabilistic"){
    # Obtain initial propensity scores
    xp   <- cbind(1, sim_X)
    beta_S <- c(-7,3,1,0.5,0.1,0.5) 
    p_S <- expit(xp %*% beta_S) # predicted propensity score for selection
    mean(p_S)
    
    # Change coefficients by different adjustment factors and see which one gets closest to the desired marginal probability
    adj_factors = t(as.matrix(seq(0,2,0.05))) # different beta adjustment factors that I'm trying
    colnames(adj_factors)=adj_factors
    p_S_adj_factors = apply(adj_factors,2,function(x) {
      beta_S_adj = c(beta_S[1]+x,beta_S[2:length(beta_S)]) # Change intercept
      expit(xp %*% beta_S_adj)}) # Try out different adjustment factors for adjusting P(S=1|X)
    mean_p_S_adj_factors = apply(p_S_adj_factors,2,mean)
    best_adj_factor = adj_factors[,which(abs(desired_p_S-mean_p_S_adj_factors)==min(abs(desired_p_S-mean_p_S_adj_factors)))] %>%
      as.numeric() # choose adjustment factor which gives you P(S=1) closest to desired_p_S
    best_adj_factor
    
    beta_S_adj = c(beta_S[1]+best_adj_factor,beta_S[2:length(beta_S)]) # Adjust beta_S to obtain P(S=1) approx = desired_p_S
    p_S_adj <- expit(xp %*% beta_S_adj)
    mean(p_S_adj)
    
    # Generate S from adjusted propensity score
    set.seed(my_seed)
    sim_S = rbinom(N,1,p_S_adj) # simulated values
    
  } else if (data_gen_scenario == "s_probabilistic_b"){
    # Obtain initial propensity scores
    xp   <- cbind(1, sim_X)
    beta_S <- c(-3.5,2,1,0.5,0.1,0.05) # s_probabilistic
    p_S <- expit(xp %*% beta_S) # predicted propensity score for selection
    mean(p_S)
    
    # Change coefficients by different adjustment factors and see which one gets closest to the desired marginal probability
    adj_factors = t(as.matrix(seq(0,2,0.05))) # different beta adjustment factors that I'm trying
    colnames(adj_factors)=adj_factors
    p_S_adj_factors = apply(adj_factors,2,function(x) {
      beta_S_adj = c(beta_S[1]+x,beta_S[2:length(beta_S)]) # Change intercept
      expit(xp %*% beta_S_adj)}) # Try out different adjustment factors for adjusting P(S=1|X)
    mean_p_S_adj_factors = apply(p_S_adj_factors,2,mean)
    best_adj_factor = adj_factors[,which(abs(desired_p_S-mean_p_S_adj_factors)==min(abs(desired_p_S-mean_p_S_adj_factors)))] %>%
      as.numeric() # choose adjustment factor which gives you P(S=1) closest to desired_p_S
    best_adj_factor
    
    beta_S_adj = c(beta_S[1]+best_adj_factor,beta_S[2:length(beta_S)]) # Adjust beta_S to obtain P(S=1) approx = desired_p_S
    p_S_adj <- expit(xp %*% beta_S_adj)
    mean(p_S_adj)
    
    # Generate S from adjusted propensity score
    set.seed(my_seed)
    sim_S = rbinom(N,1,p_S_adj) # simulated values
    
  } else if (data_gen_scenario == "s_probabilistic_c"){
    # Obtain initial propensity scores
    xp   <- cbind(1, sim_X, "sin_X1_cube_X2"=sin(sim_X[,"X1_cube"]*sim_X[,"X2"]))
    beta_S <- c(-7,2,1,0.5,0.1,0.1,0.1) # s_probabilistic
    p_S <- expit(xp %*% beta_S) # predicted propensity score for selection
    mean(p_S)
    
    # Change coefficients by different adjustment factors and see which one gets closest to the desired marginal probability
    adj_factors = t(as.matrix(seq(0,2,0.05))) # different beta adjustment factors that I'm trying
    colnames(adj_factors)=adj_factors
    p_S_adj_factors = apply(adj_factors,2,function(x) {
      beta_S_adj = c(beta_S[1]+x,beta_S[2:length(beta_S)]) # Change intercept
      expit(xp %*% beta_S_adj)}) # Try out different adjustment factors for adjusting P(S=1|X)
    mean_p_S_adj_factors = apply(p_S_adj_factors,2,mean)
    best_adj_factor = adj_factors[,which(abs(desired_p_S-mean_p_S_adj_factors)==min(abs(desired_p_S-mean_p_S_adj_factors)))] %>%
      as.numeric() # choose adjustment factor which gives you P(S=1) closest to desired_p_S
    best_adj_factor
    
    beta_S_adj = c(beta_S[1]+best_adj_factor,beta_S[2:length(beta_S)]) # Adjust beta_S to obtain P(S=1) approx = desired_p_S
    p_S_adj <- expit(xp %*% beta_S_adj)
    mean(p_S_adj)
    
    # Generate S from adjusted propensity score
    set.seed(my_seed)
    sim_S = rbinom(N,1,p_S_adj) # simulated values
    
  } else {
    sim_S = ifelse(sim_X[,1]>qnorm(X1_upper_bound),1,
                   ifelse(sim_X[,1]<qnorm(X1_lower_bound),0,
                          rbinom(N,1,5/overlap_factor*(desired_p_S-X1_upper_bound_quantile/100))))
  }
  
  
  # Get summary statistics
  sim_n_rand = sum(as.numeric(as.character(sim_S)))
  sim_n_obs = N - sim_n_rand
  
  
  #####################################################################################################
  ## Vary association of U with A and Y, and replace A,Y values with given dependence on U
  #####################################################################################################
  # Introduce unmeasured covariate U
  print("Introducing U and creating A-U dependence")
  
  ## New A that depends on unmeasured covariate U
  # Propensity score for treatment in randomized data = observed proportions on each treatment
  levels_A = c(0,1)
  set.seed(my_seed)
  sim_A_S1 = sample(levels_A,size=sim_n_rand,replace=T,prob=p_A_S1)
  sim_A_S1 = factor(sim_A_S1,levels=levels_A)
  
  # Estimate propensity score for treatment in observational data from model fit to observational data 
  if(data_gen_scenario == "s_probabilistic_c"){
    p_A_S0 = expit(-3.5+sim_X[sim_S==0,] %*% beta_AX + beta_AU*sim_U[sim_S==0])
  } else {
    p_A_S0 = expit(-0.8+sim_X[sim_S==0,] %*% beta_AX + beta_AU*sim_U[sim_S==0])
  }
  
  
  set.seed(my_seed)
  sim_A_S0 = rbinom(sim_n_obs,1,p_A_S0)
  sim_A_S0 = factor(sim_A_S0,levels=levels_A)
  
  # Check propensity scores
  # sample
  my_subset = sample(1:sim_n_obs, 1000, replace=TRUE)
  p_A_S0_sample = p_A_S0[my_subset,]
  
  plot_p_A = as.data.frame(p_A_S0_sample) %>% ggplot( aes(x=p_A_S0_sample)) +
    geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity') +
    theme_bw() + ggtitle("Propensity for A for S=0") + 
    xlab("P(A=1|X,U)")
  
  # Combine into one new variable A
  sim_p_A = cbind(rep(NA,N))
  sim_p_A[sim_S==1] = p_A_S1[1]
  sim_p_A[sim_S==0] = p_A_S0
  p_A_simple = cbind(sim_p_A,1-sim_p_A,0)
  
  sim_A = rep(NA,N)
  sim_A[sim_S==1] = as.numeric(as.character(sim_A_S1))
  sim_A[sim_S==0] = as.numeric(as.character(sim_A_S0))
  
  ####################################################################################################
  ## New Y that depends on unmeasured covariate U
  print("Generating Y and creating Y-U dependence")
  
  # Model matrices for prediction
  # Treatment
  sim_A1 = cbind("A"=rep(0,N))
  sim_A2 = cbind("A"=rep(1,N)) 
  
  set.seed(my_seed)
  if(complex_gen_model==F){
    #sim_Y  = -1.5 + beta_YA*sim_A  + sim_X %*% beta_YX + beta_YU*sim_U + 2*sim_X[,1]*sim_A  + rnorm(N,0,1)
    sim_Y1 = -1.5 + beta_YA*sim_A1 + sim_X %*% beta_YX + beta_YU[1]*sim_U + beta_YU[2]*sim_U*sim_X[,1] + 
      beta_YAX*sim_X[,1]*sim_A1 + rnorm(N,0,1)
    sim_Y2 = -1.5 + beta_YA*sim_A2 + sim_X %*% beta_YX + beta_YU[1]*sim_U + beta_YU[2]*sim_U*sim_X[,1] + 
      beta_YAX*sim_X[,1]*sim_A2 + rnorm(N,0,1)
  } else {
    #sim_Y  = -1.5 - beta_YA*sim_A  + sim_X %*% beta_YX + beta_YU*sim_U + 2*sim_X[,1]*sim_A  + 2*X1_cube + rnorm(N,0,1)
    sim_Y1 = -1.5 + beta_YA*sim_A1 + sim_X %*% beta_YX + beta_YU[1]*sim_U + beta_YU[2]*sim_U*sim_X[,1] + 
      beta_YAX[1]*sim_X[,1]*sim_A1 + beta_YAX[2]*X1_cube*sim_A1 + rnorm(N,0,1)
    sim_Y2 = -1.5 + beta_YA*sim_A2 + sim_X %*% beta_YX + beta_YU[1]*sim_U + beta_YU[2]*sim_U*sim_X[,1] + 
      beta_YAX[1]*sim_X[,1]*sim_A2 + beta_YAX[2]*X1_cube*sim_A2 + rnorm(N,0,1)
  }
  
  # Update if constant bias assumption doesn't hold
  if(data_gen_scenario=="knot"){
    sim_Y1 = sim_Y1 - 15*(sim_X[,1] < -1)*(sim_X[,1] + 1)
    sim_Y2 = sim_Y2 - 30*(sim_X[,1] < -1)*(sim_X[,1] + 1)
  }
  if(data_gen_scenario=="knot2"){
    sim_Y1 = sim_Y1 - 15*(sim_X[,1] < -0.5)*(sim_X[,1] + 0.5)
    sim_Y2 = sim_Y2 - 30*(sim_X[,1] < -0.5)*(sim_X[,1] + 0.5)
  }
  if(data_gen_scenario=="knot3"){
    sim_Y1 = sim_Y1 - 2*(sim_X[,1] < 0.5)*(sim_X[,1] - 0.5)
    sim_Y2 = sim_Y2 - 4*(sim_X[,1] < 0.5)*(sim_X[,1] - 0.5)
  }
  if(data_gen_scenario=="UX1"){
    sim_Y1 = sim_Y1 + 2*sim_U*sim_X[,1]
    sim_Y2 = sim_Y2 + 3*sim_U*sim_X[,1]
  }
  if(data_gen_scenario=="cbv"){
    sim_Y1 = sim_Y1 - 15*sim_U*(sim_X[,1] < -0.5)*(sim_X[,1] + 0.5)
    sim_Y2 = sim_Y2 - 30*sim_U*(sim_X[,1] < -0.5)*(sim_X[,1] + 0.5)
  }
  
  # Create simulated data
  sim_Y = ifelse(sim_A==0,sim_Y1,
                 ifelse(sim_A==1,sim_Y2, NA))
  sim_A = factor(sim_A,levels=c(0,1,2),labels=c("1","2","3"))
  sim_data_conf = cbind.data.frame(sim_Y, sim_S, sim_A, sim_X, sim_U,
                                   sim_Y1, sim_Y2)
  colnames(sim_data_conf) = c("Y","S","A",colnames(sim_X), "U",
                              "Y1","Y2")
  
  ## Plot: compare Y distributions pre-U and post-U
  # pre-U Y's
  set.seed(my_seed)
  # Sample
  i = sample(1:N, 1000, replace=TRUE)
  sim_Y_sample = sim_Y[i]
  sim_A_sample = sim_A[i]
  
  # Plot
  plot_Y = as.data.frame(cbind("A"=sim_A_sample,"Y"=sim_Y_sample)) %>% 
    mutate(A=as.factor(A),Y=as.numeric(as.character(Y))) %>% 
    ggplot( aes(x=Y, fill=A)) +
    geom_density(color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    theme_bw() + ggtitle("Y distribution") + 
    xlab("Y")
  
  # Check amount of external and internal validity bias
  ave_bias = rbind("Y1"=c("Pop truth"=mean(sim_Y1), # Pop truth
               "RCT truth"=mean(sim_Y1[sim_S==1]), # Truth in RCT
               "RCT observed"=mean(sim_Y[sim_A=="1" & sim_S==1]), # Observed in RCT
               "Obs truth"=mean(sim_Y1[sim_S==0]), # Truth in obs
               "Obs observed"=mean(sim_Y[sim_A=="1" & sim_S==0])), # Observed in obs
        "Y2"=c("Pop truth"=mean(sim_Y2), # Pop truth
               "RCT truth"=mean(sim_Y2[sim_S==1]), # Truth in RCT
               "RCT observed"=mean(sim_Y[sim_A=="2" & sim_S==1]), # Observed in RCT
               "Obs truth"=mean(sim_Y2[sim_S==0]), # Truth in obs
               "Obs observed"=mean(sim_Y[sim_A=="2" & sim_S==0]))) # Observed in obs
  
  #######################################################
  # Extract each object from sim_simple and save them
  sim_data_simple = cbind(sim_data_conf,"Y3"=0)
  
  # Get propensity scores
  fit_S = cbind.data.frame(sim_S,sim_X_matrix) %>% glm(formula(.), data=.,family="binomial")
  # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
  p_S_adj = predict(fit_S,type = "response")
  
  p_S_adj_simple = cbind(p_S_adj, 1-p_S_adj, 0)
  
  return(list(sim_data_conf=sim_data_simple,
              p_A=p_A_simple,
              p_S_adj=p_S_adj,
              plot_p_A=plot_p_A,
              plot_Y=plot_Y,
              ave_bias=ave_bias,
              overlap_upper_bound=X1_upper_bound,
              overlap_lower_bound=X1_lower_bound))
}
#####################################################################################################
#####################################################################################################
## Function to run simulation and compare estimators; returns target sample dataset for iterations with error
#####################################################################################################
# Note: the weighted2stageCCDS corresponds to the 2-stage CCDS presented in the manuscript. 
# The 2-stage HCCDS (hybrid CCDS) and weighted 2-stage HCCDS were not presented in the manuscript as their performance did not justify inclusion.

#####################################################################################################
## Error occurs for ksvm and with some algorithms in SuperLearner - Error in checkForRemoteErrors(val) : one node produced an error: 0 (non-NA) cases Calls: data.frame ... clusterApply -> staticClusterApply -> checkForRemoteErrors Execution halted
## Rerunning with dataset usually does not re-produce the error
run_simulation_remove_errors = function(n_target=my_n_target, # target sample size
                                 target_pop, # target population from which to draw
                                 observed_U = F, # whether unmeasured confounding exists
                                 complex_fit_models=F, # whether to use ensemble ("ensemble)", random forest ("ranger"), ksvm ("ksvm"), correctly specified ("correctly_specified"), or OLS regressions (F) for Y, A, and S
                                 SL_library=my_SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                 S_ridge=F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument (i.e., only used if complex_fit_models==F)
                                 overlap_function=F, # true overlap region (expression to evaluate) or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                 high_X1="qnorm(0.9)", # X1 threshold above which S=1 for all observations
                                 low_X1="qnorm(0.5)", # X1 threshold below which S=0 for all observations
                                 propensity_trim_threshold=0.001, # threshold at which to trim propensity scores
                                 Lu_estimators = T, # whether to include Lu et al. 2019 comparison estimators
                                 ...){ 
  
  ##########################################################################################
  ### Sample from population
  ##########################################################################################
  i = sample(1:nrow(target_pop), n_target, replace=TRUE)
  my_target_sample = target_pop[i,]
  
  out = tryCatch({
    get_estimates(my_target_sample, # target sample
                  n_target=n_target, # target sample size
                  observed_U = observed_U, # whether unmeasured confounding exists
                  complex_fit_models=complex_fit_models, # whether to use ensemble ("ensemble)", random forest ("ranger"), ksvm ("ksvm"), correctly specified ("correctly_specified"), or OLS regressions ("F") for Y, A, and S
                  SL_library=SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                  S_ridge=S_ridge, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                  overlap_function=overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                  high_X1=high_X1, # X1 threshold above which S=1 for all observations
                  low_X1=low_X1, # X1 threshold below which S=0 for all observations
                  propensity_trim_threshold=propensity_trim_threshold, # threshold at which to trim propensity scores
                  Lu_estimators=Lu_estimators, # whether to include Lu et al. 2019 comparison estimators
                  ...)
  },
  error=function(c) {
    return(my_target_sample) # target sample in which error occurred. Note: this can blow up output if a lot of errors occur
  })
return(out)
}

#####################################################################################################
#####################################################################################################
## Function to run simulation with bootstrap for confidence intervals
#####################################################################################################
#####################################################################################################
run_simulation_with_bootstrap = function(n_target=my_n_target, # target sample size
                                        target_pop, # target population from which to draw
                                        observed_U = F, # whether unmeasured confounding exists
                                        complex_fit_models=F, # whether to use ensemble ("ensemble)", random forest ("ranger"), ksvm ("ksvm"), correctly specified ("correctly_specified"), or OLS regressions ("F") for Y, A, and S
                                        SL_library=my_SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                        S_ridge=F, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                        overlap_function=F, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                        high_X1="qnorm(0.9)", # X1 threshold above which S=1 for all observations
                                        low_X1="qnorm(0.5)", # X1 threshold below which S=0 for all observations
                                        M=500, # Number of bootstrap resamples
                                        propensity_trim_threshold=0.001, # threshold at which to trim propensity scores
                                        parallelize = F, # whether to parallelize bootstrap
                                        Lu_estimators = F, # whether to include Lu et al. 2019 comparison estimators
                                        ...){
  
  ##########################################################################################
  ### Sample from population
  ##########################################################################################
  i = sample(1:nrow(target_pop), n_target, replace=TRUE)
  my_target_sample = target_pop[i,]
  
  ##########################################################################################
  ### Run bootstraps
  ##########################################################################################
  if(parallelize){
    boot_results = data.frame(t(parSapply(cl, 1:M,
                                          function(i,...){
                                            #print(i)
                                            
                                            # Bootstrap resample
                                            j = sample(1:nrow(my_target_sample), n_target, replace=TRUE)
                                            my_boot_sample = my_target_sample[j,]
                                            
                                            # Estimates
                                            out = tryCatch({
                                              get_estimates(my_boot_sample, # target sample
                                                            n_target=n_target, # target sample size
                                                            observed_U = observed_U, # whether unmeasured confounding exists
                                                            complex_fit_models=complex_fit_models, # whether to use ensemble ("ensemble)", random forest ("ranger"), ksvm ("ksvm"), correctly specified ("correctly_specified"), or OLS regressions ("F") for Y, A, and S
                                                            SL_library=SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                                            S_ridge=S_ridge, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                                            overlap_function=overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                                            high_X1=high_X1, # X1 threshold above which S=1 for all observations
                                                            low_X1=low_X1, # X1 threshold below which S=0 for all observations
                                                            propensity_trim_threshold=propensity_trim_threshold, # threshold at which to trim propensity scores
                                                            Lu_estimators=Lu_estimators,...) # whether to include Lu et al. 2019 comparison estimators
                                            },
                                            error=function(c) {
                                              if(complex_fit_models=="ensemble"){
                                                if(Lu_estimators){
                                                  rep(NA,79+11*length(SL_library)+9*4) # 79 = number of output columns for non-ensemble output + algorithm weight columns + Lu output columns
                                                } else {
                                                  rep(NA,79+11*length(SL_library)) # 79 = number of output columns for non-ensemble output + algorithm weight columns 
                                                }
                                                                                              
                                              } else {
                                                if(Lu_estimators){
                                                  rep(NA,79+9*4)
                                                } else {
                                                  rep(NA,79)
                                                }
                                              }
                                            })
                                          })))
  } else {
    boot_results = data.frame(t(sapply(1:M,
                                          function(i,...){
                                            #print(i)
                                            
                                            # Bootstrap resample
                                            j = sample(1:nrow(my_target_sample), n_target, replace=TRUE)
                                            my_boot_sample = my_target_sample[j,]
                                            
                                            # Estimates
                                            out = tryCatch({
                                              get_estimates(my_boot_sample, # target sample
                                                            n_target=n_target, # target sample size
                                                            observed_U = observed_U, # whether unmeasured confounding exists
                                                            complex_fit_models=complex_fit_models, # whether to use ensemble ("ensemble)", random forest ("ranger"), ksvm ("ksvm"), correctly specified ("correctly_specified"), or OLS regressions ("F") for Y, A, and S
                                                            SL_library=SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                                            S_ridge=S_ridge, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                                                            overlap_function=overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                                                            high_X1=high_X1, # X1 threshold above which S=1 for all observations
                                                            low_X1=low_X1, # X1 threshold below which S=0 for all observations
                                                            propensity_trim_threshold=propensity_trim_threshold, # threshold at which to trim propensity scores
                                                            Lu_estimators=Lu_estimators,...) # whether to include Lu et al. 2019 comparison estimators
                                            },
                                            error=function(c) {
                                              if(complex_fit_models=="ensemble"){
                                                if(Lu_estimators){
                                                  rep(NA,79+11*length(SL_library)+9*4) # 70 = number of output columns for non-ensemble output + algorithm weight columns + Lu output columns
                                                } else {
                                                  rep(NA,79+11*length(SL_library)) # 70 = number of output columns for non-ensemble output + algorithm weight columns
                                                }
                                                
                                              } else {
                                                if(Lu_estimators){
                                                  rep(NA,79+9*4)
                                                } else {
                                                  rep(NA,79)
                                                }
                                              }
                                            })
                                          })))
  }
  
  

  ##########################################################################################
  ### Extract Y1 and Y2 results, generate Y1 - Y2 results
  ##########################################################################################
  mean_Y1 = boot_results[grepl("mean_Y1_", names(boot_results))]
  mean_Y2 = boot_results[grepl("mean_Y2_", names(boot_results))]
  mean_Y1minusY2 = mean_Y1 - mean_Y2
  colnames(mean_Y1minusY2) = gsub("Y1", "Y1minusY2", colnames(mean_Y1minusY2))
  
  results_means = cbind.data.frame(mean_Y1, mean_Y2, mean_Y1minusY2)[,order(c(1:ncol(mean_Y1), 
                                                                              1:ncol(mean_Y2),
                                                                              1:ncol(mean_Y1minusY2)))]
  
  n_estimators = length(mean_Y1)
  
  ##########################################################################################
  ### Obtain confidence intervals, whether those bounds cover the truth, and SE
  ##########################################################################################
  truth = apply(target_pop[,c("Y1","Y2")],2,mean)
  truth["Y1minusY2"] = truth["Y1"] - truth["Y2"]
  
  bounds = apply(results_means,2, 
                 function(x) quantile(x,
                                      probs = c(0.025,0.975),
                                      na.rm = TRUE))
  coverage = NULL
  for(c in 1:n_estimators){coverage = c(coverage, 
                                        bounds[1,(3*c-2):(3*c)] < truth & 
                                          truth < bounds[2,(3*c-2):(3*c)])}
  
  boot_se = apply(results_means,2, 
                  function(x) sqrt(var(x, na.rm = TRUE)))
  
  ci_width = bounds[2,]-bounds[1,]
  
  out = list("boot_se" = boot_se, 
             "ci_width" = ci_width, 
             "coverage" = coverage, 
             "bounds" = bounds)
  return(out)
}


#####################################################################################################
#####################################################################################################
## Function to get estimates from CCDS and comparison estimators
#####################################################################################################
#####################################################################################################
get_estimates = function(target_sample = my_target_sample, # target sample
                         n_target, # target sample size
                         observed_U, # whether unmeasured confounding exists
                         complex_fit_models, # whether to use ensemble ("ensemble)", random forest ("ranger"), ksvm ("ksvm"), correctly specified ("correctly_specified"), or OLS regressions ("F") for Y, A, and S
                         SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                         S_ridge, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                         overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                         high_X1, # X1 threshold above which S=1 for all observations
                         low_X1, # X1 threshold below which S=0 for all observations
                         propensity_trim_threshold, # threshold at which to trim propensity scores
                         Lu_estimators, # whether to include Lu et al. 2019 comparison estimators
                         ...){
  
  # Extract data
  if(observed_U==T){ # whether U used as "measured" covariate
    X_target = subset(target_sample,select=c(-Y,-S,-A,-Y1,-Y2,-Y3))
    X1to4_target = subset(target_sample,select=c(X1, X2, X3, X4, U))
  } else {
    X_target = subset(target_sample,select=c(-Y,-S,-A,-U,-Y1,-Y2,-Y3))
    X1to4_target = subset(target_sample,select=c(X1, X2, X3, X4))
  }
  X_target_matrix = model.matrix(~.,X_target)[,-1]
  X1to4_target_matrix = model.matrix(~.,X1to4_target)[,-1]
  Y_target = target_sample %>% pull(Y)
  S_target = target_sample %>% pull(S)
  A_target = target_sample %>% pull(A)
  levels_A = unique(A_target) %>% sort()
  A_target_matrix = model.matrix(~ ., data.frame(A_target))[,-1] %>% set_colnames(.,c("A2","A3"))
  if(sum(A_target_matrix[,2])==0){A_target_matrix=A_target_matrix[,1]}
  X1_high_target = ifelse(X_target_matrix[,1] > eval(parse(text=high_X1)),1,0)
  X1_low_target = ifelse(X_target_matrix[,1] < eval(parse(text=low_X1)),1,0)
  
  # Counterfactual outcomes
  Y1_target = target_sample %>% pull(Y1)
  Y2_target = target_sample %>% pull(Y2)
  Y3_target = rep(NA,n_target) 
  Y2minusY1_target = Y2_target - Y1_target
  
  # Subsetting randomized and observational data
  Y_rand = Y_target[S_target==1]
  Y1_rand = Y1_target[S_target==1]
  Y2_rand = Y2_target[S_target==1]
  Y3_rand = Y3_target[S_target==1]
  A_rand = A_target[S_target==1] # 1, 2, 3
  A_rand_matrix = model.matrix(~ ., data.frame(A_rand))[,-1] %>% set_colnames(.,c("A2","A3")) # 0, 1
  if(sum(A_rand_matrix[,2])==0){A_rand_matrix=A_rand_matrix[,1]}
  X_rand = X_target[S_target==1,]
  X_rand_matrix = model.matrix(~.,X_rand)[,-1]
  X1to4_rand = X1to4_target[S_target==1,]
  X1to4_rand_matrix = model.matrix(~.,X1to4_rand)[,-1]
  X1_high_rand = X1_high_target[S_target==1]
  X1_low_rand = X1_low_target[S_target==1]
  
  Y_obs = Y_target[S_target==0]
  Y1_obs = Y1_target[S_target==0]
  Y2_obs = Y2_target[S_target==0]
  Y3_obs = Y3_target[S_target==0]
  A_obs = A_target[S_target==0]
  A_obs_matrix = model.matrix(~ ., data.frame(A_obs))[,-1] %>% set_colnames(.,c("A2","A3"))
  if(sum(A_obs_matrix[,2])==0){A_obs_matrix=A_obs_matrix[,1]}
  X_obs = X_target[S_target==0,]
  X_obs_matrix = model.matrix(~.,X_obs)[,-1]
  X1to4_obs = X1to4_target[S_target==0,]
  X1to4_obs_matrix = model.matrix(~.,X1to4_obs)[,-1]
  X1_high_obs = X1_high_target[S_target==0]
  X1_low_obs = X1_low_target[S_target==0]
  
  # Target sample: number of randomized individuals and number of observational individuals
  n_rand = sum(as.numeric(as.character(S_target)))
  n_obs = n_target - n_rand
  
  ##########################################################################################
  ### Determine overlap and non-overlap regions
  ### Use method by Nethery et al. 2019
  ##########################################################################################
  # Propensity score for selection (U observed)
  if(complex_fit_models=="ensemble"){
    fit_S = SuperLearner(Y=S_target, X=as.data.frame(X_target_matrix), family="binomial", SL.library=SL_library, verbose=F) 
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = predict(fit_S, as.data.frame(X_target_matrix), onlySL = TRUE)$pred 
    
  } else if(complex_fit_models=="ksvm"){
    fit_S = cbind.data.frame(S_target,X_target_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = predict(fit_S, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    
  } else if(complex_fit_models=="ranger"){
    # Ranger: random forest
    fit_S = cbind.data.frame(S_target,X_target_matrix) %>% ranger(formula(.), data=.,
                                                                  num.trees = 300, min.node.size=floor(nrow(X_target_matrix)*0.05), # these hyperparameters work reasonably in tests
                                                                  num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = cbind.data.frame(S_target,X_target_matrix) %>%
      predict(fit_S, .,type = "response") %>% .$predictions %>% .[,2]
    
  } else if(complex_fit_models=="correctly_specified"){
    # When just including indicators, overlap region not found for some reason 
    fit_S = cbind.data.frame(S_target,X_target_matrix,X1_high_target,X1_low_target) %>% 
      glm(formula(.), data=.,family="binomial")
    pi_S = predict(fit_S,type = "response")
    
  } else if(S_ridge==T){
    SX_target <- cbind.data.frame(S_target,X_target_matrix)
    fit_S <- cv.glmnet(X_target_matrix,S_target, data = SX_target, alpha=0,family="binomial")
    
    pi_S <- predict(fit_S,X_target_matrix, type="response") # predicted probability that y=1
    
  } else {
    fit_S = cbind.data.frame(S_target,X_target_matrix) %>% glm(formula(.), data=.,family="binomial")
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = predict(fit_S,type = "response")
    
  }
  
  # Trim pi_S for logit
  pi_S_trim = ifelse(pi_S<propensity_trim_threshold,propensity_trim_threshold,
                     ifelse(pi_S>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_S))
  
  logit = qlogis
  logit_pi_S = logit(pi_S_trim)
  
  if(overlap_function == "overlap1"){ # uses nominal propensity score
    a = (max(pi_S)-min(pi_S)) %>% abs() %>% as.numeric()*.01 # a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 1% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=pi_S, # propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function == "overlap2"){ # uses nominal propensity score
    a = (max(pi_S)-min(pi_S)) %>% abs() %>% as.numeric()*.02 # a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 2% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=pi_S, # propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function == "overlap3"){ # uses logit propensity score
    a = (max(logit_pi_S)-min(logit_pi_S)) %>% abs() %>% as.numeric()*.1# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.04*min(n_obs,n_rand),0) # less than 4% of randomized units in 10% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=logit_pi_S, # logit propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function == "overlap4"){ # uses logit propensity score
    a = (max(logit_pi_S)-min(logit_pi_S)) %>% abs() %>% as.numeric()*.02# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 2% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=logit_pi_S, # logit propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function==F){ # uses logit propensity score
    # Tuning parameters for determining region of overlap
    # Given the skewdness of the data, determining overlap on logit scale
    a = (max(logit_pi_S)-min(logit_pi_S)) %>% abs() %>% as.numeric()*.01# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = max(round(0.01*min(n_obs,n_rand),0),1) # less than 1% of randomized units in 1% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap -- slow
    overlap_target = pw_overlap(ps=logit_pi_S, # propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else { # user-specified overlap region
    overlap_target=eval(parse(text=overlap_function))
  }
  
  # Regions of overlap
  overlap_rand = overlap_target[S_target==1]
  overlap_obs = overlap_target[S_target==0]
  X_rand_overlap = X_rand[overlap_rand==1,]
  X_obs_overlap = X_obs[overlap_obs==1,]
  X_rand_overlap_matrix = model.matrix(~.,X_rand_overlap)[,-1]
  X_obs_overlap_matrix = model.matrix(~.,X_obs_overlap)[,-1]
  
  X1to4_rand_overlap = X1to4_rand[overlap_rand==1,]
  X1to4_obs_overlap = X1to4_obs[overlap_obs==1,]
  X1to4_rand_overlap_matrix = model.matrix(~.,X1to4_rand_overlap)[,-1]
  X1to4_obs_overlap_matrix = model.matrix(~.,X1to4_obs_overlap)[,-1]
  
  mean_overlap_target = mean(overlap_target) # proportion of data in region of overlap
  mean_overlap_rand = mean(overlap_rand) # proportion of randomized data in region of overlap
  mean_overlap_obs = mean(overlap_obs) # proportion of observational data in region of overlap
  
  # Bounds of overlap region
  order_pi_S = pi_S[order(pi_S)]
  order_overlap_target = overlap_target[order(pi_S)]
  order_S_target = S_target[order(S_target)]
  min_p_overlap = order_pi_S[min(which(order_overlap_target==1))] # smallest propensity score in the region of overlap
  max_p_overlap = order_pi_S[max(which(order_overlap_target==1))] # largest propensity score in the region of overlap
  
  # Remove unnecessary data
  rm(X_rand_overlap, X_obs_overlap, X1to4_rand, X1to4_rand_overlap, X1to4_obs_overlap)
  gc()
  
  ##########################################################################################
  ### Determine propensities and weights used in estimators 
  ##########################################################################################
  if(complex_fit_models=="ensemble"){
    # P(R_overlap|S=1,X) 
    # Warning message due to perfect fit: Warning messages: 1: glm.fit: algorithm did not converge 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
    fit_Roverlap_rand = SuperLearner(Y=overlap_rand, X=as.data.frame(X_rand_matrix), 
                                     family="binomial", SL.library=SL_library, verbose=F)
    pi_Roverlap_rand = predict(fit_Roverlap_rand, as.data.frame(X_target_matrix), onlySL = TRUE)$pred 
    
    # P(R_overlap|S=0,X)
    # Warning message due to perfect fit: Warning messages: 1: glm.fit: algorithm did not converge 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
    fit_Roverlap_obs = SuperLearner(Y=overlap_obs, X=as.data.frame(X_obs_matrix), 
                                    family="binomial", SL.library=SL_library, verbose=F)
    pi_Roverlap_obs = predict(fit_Roverlap_obs, as.data.frame(X_target_matrix), onlySL = TRUE)$pred 
    
    # P(A|S=1,X)
    fit_A_rand = SuperLearner(Y=A_rand_matrix, X=as.data.frame(X_rand_matrix), # Note: only works for binary A; need to create wrappers for multinomial classification: https://github.com/ecpolley/SuperLearner/issues/16
                              family="binomial", SL.library=SL_library, verbose=F)
    pi_A_rand = predict(fit_A_rand, as.data.frame(X_target_matrix), onlySL = TRUE)$pred 
    
    # P(A|S=0,X)
    fit_A_obs = SuperLearner(Y=A_obs_matrix, X=as.data.frame(X_obs_matrix), 
                             family="binomial", SL.library=SL_library, verbose=F)
    pi_A_obs = predict(fit_A_obs, as.data.frame(X_target_matrix), onlySL = TRUE)$pred 
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap = SuperLearner(Y=A_rand_matrix[overlap_rand==1], X=as.data.frame(X_rand_matrix)[overlap_rand==1,], 
                                      family="binomial", SL.library=SL_library, verbose=F)
    pi_A_rand_overlap = predict(fit_A_rand_overlap, as.data.frame(X_target_matrix), onlySL = TRUE)$pred
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap = SuperLearner(Y=A_obs_matrix[overlap_obs==1], X=as.data.frame(X_obs_matrix)[overlap_obs==1,], 
                                     family="binomial", SL.library=SL_library, verbose=F)
    pi_A_obs_overlap = predict(fit_A_obs_overlap, as.data.frame(X_target_matrix), onlySL = TRUE)$pred

  } else if(complex_fit_models=="ksvm"){
    # P(R_overlap|S=1,X)
    fit_Roverlap_rand = cbind.data.frame("overlap"=overlap_rand,X_rand_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_Roverlap_rand = predict(fit_Roverlap_rand, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    # P(R_overlap|S=0,X)
    fit_Roverlap_obs = cbind.data.frame("overlap"=overlap_obs,X_obs_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_Roverlap_obs = predict(fit_Roverlap_obs, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    # P(A|S=1,X)
    fit_A_rand = cbind.data.frame("A"=droplevels(A_rand),X_rand_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_A_rand = predict(fit_A_rand, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    # P(A|S=0,X)
    fit_A_obs = cbind.data.frame("A"=droplevels(A_obs),X_obs_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_A_obs = predict(fit_A_obs, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap = cbind.data.frame("A"=droplevels(A_rand)[overlap_rand==1],
                                          X_rand[overlap_rand==1,]) %>% ksvm(A ~ ., data=.,type="C-svc",prob.model=T)
    pi_A_rand_overlap = predict(fit_A_rand_overlap, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap = cbind.data.frame("A"=droplevels(A_obs)[overlap_obs==1],
                                         X_obs[overlap_obs==1,]) %>% ksvm(A ~ ., data=.,type="C-svc",prob.model=T)
    pi_A_obs_overlap = predict(fit_A_obs_overlap, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    
  } else if(complex_fit_models=="ranger"){
    # P(R_overlap|S=1,X)
    fit_Roverlap_rand = cbind.data.frame("overlap"=overlap_rand,X_rand_matrix) %>% ranger(formula(.), data=.,
                                                                                          num.trees = 300, min.node.size=floor(nrow(X_rand_matrix)*0.05), # these hyperparameters work reasonably in tests
                                                                                          num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_Roverlap_rand = cbind.data.frame(X_target_matrix) %>%
      predict(fit_Roverlap_rand, .,type = "response") %>% .$predictions %>% .[,2]
    
    # P(R_overlap|S=0,X)
    fit_Roverlap_obs = cbind.data.frame("overlap"=overlap_obs,X_obs_matrix) %>% ranger(formula(.), data=.,
                                                                                       num.trees = 300, min.node.size=floor(nrow(X_obs_matrix)*0.05), # these hyperparameters work reasonably in tests
                                                                                       num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_Roverlap_obs = cbind.data.frame(X_target_matrix) %>%
      predict(fit_Roverlap_obs, .,type = "response") %>% .$predictions %>% .[,2]
    
    # P(A|S=1,X)
    fit_A_rand = cbind.data.frame("A"=A_rand,X_rand_matrix) %>% ranger(formula(.), data=.,
                                                                       num.trees = 300, min.node.size=floor(nrow(X_rand_matrix)*0.05), # these hyperparameters work reasonably in tests
                                                                       num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_rand = cbind.data.frame(X_target_matrix) %>%
      predict(fit_A_rand, .,type = "response") %>% .$predictions %>% .[,2]
    
    # P(A|S=0,X)
    fit_A_obs = cbind.data.frame("A"=A_obs,X_obs_matrix) %>% ranger(formula(.), data=.,
                                                                    num.trees = 300, min.node.size=floor(nrow(X_obs_matrix)*0.05), # these hyperparameters work reasonably in tests
                                                                    num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_obs = cbind.data.frame(X_target_matrix) %>%
      predict(fit_A_obs, .,type = "response") %>% .$predictions %>% .[,2]
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap = cbind.data.frame("A"=A_rand[overlap_rand==1],
                                          X_rand[overlap_rand==1,]) %>% ranger(formula(.), data=.,
                                                                               num.trees = 300, min.node.size=floor(nrow(X_rand[overlap_rand==1,])*0.05), # these hyperparameters work reasonably in tests
                                                                               num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_rand_overlap = as.data.frame(X_target_matrix) %>%
      predict(fit_A_rand_overlap, .,type = "response") %>% .$predictions %>% .[,2]
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap = cbind.data.frame("A"=A_obs[overlap_obs==1],
                                         X_obs[overlap_obs==1,]) %>% ranger(formula(.), data=.,
                                                                            num.trees = 300, min.node.size=floor(nrow(X_obs[overlap_obs==1,])*0.05), # these hyperparameters work reasonably in tests
                                                                            num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_obs_overlap = as.data.frame(X_target_matrix) %>%
      predict(fit_A_obs_overlap, .,type = "response") %>% .$predictions %>% .[,2]
    
  } else if(complex_fit_models=="correctly_specified"){
    # P(R_overlap|S=1,X)
    # Warning message: glm.fit: algorithm did not converge - unstable because can almost perfectly predict from X1 so all predictions almost perfectly at 0/1
    fit_Roverlap_rand =  cbind.data.frame("overlap"=overlap_rand,"X1_high"=X1_high_rand) %>% 
      glm(formula(.), data=.,family="binomial")
    pi_Roverlap_rand = predict(fit_Roverlap_rand,
                               as.data.frame(X1_high_target) %>% set_colnames("X1_high"),
                               type = "response")
    
    # P(R_overlap|S=0,X)
    # Warning message: glm.fit: algorithm did not converge - unstable because can almost perfectly predict from X1 so all predictions almost perfectly at 0/1
    fit_Roverlap_obs =  cbind.data.frame("overlap"=overlap_obs,"X1_low"=X1_low_obs) %>% 
      glm(formula(.), data=.,family="binomial")
    pi_Roverlap_obs = predict(fit_Roverlap_obs,
                              as.data.frame(X1_low_target) %>% set_colnames("X1_low"),
                              type = "response")
    
    # P(A|S=1,X) - eliminiate residual noise through estimating although propensity is known
    fit_A_rand =  cbind.data.frame("A"=A_rand,X_rand_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_rand = predict(fit_A_rand,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=0,X)
    fit_A_obs =  cbind.data.frame("A"=A_obs,X_obs_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_obs = predict(fit_A_obs,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap =  cbind.data.frame("A"=A_rand[overlap_rand==1],
                                           X_rand_overlap_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_rand_overlap = predict(fit_A_rand_overlap,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap =  cbind.data.frame("A"=A_obs[overlap_obs==1],
                                          X_obs_overlap_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_obs_overlap = predict(fit_A_obs_overlap,as.data.frame(X_target_matrix),type = "response")
    
  } else {
    # P(R_overlap|S=1,X)
    # Warning message: 1: glm.fit: algorithm did not converge 2: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because can almost perfectly predict from X1 so all predictions almost perfectly at 0/1
    fit_Roverlap_rand =  cbind.data.frame("overlap"=overlap_rand,X_rand_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_Roverlap_rand = predict(fit_Roverlap_rand,as.data.frame(X_target_matrix),type = "response")
    
    # P(R_overlap|S=0,X)
    # Warning message: 1: glm.fit: algorithm did not converge 2: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because can almost perfectly predict from X1 so all predictions almost perfectly at 0/1
    fit_Roverlap_obs =  cbind.data.frame("overlap"=overlap_obs,X_obs_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_Roverlap_obs = predict(fit_Roverlap_obs,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=1,X) - eliminiate residual noise through estimating although propensity is known
    fit_A_rand =  cbind.data.frame("A"=A_rand,X_rand_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_rand = predict(fit_A_rand,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=0,X)
    fit_A_obs =  cbind.data.frame("A"=A_obs,X_obs_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_obs = predict(fit_A_obs,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap =  cbind.data.frame("A"=A_rand[overlap_rand==1],
                                           X_rand_overlap_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_rand_overlap = predict(fit_A_rand_overlap,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap =  cbind.data.frame("A"=A_obs[overlap_obs==1],
                                          X_obs_overlap_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_A_obs_overlap = predict(fit_A_obs_overlap,as.data.frame(X_target_matrix),type = "response")
    
  } 
  
  ### Calculate weights
  w1_A1 = ifelse(S_target==1 & A_target==1,1/(1-pi_A_rand),0)
  w2_A1 = ifelse(S_target==0 & A_target==1,1/(1-pi_A_obs),0)
  w3_A1 = ifelse(S_target==0 & A_target==1 & overlap_target==1,
                 1/(pi_Roverlap_obs*(1-pi_A_obs_overlap)),0)
  w4_A1 = ifelse(S_target==1 & A_target==1 & overlap_target==1,
                 (1-pi_S)/(pi_S*pi_Roverlap_rand*(1-pi_A_rand_overlap)),0) # large weights occur because X1 very strongly predictive and any values smaller than observed RCT overlap values (even if they're in the overlap region) are predicted to be non-overlap
  
  w1_A2 = ifelse(S_target==1 & A_target==2,1/pi_A_rand,0)
  w2_A2 = ifelse(S_target==0 & A_target==2,1/pi_A_obs,0)
  w3_A2 = ifelse(S_target==0 & A_target==2 & overlap_target==1,
                 1/(pi_Roverlap_obs*pi_A_obs_overlap),0)
  w4_A2 = ifelse(S_target==1 & A_target==2 & overlap_target==1,
                 (1-pi_S)/(pi_S*pi_Roverlap_rand*pi_A_rand_overlap),0)
  
  
  ## Trim weights
  # Trim propensities to ensure no near-zero values (should be few due to subsetting to overlap region; will ensure that can calculate mean and variance for weight trimming)
  pi_A_rand_trim = ifelse(pi_A_rand<propensity_trim_threshold,propensity_trim_threshold,
                          ifelse(pi_A_rand>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_A_rand))
  
  pi_A_obs_trim = ifelse(pi_A_obs<propensity_trim_threshold,propensity_trim_threshold,
                         ifelse(pi_A_obs>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_A_obs))
  
  pi_A_rand_overlap_trim = ifelse(pi_A_rand_overlap<propensity_trim_threshold,propensity_trim_threshold,
                                  ifelse(pi_A_rand_overlap>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_A_rand_overlap))
  
  pi_A_obs_overlap_trim = ifelse(pi_A_obs_overlap<propensity_trim_threshold,propensity_trim_threshold,
                                 ifelse(pi_A_obs_overlap>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_A_obs_overlap))
  
  pi_Roverlap_rand_trim = ifelse(pi_Roverlap_rand<propensity_trim_threshold,propensity_trim_threshold,
                                 ifelse(pi_Roverlap_rand>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_Roverlap_rand))
  
  pi_Roverlap_obs_trim = ifelse(pi_Roverlap_obs<propensity_trim_threshold,propensity_trim_threshold,
                                ifelse(pi_Roverlap_obs>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_Roverlap_obs))
  
  ### Trim weight denominators so their product likewise isn't smaller than the propensity threshold
  ### Note: some of the corresponding weights will be set to zero if the observation is not in the prerequisite treatment/selection/overlap group
  # Weight denominators
  w3_denom_A1 = pi_Roverlap_obs_trim*(1-pi_A_obs_overlap_trim)
  w4_denom_A1 = pi_S_trim*pi_Roverlap_rand_trim*(1-pi_A_rand_overlap_trim)
  
  w3_denom_A2 = pi_Roverlap_obs_trim*pi_A_obs_overlap_trim
  w4_denom_A2 = pi_S_trim*pi_Roverlap_rand_trim*pi_A_rand_overlap_trim
  
  # Trim
  w3_denom_A1_trim = ifelse(w3_denom_A1<propensity_trim_threshold,propensity_trim_threshold,
                         ifelse(w3_denom_A1>(1-propensity_trim_threshold),(1-propensity_trim_threshold),w3_denom_A1))
  w4_denom_A1_trim = ifelse(w4_denom_A1<propensity_trim_threshold,propensity_trim_threshold,
                         ifelse(w4_denom_A1>(1-propensity_trim_threshold),(1-propensity_trim_threshold),w4_denom_A1))
  
  w3_denom_A2_trim = ifelse(w3_denom_A2<propensity_trim_threshold,propensity_trim_threshold,
                            ifelse(w3_denom_A2>(1-propensity_trim_threshold),(1-propensity_trim_threshold),w3_denom_A2))
  w4_denom_A2_trim = ifelse(w4_denom_A2<propensity_trim_threshold,propensity_trim_threshold,
                            ifelse(w4_denom_A2>(1-propensity_trim_threshold),(1-propensity_trim_threshold),w4_denom_A2))
  
  ### Calculate trimmed weights
  w1_A1_trim = ifelse(S_target==1 & A_target==1,1/(1-pi_A_rand_trim),0)
  w2_A1_trim = ifelse(S_target==0 & A_target==1,1/(1-pi_A_obs_trim),0)
  w3_A1_trim = ifelse(S_target==0 & A_target==1 & overlap_target==1,
                 1/w3_denom_A1_trim,0)
  w4_A1_trim = ifelse(S_target==1 & A_target==1 & overlap_target==1,
                 (1-pi_S_trim)/w4_denom_A1_trim,0) # large weights occur because X1 very strongly predictive and any values smaller than observed RCT overlap values (even if they're in the overlap region) are predicted to be non-overlap
  
  w1_A2_trim = ifelse(S_target==1 & A_target==2,1/pi_A_rand_trim,0)
  w2_A2_trim = ifelse(S_target==0 & A_target==2,1/pi_A_obs_trim,0)
  w3_A2_trim = ifelse(S_target==0 & A_target==2 & overlap_target==1,
                 1/w3_denom_A2_trim,0)
  w4_A2_trim = ifelse(S_target==1 & A_target==2 & overlap_target==1,
                 (1-pi_S_trim)/w4_denom_A2_trim,0)
  
  
  # Summarize proportion of observations with trimmed weights
  n_trim_w1_A1 = mean(w1_A1 != w1_A1_trim)
  n_trim_w2_A1 = mean(w2_A1 != w2_A1_trim)
  n_trim_w3_A1 = mean(w3_A1 != w3_A1_trim)
  n_trim_w4_A1 = mean(w4_A1 != w4_A1_trim)
  
  n_trim_w1_A2 = mean(w1_A2 != w1_A2_trim)
  n_trim_w2_A2 = mean(w2_A2 != w2_A2_trim)
  n_trim_w3_A2 = mean(w3_A2 != w3_A2_trim)
  n_trim_w4_A2 = mean(w4_A2 != w4_A2_trim)
  
  ##########################################################################################
  ### Model conditional distributions used in estimators 
  ##########################################################################################
  ### Simple case: parametric model that is correctly-specified except unmeasured confounder. But since we have a linear model, it doesn't matter
  if(complex_fit_models=="ensemble"){
    # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
    fit_rand <- SuperLearner(Y=Y_rand, X=cbind.data.frame("A"=A_rand_matrix,X_rand_matrix), 
                             family="gaussian", SL.library=SL_library, verbose=F)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
    fit_obs = SuperLearner(Y=Y_obs, X=cbind.data.frame("A"=A_obs_matrix,X_obs_matrix), 
                           family="gaussian", SL.library=SL_library, verbose=F)
    
    # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
    fit_rand_overlap = SuperLearner(Y=Y_rand[overlap_rand==1], 
                                    X=cbind.data.frame("A"=A_rand_matrix,X_rand_matrix)[overlap_rand==1,], 
                                    family="gaussian", SL.library=SL_library, verbose=F)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
    fit_obs_overlap = SuperLearner(Y=Y_obs[overlap_obs==1], 
                                   X=cbind.data.frame("A"=A_obs_matrix,X_obs_matrix)[overlap_obs==1,], 
                                   family="gaussian", SL.library=SL_library, verbose=F)
    
  } else if(complex_fit_models=="ksvm"){
    # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
    fit_rand <- cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand) %>% 
      ksvm(Y ~ .,data = .)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
    fit_obs = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs) %>%
      ksvm(Y ~ .,data = .) 
    
    # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
    fit_rand_overlap = cbind.data.frame("Y"=Y_rand[overlap_rand==1],
                                        "A"=A_rand[overlap_rand==1],
                                        X_rand[overlap_rand==1,]) %>%
      ksvm(Y ~ .,data = .)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
    fit_obs_overlap = cbind.data.frame("Y"=Y_obs[overlap_obs==1],
                                       "A"=A_obs[overlap_obs==1],
                                       X_obs[overlap_obs==1,]) %>%
      ksvm(Y ~ .,data = .)
    
  } else if(complex_fit_models=="ranger"){
    # ranger: random forests
    # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
    fit_rand <- cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand) %>%
      ranger(Y ~ .,data = .,
             num.trees = 300, min.node.size=floor(nrow(X_rand)*0.05),
             num.threads=1, respect.unordered.factors=T)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
    fit_obs = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs) %>%
      ranger(Y ~ .,data = .,
             num.trees = 300, min.node.size=floor(nrow(X_obs)*0.05),
             num.threads=1, respect.unordered.factors=T)
    
    # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
    fit_rand_overlap = cbind.data.frame("Y"=Y_rand[overlap_rand==1],
                                        "A"=A_rand[overlap_rand==1],
                                        X_rand[overlap_rand==1,]) %>%
      ranger(Y ~ .,data = .,
             num.trees = 300, min.node.size=floor(nrow(X_rand[overlap_rand==1,])*0.05),
             num.threads=1, respect.unordered.factors=T)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
    fit_obs_overlap = cbind.data.frame("Y"=Y_obs[overlap_obs==1],
                                       "A"=A_obs[overlap_obs==1],
                                       X_obs[overlap_obs==1,]) %>%
      ranger(Y ~ .,data = .,
             num.trees = 300, min.node.size=floor(nrow(X_obs[overlap_obs==1,])*0.05),
             num.threads=1, respect.unordered.factors=T)
    
    
  } else if(complex_fit_models=="correctly_specified" & is.null(target_sample$X1_cube)){ # Linear model fits are the only ones missing X1_cube from the correctly_specified target sample data
    # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
    fit_rand = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
      lm(Y ~ . + X1:A,data=.) 
    # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
    fit_obs = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
      lm(Y ~ . + X1:A,data=.) 
    
    # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
    fit_rand_overlap = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
      lm(Y ~ . + X1:A,data=., subset=(overlap_rand==1)) 
    
    # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
    fit_obs_overlap = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
      lm(Y ~ . + X1:A,data=., subset=(overlap_obs==1)) 
    
  } else if(complex_fit_models=="correctly_specified" & is.null(target_sample$X1_knot)){ # Overfit by potentially only X1_cube term
    # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
    fit_rand = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
      lm(Y ~ . + X1:A + X1_cube:A,data=.) 
    # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
    fit_obs = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
      lm(Y ~ . + X1:A + X1_cube:A,data=.) 
    
    # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
    fit_rand_overlap = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
      lm(Y ~ . + X1:A + X1_cube:A,data=., subset=(overlap_rand==1)) 
    
    # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
    fit_obs_overlap = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
      lm(Y ~ . + X1:A + X1_cube:A,data=., subset=(overlap_obs==1)) 
    
  } else { # Include interactions between all covariates and A
    # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
    fit_rand = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
      lm(Y ~ . + .:A,data=.) 
    
    # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
    fit_obs = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
      lm(Y ~ . + .:A,data=.) 
    
    # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
    fit_rand_overlap = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
      lm(Y ~ . + .:A,data=., subset=(overlap_rand==1)) 
    
    # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
    fit_obs_overlap = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
      lm(Y ~ . + .:A,data=., subset=(overlap_obs==1)) 
    
  }  
  
  
  
  
  ##########################################################################################
  ### Set up target data to estimate potential outcomes
  ##########################################################################################
  # Set up target data with the same variable names for prediction
  #    A2    A3     X     
  data_A1_target = cbind.data.frame("A"=levels_A[1],X_target) # target data with A=1 for all observations
  data_A2_target = cbind.data.frame("A"=levels_A[2],X_target) # target data with A=2 for all observations
  data_A3_target = cbind.data.frame("A"=levels_A[3],X_target) # target data with A=3 for all observations
  
  if(complex_fit_models=="ensemble"){
    data_A1_target_matrix = cbind.data.frame("A"=as.numeric(levels_A[1])-1,X_target_matrix) # target data with A=0 for all observations
    data_A2_target_matrix = cbind.data.frame("A"=as.numeric(levels_A[2])-1,X_target_matrix) # target data with A=1 for all observations
    data_A3_target_matrix = cbind.data.frame("A"=as.numeric(levels_A[3])-1,X_target_matrix) # target data with A=2 for all observations
  } else {
    data_A1_target_matrix = cbind.data.frame("A"=levels_A[1],X_target_matrix) # target data with A=1 for all observations
    data_A2_target_matrix = cbind.data.frame("A"=levels_A[2],X_target_matrix) # target data with A=2 for all observations
    data_A3_target_matrix = cbind.data.frame("A"=levels_A[3],X_target_matrix) # target data with A=3 for all observations
  }

  #               X      
  data_X_obs = as.data.frame(cbind(X_obs))
  data_X_obs_matrix = model.matrix(~.,data_X_obs)[,-1]
  
  data_X1to4_obs = as.data.frame(cbind(X1to4_obs))
  data_X1to4_obs_matrix = model.matrix(~.,data_X1to4_obs)[,-1]
  
  
  ## A1
  # Needed for CCDS
  # Estimator name indicates data being used (e.g., obs_overlap) and model from which potential outcomes will be estimated (e.g., for_fit_rand). The latter is needed so variable names match for predicting from the model
  data_A1_rand = data_A1_target[S_target==1,]
  data_A1_obs = data_A1_target[S_target==0,]
  
  data_A1_rand_matrix = data_A1_target_matrix[S_target==1,]
  data_A1_obs_matrix = data_A1_target_matrix[S_target==0,]
  
  # Needed for 2-stage WD
  data_A1_rand_overlap = data_A1_target[S_target==1 & overlap_target==1,]
  data_A1_rand_matrix_overlap = data_A1_target_matrix[S_target==1 & overlap_target==1,]
  
  ## A2
  # Needed for CCDS
  data_A2_rand = data_A2_target[S_target==1,]
  data_A2_obs = data_A2_target[S_target==0,]
  
  data_A2_rand_matrix = data_A2_target_matrix[S_target==1,]
  data_A2_obs_matrix = data_A2_target_matrix[S_target==0,]
  
  # Needed for 2-stage WD
  data_A2_rand_overlap = data_A2_target[S_target==1 & overlap_target==1,]
  data_A2_rand_matrix_overlap = data_A2_target_matrix[S_target==1 & overlap_target==1,]
  
  ## A3
  # Needed for CCDS
  data_A3_rand = data_A3_target[S_target==1,]
  data_A3_obs = data_A3_target[S_target==0,]
  
  data_A3_rand_matrix = data_A3_target_matrix[S_target==1,]
  data_A3_obs_matrix = data_A3_target_matrix[S_target==0,]
  
  # Needed for 2-stage WD
  data_A3_rand_overlap = data_A3_target[S_target==1 & overlap_target==1,]
  data_A3_rand_matrix_overlap = data_A3_target_matrix[S_target==1 & overlap_target==1,]
  
  
  ##########################################################################################
  ### Estimate potential outcomes: CCDS 
  ##########################################################################################
  if(complex_fit_models=="ensemble"){
    # Y1
    Y1_rand_pred = predict(fit_rand,data_A1_rand_matrix, onlySL = TRUE)$pred 
    Y1_obs_pred = predict(fit_obs,data_A1_obs_matrix, onlySL = TRUE)$pred
    Y1_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A1_obs_matrix, onlySL = TRUE)$pred
    Y1_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A1_obs_matrix, onlySL = TRUE)$pred
    Y1_obs_pred_bias_CCDS = Y1_obs_pred_bias_CCDSa - Y1_obs_pred_bias_CCDSb  # bias estimate for obs data
    Y1_obs_pred_debiased_CCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Y1_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A1_rand_matrix, onlySL = TRUE)$pred 
    
    # Y2
    Y2_rand_pred = predict(fit_rand,data_A2_rand_matrix, onlySL = TRUE)$pred
    Y2_obs_pred = predict(fit_obs,data_A2_obs_matrix, onlySL = TRUE)$pred
    Y2_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A2_obs_matrix, onlySL = TRUE)$pred
    Y2_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A2_obs_matrix, onlySL = TRUE)$pred
    Y2_obs_pred_bias_CCDS = Y2_obs_pred_bias_CCDSa - Y2_obs_pred_bias_CCDSb  # bias estimate for obs data
    Y2_obs_pred_debiased_CCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Y2_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A2_rand_matrix, onlySL = TRUE)$pred 
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred = predict(fit_rand,data_A3_rand_matrix, onlySL = TRUE)$pred
      Y3_obs_pred = predict(fit_obs,data_A3_obs_matrix, onlySL = TRUE)$pred
      Y3_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A3_obs_matrix, onlySL = TRUE)$pred
      Y3_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A3_obs_matrix, onlySL = TRUE)$pred
      Y3_obs_pred_bias_CCDS = Y3_obs_pred_bias_CCDSa - Y3_obs_pred_bias_CCDSb  # bias estimate for obs data
      Y3_obs_pred_debiased_CCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_CCDS # bias estimate for obs data
      
      # Used in CCDS AIPW
      Y3_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A3_rand_matrix, onlySL = TRUE)$pred
    }
    
  } else if(complex_fit_models=="ranger"){
    # Y1
    Y1_rand_pred = predict(fit_rand,data_A1_rand)$predictions
    Y1_obs_pred = predict(fit_obs,data_A1_obs)$predictions
    Y1_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A1_obs)$predictions
    Y1_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A1_obs)$predictions
    Y1_obs_pred_bias_CCDS = Y1_obs_pred_bias_CCDSa - Y1_obs_pred_bias_CCDSb  # bias estimate for obs data
    Y1_obs_pred_debiased_CCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Y1_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A1_rand)$predictions 
    
    # Y2
    Y2_rand_pred = predict(fit_rand,data_A2_rand)$predictions
    Y2_obs_pred = predict(fit_obs,data_A2_obs)$predictions
    Y2_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A2_obs)$predictions
    Y2_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A2_obs)$predictions
    Y2_obs_pred_bias_CCDS = Y2_obs_pred_bias_CCDSa - Y2_obs_pred_bias_CCDSb  # bias estimate for obs data
    Y2_obs_pred_debiased_CCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Y2_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A2_rand_matrix)$predictions 
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred = predict(fit_rand,data_A3_rand)$predictions
      Y3_obs_pred = predict(fit_obs,data_A3_obs)$predictions
      Y3_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A3_obs)$predictions
      Y3_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A3_obs)$predictions
      Y3_obs_pred_bias_CCDS = Y3_obs_pred_bias_CCDSa - Y3_obs_pred_bias_CCDSb  # bias estimate for obs data
      Y3_obs_pred_debiased_CCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_CCDS # bias estimate for obs data
      
      # Used in CCDS AIPW
      Y3_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A3_rand_matrix)$predictions
    }
    
  } else {
    # Warnings arise because of collinearity due to data generating mechanism
    # Y1
    Y1_rand_pred = predict(fit_rand,data_A1_rand_matrix)
    Y1_obs_pred = predict(fit_obs,data_A1_obs_matrix)
    Y1_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A1_obs_matrix)
    Y1_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A1_obs_matrix) 
    Y1_obs_pred_bias_CCDS = Y1_obs_pred_bias_CCDSa - Y1_obs_pred_bias_CCDSb  # bias estimate for obs data
    Y1_obs_pred_debiased_CCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Y1_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A1_rand_matrix) 
    
    # Y2
    Y2_rand_pred = predict(fit_rand,data_A2_rand_matrix)
    Y2_obs_pred = predict(fit_obs,data_A2_obs_matrix)
    Y2_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A2_obs_matrix)
    Y2_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A2_obs_matrix) 
    Y2_obs_pred_bias_CCDS = Y2_obs_pred_bias_CCDSa - Y2_obs_pred_bias_CCDSb  # bias estimate for obs data
    Y2_obs_pred_debiased_CCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Y2_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A2_rand_matrix) 
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred = predict(fit_rand,data_A3_rand_matrix)
      Y3_obs_pred = predict(fit_obs,data_A3_obs_matrix)
      Y3_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A3_obs_matrix)
      Y3_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A3_obs_matrix) 
      Y3_obs_pred_bias_CCDS = Y3_obs_pred_bias_CCDSa - Y3_obs_pred_bias_CCDSb  # bias estimate for obs data
      Y3_obs_pred_debiased_CCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_CCDS # bias estimate for obs data
      
      # Used in CCDS AIPW
      Y3_rand_pred_bias_CCDSb = predict(fit_rand_overlap,data_A3_rand_matrix) 
    }
    
  }
  
  # Target population estimates
  Y1_target_pred_CCDS = rep(NA,n_target)
  Y1_target_pred_CCDS[S_target==1] = Y1_rand_pred
  Y1_target_pred_CCDS[S_target==0] = Y1_obs_pred_debiased_CCDS
  
  Y2_target_pred_CCDS = rep(NA,n_target)
  Y2_target_pred_CCDS[S_target==1] = Y2_rand_pred
  Y2_target_pred_CCDS[S_target==0] = Y2_obs_pred_debiased_CCDS
  
  Y3_target_pred_CCDS = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_CCDS[S_target==1] = Y3_rand_pred
    Y3_target_pred_CCDS[S_target==0] = Y3_obs_pred_debiased_CCDS
  }
  
  # Mean potential outcome estimates
  mean_Y1_target_pred_CCDS = mean(Y1_target_pred_CCDS)
  mean_Y2_target_pred_CCDS = mean(Y2_target_pred_CCDS)
  mean_Y3_target_pred_CCDS = mean(Y3_target_pred_CCDS)
  
  # True mean potential outcomes
  mean_Y1_target = mean(Y1_target)
  mean_Y2_target = mean(Y2_target)
  mean_Y3_target = mean(Y3_target)
  
  ##########################################################################################
  ### Estimate potential outcomes: 2-stage CCDS and weighted 2-stage CCDS
  ##########################################################################################
  ### Weight for 2nd stage bias term
  ## Trim weight denominators so their product likewise isn't smaller than the propensity threshold
  # Weight denominators
  w_bias_denom = pi_Roverlap_rand_trim*pi_S_trim
  
  # Trim
  w_bias_denom_trim = ifelse(w_bias_denom<propensity_trim_threshold,propensity_trim_threshold,
                             ifelse(w_bias_denom>(1-propensity_trim_threshold),(1-propensity_trim_threshold),w_bias_denom))
  
  # Weight
  w_bias = ifelse(S_target==1 & overlap_target==1,(1-pi_S)/(pi_Roverlap_rand*pi_S),0)
  w_bias_trim = ifelse(S_target==1 & overlap_target==1,(1-pi_S_trim)/w_bias_denom_trim,0)
  
  # Amount of trimming
  n_trim_w_bias = mean(w_bias != w_bias_trim)
  
  # Normalize weights
  w_bias_norm = w_bias_trim/sum(w_bias_trim)
  w_bias_norm_rand_overlap = w_bias_norm[S_target==1 & overlap_target==1]
  
  ### hat{phi}_1 with second-stage projection of bias estimates onto X: predict counterfactuals for target population
  ### Uses randomized data in overlap region for prediction of bias term
  if(complex_fit_models=="ensemble"){
    # Y1
    Y1_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A1_rand_matrix_overlap, onlySL = TRUE)$pred
    Y1_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A1_rand_matrix_overlap, onlySL = TRUE)$pred
    Y1_rand_pred_bias_overlap_2stageCCDS = Y1_rand_pred_bias_overlap_2stageCCDS_a-
      Y1_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y1_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y1_rand_overlap_2stageCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=.)
    
    Y1_obs_pred_bias_2stageCCDS = predict(fit_Y1_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_2stageCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_2stageCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y1_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y1_obs_pred_bias_weighted2stageCCDS = predict(fit_Y1_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_weighted2stageCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    
    # Y2
    Y2_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A2_rand_matrix_overlap, onlySL = TRUE)$pred
    Y2_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A2_rand_matrix_overlap, onlySL = TRUE)$pred
    Y2_rand_pred_bias_overlap_2stageCCDS = Y2_rand_pred_bias_overlap_2stageCCDS_a-
      Y2_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y2_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y2_rand_overlap_2stageCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=.)
    
    Y2_obs_pred_bias_2stageCCDS = predict(fit_Y2_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_2stageCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_2stageCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y2_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y2_obs_pred_bias_weighted2stageCCDS = predict(fit_Y2_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_weighted2stageCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A3_rand_matrix_overlap, onlySL = TRUE)$pred
      Y3_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A3_rand_matrix_overlap, onlySL = TRUE)$pred
      Y3_rand_pred_bias_overlap_2stageCCDS = Y3_rand_pred_bias_overlap_2stageCCDS_a-
        Y3_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
      
      # Estimate restricted mean model E(Y3_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
      fit_Y3_rand_overlap_2stageCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=.)
      
      Y3_obs_pred_bias_2stageCCDS = predict(fit_Y3_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_2stageCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_2stageCCDS # bias estimate for obs data
      
      # Weighted version
      fit_Y3_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
      
      Y3_obs_pred_bias_weighted2stageCCDS = predict(fit_Y3_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_weighted2stageCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    }
    
  } else if(complex_fit_models=="ranger"){
    # Y1
    Y1_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A1_rand_overlap)$predictions
    Y1_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A1_rand_overlap)$predictions
    Y1_rand_pred_bias_overlap_2stageCCDS = Y1_rand_pred_bias_overlap_2stageCCDS_a-
      Y1_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y1_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y1_rand_overlap_2stageCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=.)
    
    Y1_obs_pred_bias_2stageCCDS = predict(fit_Y1_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_2stageCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_2stageCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y1_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y1_obs_pred_bias_weighted2stageCCDS = predict(fit_Y1_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_weighted2stageCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    
    # Y2
    Y2_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A2_rand_overlap)$predictions
    Y2_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A2_rand_overlap)$predictions
    Y2_rand_pred_bias_overlap_2stageCCDS = Y2_rand_pred_bias_overlap_2stageCCDS_a-
      Y2_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y2_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y2_rand_overlap_2stageCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=.)
    
    Y2_obs_pred_bias_2stageCCDS = predict(fit_Y2_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_2stageCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_2stageCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y2_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y2_obs_pred_bias_weighted2stageCCDS = predict(fit_Y2_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_weighted2stageCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A3_rand_overlap)$predictions
      Y3_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A3_rand_overlap)$predictions
      Y3_rand_pred_bias_overlap_2stageCCDS = Y3_rand_pred_bias_overlap_2stageCCDS_a-
        Y3_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
      
      # Estimate restricted mean model E(Y3_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
      fit_Y3_rand_overlap_2stageCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=.)
      
      Y3_obs_pred_bias_2stageCCDS = predict(fit_Y3_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_2stageCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_2stageCCDS # bias estimate for obs data
      
      # Weighted version
      fit_Y3_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
      
      Y3_obs_pred_bias_weighted2stageCCDS = predict(fit_Y3_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_weighted2stageCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    }
    
  } else {
    # Y1
    Y1_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A1_rand_matrix_overlap)
    Y1_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A1_rand_matrix_overlap)
    Y1_rand_pred_bias_overlap_2stageCCDS = Y1_rand_pred_bias_overlap_2stageCCDS_a-
      Y1_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y1_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y1_rand_overlap_2stageCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>% 
      lm(Y1_bias ~ ., data=.)
    
    Y1_obs_pred_bias_2stageCCDS = predict(fit_Y1_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_2stageCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_2stageCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y1_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y1_obs_pred_bias_weighted2stageCCDS = predict(fit_Y1_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_weighted2stageCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    
    # Y2
    Y2_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A2_rand_matrix_overlap)
    Y2_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A2_rand_matrix_overlap)
    Y2_rand_pred_bias_overlap_2stageCCDS = Y2_rand_pred_bias_overlap_2stageCCDS_a-
      Y2_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y2_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y2_rand_overlap_2stageCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>% 
      lm(Y2_bias ~ ., data=.)
    
    Y2_obs_pred_bias_2stageCCDS = predict(fit_Y2_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_2stageCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_2stageCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y2_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y2_obs_pred_bias_weighted2stageCCDS = predict(fit_Y2_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_weighted2stageCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred_bias_overlap_2stageCCDS_a = predict(fit_obs_overlap,data_A3_rand_matrix_overlap)
      Y3_rand_pred_bias_overlap_2stageCCDS_b = predict(fit_rand_overlap,data_A3_rand_matrix_overlap)
      Y3_rand_pred_bias_overlap_2stageCCDS = Y3_rand_pred_bias_overlap_2stageCCDS_a-
        Y3_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
      
      # Estimate restricted mean model E(Y3_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the observational data in the overlap region
      fit_Y3_rand_overlap_2stageCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>% 
        lm(Y3_bias ~ ., data=.)
      
      Y3_obs_pred_bias_2stageCCDS = predict(fit_Y3_rand_overlap_2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_2stageCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_2stageCCDS # bias estimate for obs data
      
      # Weighted version
      fit_Y3_rand_overlap_weighted2stageCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
      
      Y3_obs_pred_bias_weighted2stageCCDS = predict(fit_Y3_rand_overlap_weighted2stageCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_weighted2stageCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data
    }
    
  }
  
  # Target population estimates
  Y1_target_pred_2stageCCDS = rep(NA,n_target)
  Y1_target_pred_2stageCCDS[S_target==1] = Y1_rand_pred
  Y1_target_pred_2stageCCDS[S_target==0] = Y1_obs_pred_debiased_2stageCCDS
  
  Y2_target_pred_2stageCCDS = rep(NA,n_target)
  Y2_target_pred_2stageCCDS[S_target==1] = Y2_rand_pred
  Y2_target_pred_2stageCCDS[S_target==0] = Y2_obs_pred_debiased_2stageCCDS
  
  Y3_target_pred_2stageCCDS = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_2stageCCDS[S_target==1] = Y3_rand_pred
    Y3_target_pred_2stageCCDS[S_target==0] = Y3_obs_pred_debiased_2stageCCDS
  }
  
  # Mean potential outcome estimates
  mean_Y1_target_pred_2stageCCDS = mean(Y1_target_pred_2stageCCDS)
  mean_Y2_target_pred_2stageCCDS = mean(Y2_target_pred_2stageCCDS)
  mean_Y3_target_pred_2stageCCDS = mean(Y3_target_pred_2stageCCDS)
  
  
  # Weighted target population estimates
  Y1_target_pred_weighted2stageCCDS = rep(NA,n_target)
  Y1_target_pred_weighted2stageCCDS[S_target==1] = Y1_rand_pred
  Y1_target_pred_weighted2stageCCDS[S_target==0] = Y1_obs_pred_debiased_weighted2stageCCDS
  
  Y2_target_pred_weighted2stageCCDS = rep(NA,n_target)
  Y2_target_pred_weighted2stageCCDS[S_target==1] = Y2_rand_pred
  Y2_target_pred_weighted2stageCCDS[S_target==0] = Y2_obs_pred_debiased_weighted2stageCCDS
  
  Y3_target_pred_weighted2stageCCDS = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_weighted2stageCCDS[S_target==1] = Y3_rand_pred
    Y3_target_pred_weighted2stageCCDS[S_target==0] = Y3_obs_pred_debiased_weighted2stageCCDS
  }
  
  # Mean potential outcome estimates
  mean_Y1_target_pred_weighted2stageCCDS = mean(Y1_target_pred_weighted2stageCCDS)
  mean_Y2_target_pred_weighted2stageCCDS = mean(Y2_target_pred_weighted2stageCCDS)
  mean_Y3_target_pred_weighted2stageCCDS = mean(Y3_target_pred_weighted2stageCCDS)
  
  
  
  ##########################################################################################
  ### Estimate potential outcomes: 2-stage hybrid CCDS and weighted 2-stage hybrid CCDS
  ##########################################################################################
  ### Use model fits from entire RCT and entire obs data but estimate bias term via predictions 
  ###   for randomized data in overlap region
  if(complex_fit_models=="ensemble"){
    # Y1
    Y1_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A1_rand_matrix_overlap, onlySL = TRUE)$pred
    Y1_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A1_rand_matrix_overlap, onlySL = TRUE)$pred
    Y1_rand_pred_bias_overlap_2stageHCCDS = Y1_rand_pred_bias_overlap_2stageHCCDS_a-
      Y1_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y1_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y1_rand_overlap_2stageHCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=.)
    
    Y1_obs_pred_bias_2stageHCCDS = predict(fit_Y1_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_2stageHCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_2stageHCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y1_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y1_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y1_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_weighted2stageHCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    
    # Y2
    Y2_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A2_rand_matrix_overlap, onlySL = TRUE)$pred
    Y2_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A2_rand_matrix_overlap, onlySL = TRUE)$pred
    Y2_rand_pred_bias_overlap_2stageHCCDS = Y2_rand_pred_bias_overlap_2stageHCCDS_a-
      Y2_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y2_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y2_rand_overlap_2stageHCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=.)
    
    Y2_obs_pred_bias_2stageHCCDS = predict(fit_Y2_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_2stageHCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_2stageHCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y2_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y2_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y2_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_weighted2stageHCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A3_rand_matrix_overlap, onlySL = TRUE)$pred
      Y3_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A3_rand_matrix_overlap, onlySL = TRUE)$pred
      Y3_rand_pred_bias_overlap_2stageHCCDS = Y3_rand_pred_bias_overlap_2stageHCCDS_a-
        Y3_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
      
      # Estimate restricted mean model E(Y3_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
      fit_Y3_rand_overlap_2stageHCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=.)
      
      Y3_obs_pred_bias_2stageHCCDS = predict(fit_Y3_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_2stageHCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_2stageHCCDS # bias estimate for obs data
      
      # Weighted version
      fit_Y3_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
      
      Y3_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y3_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_weighted2stageHCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    }
    
  } else if(complex_fit_models=="ranger"){
    # Y1
    Y1_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A1_rand_overlap)$predictions
    Y1_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A1_rand_overlap)$predictions
    Y1_rand_pred_bias_overlap_2stageHCCDS = Y1_rand_pred_bias_overlap_2stageHCCDS_a-
      Y1_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y1_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y1_rand_overlap_2stageHCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=.)
    
    Y1_obs_pred_bias_2stageHCCDS = predict(fit_Y1_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_2stageHCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_2stageHCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y1_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y1_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y1_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_weighted2stageHCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    
    # Y2
    Y2_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A2_rand_overlap)$predictions
    Y2_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A2_rand_overlap)$predictions
    Y2_rand_pred_bias_overlap_2stageHCCDS = Y2_rand_pred_bias_overlap_2stageHCCDS_a-
      Y2_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y2_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y2_rand_overlap_2stageHCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=.)
    
    Y2_obs_pred_bias_2stageHCCDS = predict(fit_Y2_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_2stageHCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_2stageHCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y2_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y2_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y2_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_weighted2stageHCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A3_rand_overlap)$predictions
      Y3_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A3_rand_overlap)$predictions
      Y3_rand_pred_bias_overlap_2stageHCCDS = Y3_rand_pred_bias_overlap_2stageHCCDS_a-
        Y3_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
      
      # Estimate restricted mean model E(Y3_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
      fit_Y3_rand_overlap_2stageHCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=.)
      
      Y3_obs_pred_bias_2stageHCCDS = predict(fit_Y3_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_2stageHCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_2stageHCCDS # bias estimate for obs data
      
      # Weighted version
      fit_Y3_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
      
      Y3_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y3_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_weighted2stageHCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    }
    
  } else {
    # Y1
    Y1_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A1_rand_matrix_overlap)
    Y1_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A1_rand_matrix_overlap)
    Y1_rand_pred_bias_overlap_2stageHCCDS = Y1_rand_pred_bias_overlap_2stageHCCDS_a-
      Y1_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y1_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y1_rand_overlap_2stageHCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>% 
      lm(Y1_bias ~ ., data=.)
    
    Y1_obs_pred_bias_2stageHCCDS = predict(fit_Y1_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_2stageHCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_2stageHCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y1_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y1_bias"=Y1_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y1_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y1_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y1_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y1_obs_pred_debiased_weighted2stageHCCDS = Y1_obs_pred - # obs estimate
      Y1_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    
    # Y2
    Y2_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A2_rand_matrix_overlap)
    Y2_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A2_rand_matrix_overlap)
    Y2_rand_pred_bias_overlap_2stageHCCDS = Y2_rand_pred_bias_overlap_2stageHCCDS_a-
      Y2_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
    
    # Estimate restricted mean model E(Y2_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
    fit_Y2_rand_overlap_2stageHCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>% 
      lm(Y2_bias ~ ., data=.)
    
    Y2_obs_pred_bias_2stageHCCDS = predict(fit_Y2_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_2stageHCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_2stageHCCDS # bias estimate for obs data
    
    # Weighted version
    fit_Y2_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y2_bias"=Y2_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
      lm(Y2_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
    
    Y2_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y2_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
    
    Y2_obs_pred_debiased_weighted2stageHCCDS = Y2_obs_pred - # obs estimate
      Y2_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    
    # Y3
    if(!is.na(levels_A[3])){
      Y3_rand_pred_bias_overlap_2stageHCCDS_a = predict(fit_obs,data_A3_rand_matrix_overlap)
      Y3_rand_pred_bias_overlap_2stageHCCDS_b = predict(fit_rand,data_A3_rand_matrix_overlap)
      Y3_rand_pred_bias_overlap_2stageHCCDS = Y3_rand_pred_bias_overlap_2stageHCCDS_a-
        Y3_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
      
      # Estimate restricted mean model E(Y3_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
      fit_Y3_rand_overlap_2stageHCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>% 
        lm(Y3_bias ~ ., data=.)
      
      Y3_obs_pred_bias_2stageHCCDS = predict(fit_Y3_rand_overlap_2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_2stageHCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_2stageHCCDS # bias estimate for obs data
      
      # Weighted version
      fit_Y3_rand_overlap_weighted2stageHCCDS = cbind.data.frame("Y3_bias"=Y3_rand_pred_bias_overlap_2stageHCCDS,X1to4_rand_overlap_matrix) %>%
        lm(Y3_bias ~ ., data=., weights=w_bias_norm_rand_overlap)
      
      Y3_obs_pred_bias_weighted2stageHCCDS = predict(fit_Y3_rand_overlap_weighted2stageHCCDS,as.data.frame(data_X1to4_obs_matrix)) # bias estimate for obs data
      
      Y3_obs_pred_debiased_weighted2stageHCCDS = Y3_obs_pred - # obs estimate
        Y3_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
    }
  }
  
  # Target population estimates
  Y1_target_pred_2stageHCCDS = rep(NA,n_target)
  Y1_target_pred_2stageHCCDS[S_target==1] = Y1_rand_pred
  Y1_target_pred_2stageHCCDS[S_target==0] = Y1_obs_pred_debiased_2stageHCCDS
  
  Y2_target_pred_2stageHCCDS = rep(NA,n_target)
  Y2_target_pred_2stageHCCDS[S_target==1] = Y2_rand_pred
  Y2_target_pred_2stageHCCDS[S_target==0] = Y2_obs_pred_debiased_2stageHCCDS
  
  Y3_target_pred_2stageHCCDS = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_2stageHCCDS[S_target==1] = Y3_rand_pred
    Y3_target_pred_2stageHCCDS[S_target==0] = Y3_obs_pred_debiased_2stageHCCDS
  }
  
  # Mean potential outcome estimates
  mean_Y1_target_pred_2stageHCCDS = mean(Y1_target_pred_2stageHCCDS)
  mean_Y2_target_pred_2stageHCCDS = mean(Y2_target_pred_2stageHCCDS)
  mean_Y3_target_pred_2stageHCCDS = mean(Y3_target_pred_2stageHCCDS)
  
  
  # Weighted target population estimates
  Y1_target_pred_weighted2stageHCCDS = rep(NA,n_target)
  Y1_target_pred_weighted2stageHCCDS[S_target==1] = Y1_rand_pred
  Y1_target_pred_weighted2stageHCCDS[S_target==0] = Y1_obs_pred_debiased_weighted2stageHCCDS
  
  Y2_target_pred_weighted2stageHCCDS = rep(NA,n_target)
  Y2_target_pred_weighted2stageHCCDS[S_target==1] = Y2_rand_pred
  Y2_target_pred_weighted2stageHCCDS[S_target==0] = Y2_obs_pred_debiased_weighted2stageHCCDS
  
  Y3_target_pred_weighted2stageHCCDS = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_weighted2stageHCCDS[S_target==1] = Y3_rand_pred
    Y3_target_pred_weighted2stageHCCDS[S_target==0] = Y3_obs_pred_debiased_weighted2stageHCCDS
  }
  
  # Mean potential outcome estimates
  mean_Y1_target_pred_weighted2stageHCCDS = mean(Y1_target_pred_weighted2stageHCCDS)
  mean_Y2_target_pred_weighted2stageHCCDS = mean(Y2_target_pred_weighted2stageHCCDS)
  mean_Y3_target_pred_weighted2stageHCCDS = mean(Y3_target_pred_weighted2stageHCCDS)
  
  
  
  
  ##########################################################################################
  ### Estimate potential outcomes: CCDS-IPW 
  ##########################################################################################
  # Y1
  mean_Y1_rand_pred_CCDS_IPW = sum(w1_A1_trim*Y_target)/sum(w1_A1_trim)
  mean_Y1_obs_pred_CCDS_IPW = sum(w2_A1_trim*Y_target)/sum(w2_A1_trim)
  mean_Y1_obs_pred_bias_CCDS_IPWa = sum(w3_A1_trim*Y_target)/sum(w3_A1_trim)
  mean_Y1_obs_pred_bias_CCDS_IPWb = sum(w4_A1_trim*Y_target)/sum(w4_A1_trim) 
  mean_Y1_obs_pred_bias_CCDS_IPW = mean_Y1_obs_pred_bias_CCDS_IPWa - mean_Y1_obs_pred_bias_CCDS_IPWb  # bias estimate for obs data
  mean_Y1_obs_pred_CCDS_IPW_debiased = mean_Y1_obs_pred_CCDS_IPW - # obs estimate
    mean_Y1_obs_pred_bias_CCDS_IPW # bias estimate for obs data
  
  # Y2
  mean_Y2_rand_pred_CCDS_IPW = sum(w1_A2_trim*Y_target)/sum(w1_A2_trim)
  mean_Y2_obs_pred_CCDS_IPW = sum(w2_A2_trim*Y_target)/sum(w2_A2_trim)
  mean_Y2_obs_pred_bias_CCDS_IPWa = sum(w3_A2_trim*Y_target)/sum(w3_A2_trim)
  mean_Y2_obs_pred_bias_CCDS_IPWb = sum(w4_A2_trim*Y_target)/sum(w4_A2_trim) 
  mean_Y2_obs_pred_bias_CCDS_IPW = mean_Y2_obs_pred_bias_CCDS_IPWa - mean_Y2_obs_pred_bias_CCDS_IPWb  # bias estimate for obs data
  mean_Y2_obs_pred_CCDS_IPW_debiased = mean_Y2_obs_pred_CCDS_IPW - # obs estimate
    mean_Y2_obs_pred_bias_CCDS_IPW # bias estimate for obs data
  
  # Y3
  mean_Y3_rand_pred_CCDS_IPW = mean_Y3_obs_pred_CCDS_IPW_debiased = NA
  if(!is.na(levels_A[3])){
    mean_Y3_rand_pred_CCDS_IPW = sum(w1_A3_trim*Y_target)/sum(w1_A3_trim)
    mean_Y3_obs_pred_CCDS_IPW = sum(w2_A3_trim*Y_target)/sum(w2_A3_trim)
    mean_Y3_obs_pred_bias_CCDS_IPWa = sum(w3_A3_trim*Y_target)/sum(w3_A3_trim)
    mean_Y3_obs_pred_bias_CCDS_IPWb = sum(w4_A3_trim*Y_target)/sum(w4_A3_trim) 
    mean_Y3_obs_pred_bias_CCDS_IPW = mean_Y3_obs_pred_bias_CCDS_IPWa - mean_Y3_obs_pred_bias_CCDS_IPWb  # bias estimate for obs data
    mean_Y3_obs_pred_CCDS_IPW_debiased = mean_Y3_obs_pred_CCDS_IPW - # obs estimate
      mean_Y3_obs_pred_bias_CCDS_IPW # bias estimate for obs data
  }
  
  # Target population estimates
  mean_Y1_target_pred_CCDS_IPW = n_rand/n_target*mean_Y1_rand_pred_CCDS_IPW+
    n_obs/n_target*mean_Y1_obs_pred_CCDS_IPW_debiased
  mean_Y2_target_pred_CCDS_IPW = n_rand/n_target*mean_Y2_rand_pred_CCDS_IPW+
    n_obs/n_target*mean_Y2_obs_pred_CCDS_IPW_debiased
  mean_Y3_target_pred_CCDS_IPW = n_rand/n_target*mean_Y3_rand_pred_CCDS_IPW+
    n_obs/n_target*mean_Y3_obs_pred_CCDS_IPW_debiased
  
  ##########################################################################################
  ### Estimate potential outcomes: CCDS-AIPW 
  ##########################################################################################
  # Y1
  Y1_rand_pred_CCDS_AIPW = n_rand/n_target*w1_A1_trim[S_target==1]/sum(w1_A1_trim)*(Y_rand-Y1_rand_pred) + Y1_rand_pred # pre-subsetting to RCT data; weights then subset to RCT A=1 data
  Y1_obs_pred_CCDS_AIPW =  n_obs/n_target*w2_A1_trim[S_target==0]/sum(w2_A1_trim)*(Y_obs-Y1_obs_pred) + Y1_obs_pred # pre-subsetting to obs data
  Y1_obs_pred_bias_CCDS_AIPWa = n_obs/n_target*w3_A1_trim[S_target==0]/sum(w3_A1_trim)*(Y_obs-Y1_obs_pred_bias_CCDSa) + Y1_obs_pred_bias_CCDSa # pre-subsetting to obs data
  Y1_pred_bias_CCDS_AIPWb = ifelse(S_target==1,n_obs/n_target*w4_A1_trim[S_target==1]/sum(w4_A1_trim)*(Y_rand-Y1_rand_pred_bias_CCDSb),
                                   Y1_obs_pred_bias_CCDSb) # combo of RCT (first component) and obs (second component) data
  
  
  # Y2
  Y2_rand_pred_CCDS_AIPW = n_rand/n_target*w1_A2_trim[S_target==1]/sum(w1_A2_trim)*(Y_rand-Y2_rand_pred) + Y2_rand_pred # pre-subsetting to RCT data; weights then subset to RCT A=1 data
  Y2_obs_pred_CCDS_AIPW =  n_obs/n_target*w2_A2_trim[S_target==0]/sum(w2_A2_trim)*(Y_obs-Y2_obs_pred) + Y2_obs_pred # pre-subsetting to obs data
  Y2_obs_pred_bias_CCDS_AIPWa = n_obs/n_target*w3_A2_trim[S_target==0]/sum(w3_A2_trim)*(Y_obs-Y2_obs_pred_bias_CCDSa) + Y2_obs_pred_bias_CCDSa # pre-subsetting to obs data
  Y2_pred_bias_CCDS_AIPWb = ifelse(S_target==1,n_obs/n_target*w4_A1_trim[S_target==1]/sum(w4_A1_trim)*(Y_rand-Y2_rand_pred_bias_CCDSb),
                                   Y2_obs_pred_bias_CCDSb)
  # Y3
  Y3_rand_pred_CCDS_AIPW = Y3_obs_pred_CCDS_AIPW = Y3_obs_pred_bias_CCDS_AIPWa = Y3_pred_bias_CCDS_AIPWb = NA
  if(!is.na(levels_A[3])){
    Y3_rand_pred_CCDS_AIPW = n_rand/n_target*w1_A3_trim[S_target==1]/sum(w1_A3_trim)*(Y_rand-Y3_rand_pred) + Y3_rand_pred # pre-subsetting to RCT data; weights then subset to RCT A=1 data
    Y3_obs_pred_CCDS_AIPW =  n_obs/n_target*w2_A3_trim[S_target==0]/sum(w2_A3_trim)*(Y_obs-Y3_obs_pred) + Y3_obs_pred # pre-subsetting to obs data
    Y3_obs_pred_bias_CCDS_AIPWa = n_obs/n_target*w3_A3_trim[S_target==0]/sum(w3_A3_trim)*(Y_obs-Y3_obs_pred_bias_CCDSa) + Y3_obs_pred_bias_CCDSa # pre-subsetting to obs data
    Y3_pred_bias_CCDS_AIPWb = ifelse(S_target==1,n_obs/n_target*w4_A1_trim[S_target==1]/sum(w4_A1_trim)*(Y_rand-Y3_rand_pred_bias_CCDSb),
                                     Y3_obs_pred_bias_CCDSb) # combo of RCT (first component) and obs (second component) data
  }
  
  # Target population estimates: mean potential outcome estimates
  mean_Y1_target_pred_CCDS_AIPW = 1/n_target*(sum(Y1_rand_pred_CCDS_AIPW)+sum(Y1_obs_pred_CCDS_AIPW)-
                                                sum(Y1_obs_pred_bias_CCDS_AIPWa)+sum(Y1_pred_bias_CCDS_AIPWb))
  mean_Y2_target_pred_CCDS_AIPW = 1/n_target*(sum(Y2_rand_pred_CCDS_AIPW)+sum(Y2_obs_pred_CCDS_AIPW)-
                                                sum(Y2_obs_pred_bias_CCDS_AIPWa)+sum(Y2_pred_bias_CCDS_AIPWb))
  mean_Y3_target_pred_CCDS_AIPW = 1/n_target*(sum(Y3_rand_pred_CCDS_AIPW)+sum(Y3_obs_pred_CCDS_AIPW)-
                                                sum(Y3_obs_pred_bias_CCDS_AIPWa)+sum(Y3_pred_bias_CCDS_AIPWb))
  
  ##########################################################################################
  ### Estimate potential outcomes: obs/RCT
  ##########################################################################################
  ### Observational/randomized model predictions for observational/randomized data, respectively (ignoring unmeasured confounding)
  # Y1
  Y1_rand_pred_obs_rand = Y1_rand_pred # RCT estimate same as with CCDS-OR
  Y1_obs_pred_obs_rand = Y1_obs_pred  # obs estimate same as biased obs estimate in CCDS-OR
  
  # Y2
  Y2_rand_pred_obs_rand = Y2_rand_pred # RCT estimate same as with CCDS-OR
  Y2_obs_pred_obs_rand = Y2_obs_pred  # obs estimate same as biased obs estimate in CCDS-OR
  
  # Y3
  if(!is.na(levels_A[3])){
    Y3_rand_pred_obs_rand = Y3_rand_pred # RCT estimate same as with CCDS-OR
    Y3_obs_pred_obs_rand = Y3_obs_pred  # obs estimate same as biased obs estimate in CCDS-OR
  }
  
  # Target population estimates
  Y1_target_pred_obs_rand = rep(NA,n_target)
  Y1_target_pred_obs_rand[S_target==1] = Y1_rand_pred_obs_rand
  Y1_target_pred_obs_rand[S_target==0] = Y1_obs_pred_obs_rand
  
  Y2_target_pred_obs_rand = rep(NA,n_target)
  Y2_target_pred_obs_rand[S_target==1] = Y2_rand_pred_obs_rand
  Y2_target_pred_obs_rand[S_target==0] = Y2_obs_pred_obs_rand
  
  Y3_target_pred_obs_rand = rep(NA,n_target)
  if(!is.na(levels_A[3])) {
    Y3_target_pred_obs_rand[S_target==1] = Y3_rand_pred_obs_rand
    Y3_target_pred_obs_rand[S_target==0] = Y3_obs_pred_obs_rand
  }
  
  # Mean potential outcome estimates
  mean_Y1_target_pred_obs_rand = mean(Y1_target_pred_obs_rand)
  mean_Y2_target_pred_obs_rand = mean(Y2_target_pred_obs_rand)
  mean_Y3_target_pred_obs_rand = mean(Y3_target_pred_obs_rand)
  
  
  ##########################################################################################
  ### Estimate potential outcomes: RCT
  ##########################################################################################
  ### Randomized model predictions for all data (ignoring positivity violations)
  # Target population estimates
  if(complex_fit_models=="ensemble"){
    Y1_target_pred_rand = predict(fit_rand,data_A1_target_matrix, onlySL = TRUE)$pred
    Y2_target_pred_rand = predict(fit_rand,data_A2_target_matrix, onlySL = TRUE)$pred
    if(!is.na(levels_A[3])){
      Y3_target_pred_rand = predict(fit_rand,data_A3_target_matrix, onlySL = TRUE)$pred
    }
  } else if(complex_fit_models=="ranger"){
    Y1_target_pred_rand = predict(fit_rand,data_A1_target)$predictions
    Y2_target_pred_rand = predict(fit_rand,data_A2_target)$predictions
    if(!is.na(levels_A[3])){
      Y3_target_pred_rand = predict(fit_rand,data_A3_target)$predictions
    }
  } else {
    Y1_target_pred_rand = predict(fit_rand,data_A1_target_matrix)
    Y2_target_pred_rand = predict(fit_rand,data_A2_target_matrix)
    if(!is.na(levels_A[3])){
      Y3_target_pred_rand = predict(fit_rand,data_A3_target_matrix)
    }
    
  }
  
  # Mean potential outcome estimates
  mean_Y1_target_pred_rand = mean(Y1_target_pred_rand)
  mean_Y2_target_pred_rand = mean(Y2_target_pred_rand)
  if(!is.na(levels_A[3])){
    mean_Y3_target_pred_rand = mean(Y3_target_pred_rand)
  } else {
    mean_Y3_target_pred_rand = NA
  }
  
  
  
  ##########################################################################################
  ### Estimate potential outcomes: 2-stage whole data (Kallus-like)
  ##########################################################################################
  ### Key differences from 2-stage WD: they standardize to X_rand to predict bias vs we standardize to X_overlap
  # (1) Use S=0 data to create CATE (Y_a for us) model, fit_obs
  # (2) Predict CATE (Y_a for us) for S=1 overlap data using model from (1)
  if(complex_fit_models=="ensemble"){
    Y1_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A1_rand_matrix, onlySL = TRUE)$pred
    Y2_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A2_rand_matrix, onlySL = TRUE)$pred
    if(!is.na(levels_A[3])){
      Y3_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A3_rand_matrix, onlySL = TRUE)$pred
    }
    
  } else if(complex_fit_models=="ranger"){
    Y1_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A1_rand)$predictions
    Y2_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A2_rand)$predictions
    if(!is.na(levels_A[3])){
      Y3_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A3_rand)$predictions
    }
    
  } else {
    Y1_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A1_rand_matrix)
    Y2_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A2_rand_matrix)
    if(!is.na(levels_A[3])){
      Y3_rand_pred_from_fit_obs_2stageWD = predict(fit_obs,data_A3_rand_matrix)
    }
    
  }
  
  # (3) Use S=1 data to estimate CATE (Y_a for us) for S=1 (we model this relationship)
  Y1_rand_pred_2stageWD = Y1_rand_pred
  Y2_rand_pred_2stageWD = Y2_rand_pred
  
  
  # (4) Estimate bias: eta_pred = (3) - (2)
  Y1_eta_pred_2stageWD = Y1_rand_pred_2stageWD - Y1_rand_pred_from_fit_obs_2stageWD
  Y2_eta_pred_2stageWD = Y2_rand_pred_2stageWD - Y2_rand_pred_from_fit_obs_2stageWD
  
  # (5) Create linear estimator for eta_hat: fit linear model
  fit_Y1_eta_pred_2stageWD = cbind.data.frame(Y1_eta_pred_2stageWD,X1to4_rand_matrix) %>% lm(Y1_eta_pred_2stageWD ~ ., data=.) 
  fit_Y2_eta_pred_2stageWD = cbind.data.frame(Y2_eta_pred_2stageWD,X1to4_rand_matrix) %>% lm(Y2_eta_pred_2stageWD ~ ., data=.) 
  
  # Repeat for Y3 if it exists
  if(!is.na(levels_A[3])){
    Y3_rand_pred_2stageWD = Y3_rand_pred
    Y3_eta_pred_2stageWD = Y3_rand_pred_2stageWD - Y3_rand_pred_from_fit_obs_2stageWD
    fit_Y3_eta_pred_2stageWD = cbind.data.frame(Y3_eta_pred_2stageWD,X1to4_rand_matrix) %>% lm(Y3_eta_pred_2stageWD ~ ., data=.) 
  }
  
  # Weighted version
  ### Weight for 2nd stage bias term
  w_bias_wd = ifelse(S_target==1,(1-pi_S)/(pi_S),0)
  
  # Trim weights
  point_trim_w_bias_wd =  mean(w_bias_wd)+5*sqrt(var(w_bias_wd))
  n_trim_w_bias_wd = mean(w_bias_wd>point_trim_w_bias_wd)
  w_bias_wd_trim = ifelse(w_bias_wd>point_trim_w_bias_wd,max(w_bias_wd[w_bias_wd<point_trim_w_bias_wd]),w_bias_wd)
  
  # Normalize weights
  w_bias_wd_norm = w_bias_wd_trim/sum(w_bias_wd_trim)
  w_bias_wd_norm_rand = w_bias_wd_norm[S_target==1]
 
  # Weighted regression to estimate bias
  fit_Y1_eta_pred_weighted2stageWD = cbind.data.frame(Y1_eta_pred_2stageWD,X1to4_rand_matrix) %>% 
    lm(Y1_eta_pred_2stageWD ~ ., data=., weights=w_bias_wd_norm_rand) 
  fit_Y2_eta_pred_weighted2stageWD = cbind.data.frame(Y2_eta_pred_2stageWD,X1to4_rand_matrix) %>% 
    lm(Y2_eta_pred_2stageWD ~ ., data=., weights=w_bias_wd_norm_rand) 
  
  if(!is.na(levels_A[3])){
    Y3_rand_pred_weighted2stageWD = Y3_rand_pred
    Y3_eta_pred_weighted2stageWD = Y3_rand_pred_weighted2stageWD - Y3_rand_pred_from_fit_obs_weighted2stageWD
    fit_Y3_eta_pred_weighted2stageWD = cbind.data.frame(Y3_eta_pred_weighted2stageWD,X1to4_rand_matrix) %>% 
      lm(Y3_eta_pred_weighted2stageWD ~ ., data=., weights=w_bias_wd_norm_rand) 
  }
  
  # (6) Debias CATE (Y_a for us) prediction for S=0: use (1) to predict for X_obs + eta_hat(X_obs)
  Y1_obs_pred_2stageWD = Y1_obs_pred # biased observational data prediction
  Y1_obs_pred_bias_2stageWD = predict(fit_Y1_eta_pred_2stageWD,as.data.frame(data_X1to4_obs_matrix)) # estimate of bias from linear model in (5) 
  Y1_obs_pred_debiased_2stageWD = Y1_obs_pred_2stageWD + Y1_obs_pred_bias_2stageWD # debiased observational data estimate
  
  Y2_obs_pred_2stageWD = Y2_obs_pred # biased observational data prediction
  Y2_obs_pred_bias_2stageWD = predict(fit_Y2_eta_pred_2stageWD,as.data.frame(data_X1to4_obs_matrix)) # estimate of bias from linear model in (5) 
  Y2_obs_pred_debiased_2stageWD = Y2_obs_pred_2stageWD + Y2_obs_pred_bias_2stageWD # debiased observational data estimate
  
  if(!is.na(levels_A[3])){
    Y3_obs_pred_2stageWD = Y3_obs_pred # biased observational data prediction
    Y3_obs_pred_bias_2stageWD = predict(fit_Y3_eta_pred_2stageWD,as.data.frame(data_X1to4_obs_matrix)) # estimate of bias from linear model in (5) 
    Y3_obs_pred_debiased_2stageWD = Y3_obs_pred_2stageWD + Y3_obs_pred_bias_2stageWD # debiased observational data estimate
    
  }
  
  # Weighted version
  Y1_obs_pred_weighted2stageWD = Y1_obs_pred # biased observational data prediction
  Y1_obs_pred_bias_weighted2stageWD = predict(fit_Y1_eta_pred_weighted2stageWD,as.data.frame(data_X1to4_obs_matrix)) # estimate of bias from linear model in (5) 
  Y1_obs_pred_debiased_weighted2stageCCDSWD = Y1_obs_pred_weighted2stageWD + Y1_obs_pred_bias_weighted2stageWD # debiased observational data estimate
  
  Y2_obs_pred_weighted2stageWD = Y2_obs_pred # biased observational data prediction
  Y2_obs_pred_bias_weighted2stageWD = predict(fit_Y2_eta_pred_weighted2stageWD,as.data.frame(data_X1to4_obs_matrix)) # estimate of bias from linear model in (5) 
  Y2_obs_pred_debiased_weighted2stageCCDSWD = Y2_obs_pred_weighted2stageWD + Y2_obs_pred_bias_weighted2stageWD # debiased observational data estimate
  
  if(!is.na(levels_A[3])){
    Y3_obs_pred_weighted2stageWD = Y3_obs_pred # biased observational data prediction
    Y3_obs_pred_bias_weighted2stageWD = predict(fit_Y3_eta_pred_weighted2stageWD,as.data.frame(data_X1to4_obs_matrix)) # estimate of bias from linear model in (5) 
    Y3_obs_pred_debiased_weighted2stageCCDSWD = Y3_obs_pred_weighted2stageWD + Y3_obs_pred_bias_weighted2stageWD # debiased observational data estimate
    
  }
  
  # (7) Create target population estimate by averaging across predictions for all the target sample data
  Y1_target_pred_2stageWD = rep(NA,n_target)
  Y1_target_pred_2stageWD[S_target==1] = Y1_rand_pred_2stageWD
  Y1_target_pred_2stageWD[S_target==0] = Y1_obs_pred_debiased_2stageWD
  
  Y2_target_pred_2stageWD = rep(NA,n_target)
  Y2_target_pred_2stageWD[S_target==1] = Y2_rand_pred_2stageWD
  Y2_target_pred_2stageWD[S_target==0] = Y2_obs_pred_debiased_2stageWD
  
  Y3_target_pred_2stageWD = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_2stageWD[S_target==1] = Y3_rand_pred_2stageWD
    Y3_target_pred_2stageWD[S_target==0] = Y3_obs_pred_debiased_2stageWD
  }
  mean_Y1_target_pred_2stageWD = mean(Y1_target_pred_2stageWD)
  mean_Y2_target_pred_2stageWD = mean(Y2_target_pred_2stageWD)
  mean_Y3_target_pred_2stageWD = mean(Y3_target_pred_2stageWD)
  
  # Weighted version
  Y1_target_pred_weighted2stageWD = rep(NA,n_target)
  Y1_target_pred_weighted2stageWD[S_target==1] = Y1_rand_pred_2stageWD
  Y1_target_pred_weighted2stageWD[S_target==0] = Y1_obs_pred_debiased_weighted2stageCCDSWD
  
  Y2_target_pred_weighted2stageWD = rep(NA,n_target)
  Y2_target_pred_weighted2stageWD[S_target==1] = Y2_rand_pred_2stageWD
  Y2_target_pred_weighted2stageWD[S_target==0] = Y2_obs_pred_debiased_weighted2stageCCDSWD
  
  Y3_target_pred_weighted2stageWD = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_weighted2stageWD[S_target==1] = Y3_rand_pred_2stageWD
    Y3_target_pred_weighted2stageWD[S_target==0] = Y3_obs_pred_debiased_weighted2stageCCDSWD
  }
  mean_Y1_target_pred_weighted2stageWD = mean(Y1_target_pred_weighted2stageWD)
  mean_Y2_target_pred_weighted2stageWD = mean(Y2_target_pred_weighted2stageWD)
  mean_Y3_target_pred_weighted2stageWD = mean(Y3_target_pred_weighted2stageWD)
  
  if(Lu_estimators){
    ##########################################################################################
    ### Estimate potential outcomes: Lu et al. 2019 estimators for PTSM and STSM
    ##########################################################################################
    lu_estimates = gen_estimates_lu2019_match_fit(datc = target_sample, 
                                                  X_target = X_target, 
                                                  X_target_matrix = X_target_matrix, 
                                                  Y_target = Y_target, 
                                                  S_target = S_target, 
                                                  A_target = A_target, 
                                                  X1_high_target = X1_high_target, 
                                                  X1_low_target = X1_low_target,
                                                  complex_fit_models = complex_fit_models,
                                                  SL_library = SL_library,
                                                  K=5)
  }
  
  
  ##########################################################################################
  # Evaluate bias in region of overlap and non-overlap
  # Y1
  mean_bias_Y1_obs = (mean(Y1_obs_pred_debiased_CCDS)-mean(Y1_obs))/-mean(Y1_obs)
  mean_bias_Y1_obs_observed = (mean(Y1_obs_pred)-mean(Y_obs[A_obs==1]))/mean(Y_obs[A_obs==1])
  mean_bias_Y1_obs_overlap = (mean(Y1_obs_pred_debiased_CCDS[overlap_obs==1])-
                                mean(Y1_obs[overlap_obs==1]))/mean(Y1_obs[overlap_obs==1])
  mean_bias_Y1_obs_nonoverlap = (mean(Y1_obs_pred_debiased_CCDS[overlap_obs==0])-
                                   mean(Y1_obs[overlap_obs==0]))/mean(Y1_obs[overlap_obs==0])
  
  mean_bias_Y1_rand = (mean(Y1_rand_pred)-mean(Y1_rand))/mean(Y1_rand)
  mean_bias_Y1_rand_observed = (mean(Y1_rand_pred)-mean(Y_rand[A_rand==1]))/mean(Y_rand[A_rand==1])
  mean_bias_Y1_rand_overlap = (mean(Y1_rand_pred[overlap_rand==1])-
                                 mean(Y1_rand[overlap_rand==1]))/mean(Y1_rand[overlap_rand==1])
  mean_bias_Y1_rand_overlap_observed = (mean(Y1_rand_pred[overlap_rand==1])-
                                          mean(Y_rand[A_rand==1 & overlap_rand==1]))/mean(Y_rand[A_rand==1 &overlap_rand==1])
  mean_bias_Y1_rand_nonoverlap = (mean(Y1_rand_pred[overlap_rand==0])-
                                    mean(Y1_rand[overlap_rand==0]))/mean(Y1_rand[overlap_rand==0])
  
  
  # Y2
  mean_bias_Y2_obs = (mean(Y2_obs_pred_debiased_CCDS)-mean(Y2_obs))/mean(Y2_obs)
  mean_bias_Y2_obs_observed = (mean(Y2_obs_pred)-mean(Y_obs[A_obs==2]))/mean(Y_obs[A_obs==2])
  mean_bias_Y2_obs_overlap = (mean(Y2_obs_pred_debiased_CCDS[overlap_obs==1])-
                                mean(Y2_obs[overlap_obs==1]))/mean(Y2_obs[overlap_obs==1])
  mean_bias_Y2_obs_nonoverlap = (mean(Y2_obs_pred_debiased_CCDS[overlap_obs==0])-
                                   mean(Y2_obs[overlap_obs==0]))/mean(Y2_obs[overlap_obs==0])
  
  mean_bias_Y2_rand = (mean(Y2_rand_pred)-mean(Y2_rand))/mean(Y2_rand)
  mean_bias_Y2_rand_observed = (mean(Y2_rand_pred)-mean(Y_rand[A_rand==2]))/mean(Y_rand[A_rand==2])
  mean_bias_Y2_rand_overlap = (mean(Y2_rand_pred[overlap_rand==1])-
                                 mean(Y2_rand[overlap_rand==1]))/mean(Y2_rand[overlap_rand==1])
  mean_bias_Y2_rand_overlap_observed = (mean(Y2_rand_pred[overlap_rand==1])-
                                          mean(Y_rand[A_rand==1 & overlap_rand==1]))/mean(Y_rand[A_rand==1 &overlap_rand==1])
  mean_bias_Y2_rand_nonoverlap = (mean(Y2_rand_pred[overlap_rand==0])-
                                    mean(Y2_rand[overlap_rand==0]))/mean(Y2_rand[overlap_rand==0])
  
  mean_percent_biases =c(mean_bias_Y1_obs=mean_bias_Y1_obs,
                         mean_bias_Y1_obs_observed=mean_bias_Y1_obs_observed,
                         mean_bias_Y1_obs_overlap=mean_bias_Y1_obs_overlap,
                         mean_bias_Y1_obs_nonoverlap=mean_bias_Y1_obs_nonoverlap,
                         mean_bias_Y1_rand=mean_bias_Y1_rand,
                         mean_bias_Y1_rand_observed=mean_bias_Y1_rand_observed,
                         mean_bias_Y1_rand_overlap=mean_bias_Y1_rand_overlap,
                         mean_bias_Y1_rand_overlap_observed=mean_bias_Y1_rand_overlap_observed,
                         mean_bias_Y1_rand_nonoverlap=mean_bias_Y1_rand_nonoverlap,
                         mean_bias_Y2_obs=mean_bias_Y2_obs,
                         mean_bias_Y2_obs_observed=mean_bias_Y2_obs_observed,
                         mean_bias_Y2_obs_overlap=mean_bias_Y2_obs_overlap,
                         mean_bias_Y2_obs_nonoverlap=mean_bias_Y2_obs_nonoverlap,
                         mean_bias_Y2_rand=mean_bias_Y2_rand,
                         mean_bias_Y2_rand_observed=mean_bias_Y2_rand_observed,
                         mean_bias_Y2_rand_overlap=mean_bias_Y2_rand_overlap,
                         mean_bias_Y2_rand_overlap_observed=mean_bias_Y2_rand_overlap_observed,
                         mean_bias_Y2_rand_nonoverlap=mean_bias_Y2_rand_nonoverlap)
  
  # Summarize number of trimmed weights
  n_trim = c(n_trim_w1_A1=n_trim_w1_A1,
             n_trim_w2_A1=n_trim_w2_A1,
             n_trim_w3_A1=n_trim_w3_A1,
             n_trim_w4_A1=n_trim_w4_A1,
             n_trim_w1_A2=n_trim_w1_A2,
             n_trim_w2_A2=n_trim_w2_A2,
             n_trim_w3_A2=n_trim_w3_A2,
             n_trim_w4_A2=n_trim_w4_A2,
             n_trim_w_bias=n_trim_w_bias,
             n_trim_w_bias_wd=n_trim_w_bias_wd)
  
  # Summarize mean trimmed weights
  weight_mean = c(w1_A1 = mean(w1_A1),
                  w2_A1 = mean(w2_A1),
                  w3_A1 = mean(w3_A1),
                  w4_A1 = mean(w4_A1),
                  w1_A2 = mean(w1_A2),
                  w2_A2 = mean(w2_A2),
                  w3_A2 = mean(w3_A2),
                  w4_A2 = mean(w4_A2),
                  w_bias = mean(w_bias),
                  w_bias_wd = mean(w_bias_wd))
  
  # Summarize algorithm weights if using ensemble
  if(complex_fit_models=="ensemble"){
    alg_weights_S = fit_S$coef
    alg_weights_Roverlap_rand = fit_Roverlap_rand$coef
    alg_weights_Roverlap_obs = fit_Roverlap_obs$coef
    alg_weights_A_rand = fit_A_rand$coef
    alg_weights_A_obs = fit_A_obs$coef
    alg_weights_A_rand_overlap = fit_A_rand_overlap$coef
    alg_weights_A_obs_overlap = fit_A_obs_overlap$coef
    alg_weights_rand = fit_rand$coef
    alg_weights_obs = fit_obs$coef
    alg_weights_rand_overlap = fit_rand_overlap$coef
    alg_weights_obs_overlap = fit_obs_overlap$coef
    
    names(alg_weights_S) = paste0("S_",names(alg_weights_S))
    names(alg_weights_Roverlap_rand) = paste0("Roverlap_rand_",names(alg_weights_Roverlap_rand))
    names(alg_weights_Roverlap_obs) = paste0("Roverlap_obs_",names(alg_weights_Roverlap_obs))
    names(alg_weights_A_rand) = paste0("A_rand_",names(alg_weights_A_rand))
    names(alg_weights_A_obs) = paste0("A_obs_",names(alg_weights_A_obs))
    names(alg_weights_A_rand_overlap) = paste0("A_rand_overlap_",names(alg_weights_A_rand_overlap))
    names(alg_weights_A_obs_overlap) = paste0("A_obs_overlap_",names(alg_weights_A_obs_overlap))
    names(alg_weights_rand) = paste0("rand_",names(alg_weights_rand))
    names(alg_weights_obs) = paste0("obs_",names(alg_weights_obs))
    names(alg_weights_rand_overlap) = paste0("rand_overlap_",names(alg_weights_rand_overlap))
    names(alg_weights_obs_overlap) = paste0("obs_overlap_",names(alg_weights_obs_overlap))
  }
  
  ### Output
  output = c(mean_Y1_target,
             mean_Y2_target,
             mean_Y3_target,
             
             mean_Y1_target_pred_CCDS,
             mean_Y2_target_pred_CCDS,
             mean_Y3_target_pred_CCDS,
             
             mean_Y1_target_pred_2stageCCDS,
             mean_Y2_target_pred_2stageCCDS,
             mean_Y3_target_pred_2stageCCDS,
             
             mean_Y1_target_pred_weighted2stageCCDS,
             mean_Y2_target_pred_weighted2stageCCDS,
             mean_Y3_target_pred_weighted2stageCCDS,
             
             mean_Y1_target_pred_2stageHCCDS,
             mean_Y2_target_pred_2stageHCCDS,
             mean_Y3_target_pred_2stageHCCDS,
             
             mean_Y1_target_pred_weighted2stageHCCDS,
             mean_Y2_target_pred_weighted2stageHCCDS,
             mean_Y3_target_pred_weighted2stageHCCDS,
             
             mean_Y1_target_pred_CCDS_IPW,
             mean_Y2_target_pred_CCDS_IPW,
             mean_Y3_target_pred_CCDS_IPW,
             
             mean_Y1_target_pred_CCDS_AIPW,
             mean_Y2_target_pred_CCDS_AIPW,
             mean_Y3_target_pred_CCDS_AIPW,
             
             mean_Y1_target_pred_obs_rand,
             mean_Y2_target_pred_obs_rand,
             mean_Y3_target_pred_obs_rand,
             
             mean_Y1_target_pred_rand,
             mean_Y2_target_pred_rand,
             mean_Y3_target_pred_rand,
             
             mean_Y1_target_pred_2stageWD,
             mean_Y2_target_pred_2stageWD,
             mean_Y3_target_pred_2stageWD,
             
             mean_Y1_target_pred_weighted2stageWD,
             mean_Y2_target_pred_weighted2stageWD,
             mean_Y3_target_pred_weighted2stageWD)
  
  if(Lu_estimators){
    output = c(output,
               lu_estimates,
               
               mean_overlap_target,
               mean_overlap_obs,
               mean_overlap_rand,
               min_p_overlap, 
               max_p_overlap,
               
               n_trim,
               weight_mean,
               
               mean_percent_biases)
    
    names(output) = c("mean_Y1_target",
                      "mean_Y2_target",
                      "mean_Y3_target",
                      
                      "mean_Y1_target_pred_CCDS",
                      "mean_Y2_target_pred_CCDS",
                      "mean_Y3_target_pred_CCDS",
                      
                      "mean_Y1_target_pred_2stageCCDS",
                      "mean_Y2_target_pred_2stageCCDS",
                      "mean_Y3_target_pred_2stageCCDS",
                      
                      "mean_Y1_target_pred_weighted2stageCCDS",
                      "mean_Y2_target_pred_weighted2stageCCDS",
                      "mean_Y3_target_pred_weighted2stageCCDS",
                      
                      "mean_Y1_target_pred_2stageHCCDS",
                      "mean_Y2_target_pred_2stageHCCDS",
                      "mean_Y3_target_pred_2stageHCCDS",
                      
                      "mean_Y1_target_pred_weighted2stageHCCDS",
                      "mean_Y2_target_pred_weighted2stageHCCDS",
                      "mean_Y3_target_pred_weighted2stageHCCDS",
                      
                      "mean_Y1_target_pred_CCDS_IPW",
                      "mean_Y2_target_pred_CCDS_IPW",
                      "mean_Y3_target_pred_CCDS_IPW",
                      
                      "mean_Y1_target_pred_CCDS_AIPW",
                      "mean_Y2_target_pred_CCDS_AIPW",
                      "mean_Y3_target_pred_CCDS_AIPW",
                      
                      "mean_Y1_target_pred_obs_rand",
                      "mean_Y2_target_pred_obs_rand",
                      "mean_Y3_target_pred_obs_rand",
                      
                      "mean_Y1_target_pred_rand",
                      "mean_Y2_target_pred_rand",
                      "mean_Y3_target_pred_rand",
                      
                      "mean_Y1_target_pred_2stageWD",
                      "mean_Y2_target_pred_2stageWD",
                      "mean_Y3_target_pred_2stageWD",
                      
                      "mean_Y1_target_pred_weighted2stageWD",
                      "mean_Y2_target_pred_weighted2stageWD",
                      "mean_Y3_target_pred_weighted2stageWD",
                      
                      names(lu_estimates),
                      
                      "mean_overlap_target",
                      "mean_overlap_obs",
                      "mean_overlap_rand",
                      "min_p_overlap", 
                      "max_p_overlap",
                      
                      names(n_trim),
                      names(weight_mean),
                      names(mean_percent_biases))
    
  } else {
    output = c(output,
               mean_overlap_target,
               mean_overlap_obs,
               mean_overlap_rand,
               min_p_overlap, 
               max_p_overlap,
               
               n_trim,
               weight_mean,
               
               mean_percent_biases)
    
    names(output) = c("mean_Y1_target",
                      "mean_Y2_target",
                      "mean_Y3_target",
                      
                      "mean_Y1_target_pred_CCDS",
                      "mean_Y2_target_pred_CCDS",
                      "mean_Y3_target_pred_CCDS",
                      
                      "mean_Y1_target_pred_2stageCCDS",
                      "mean_Y2_target_pred_2stageCCDS",
                      "mean_Y3_target_pred_2stageCCDS",
                      
                      "mean_Y1_target_pred_weighted2stageCCDS",
                      "mean_Y2_target_pred_weighted2stageCCDS",
                      "mean_Y3_target_pred_weighted2stageCCDS",
                      
                      "mean_Y1_target_pred_2stageHCCDS",
                      "mean_Y2_target_pred_2stageHCCDS",
                      "mean_Y3_target_pred_2stageHCCDS",
                      
                      "mean_Y1_target_pred_weighted2stageHCCDS",
                      "mean_Y2_target_pred_weighted2stageHCCDS",
                      "mean_Y3_target_pred_weighted2stageHCCDS",
                      
                      "mean_Y1_target_pred_CCDS_IPW",
                      "mean_Y2_target_pred_CCDS_IPW",
                      "mean_Y3_target_pred_CCDS_IPW",
                      
                      "mean_Y1_target_pred_CCDS_AIPW",
                      "mean_Y2_target_pred_CCDS_AIPW",
                      "mean_Y3_target_pred_CCDS_AIPW",
                      
                      "mean_Y1_target_pred_obs_rand",
                      "mean_Y2_target_pred_obs_rand",
                      "mean_Y3_target_pred_obs_rand",
                      
                      "mean_Y1_target_pred_rand",
                      "mean_Y2_target_pred_rand",
                      "mean_Y3_target_pred_rand",
                      
                      "mean_Y1_target_pred_2stageWD",
                      "mean_Y2_target_pred_2stageWD",
                      "mean_Y3_target_pred_2stageWD",
                      
                      "mean_Y1_target_pred_weighted2stageWD",
                      "mean_Y2_target_pred_weighted2stageWD",
                      "mean_Y3_target_pred_weighted2stageWD",
                      
                      "mean_overlap_target",
                      "mean_overlap_obs",
                      "mean_overlap_rand",
                      "min_p_overlap", 
                      "max_p_overlap",
                      
                      names(n_trim),
                      names(weight_mean),
                      names(mean_percent_biases))
  }

  
  if(complex_fit_models=="ensemble"){
    output = c(output,
               alg_weights_S,
               alg_weights_Roverlap_rand,
               alg_weights_Roverlap_obs,
               alg_weights_A_rand,
               alg_weights_A_obs,
               alg_weights_A_rand_overlap,
               alg_weights_A_obs_overlap,
               alg_weights_rand,
               alg_weights_obs,
               alg_weights_rand_overlap,
               alg_weights_obs_overlap)
  }
  
  # Remove fit objects
  rm(fit_S, fit_Roverlap_rand, fit_Roverlap_obs, fit_A_rand, fit_A_obs, fit_A_rand_overlap, 
     fit_A_obs_overlap, fit_rand, fit_obs, fit_rand_overlap, fit_obs_overlap)
  gc()
  
  return(output)
}



#####################################################################################################
#####################################################################################################
## Function to create summary tables of simulation results
#####################################################################################################
#####################################################################################################
get_summary_table = function(sim_results, # output of run_simulation() function
                             sim_data){ # target population dataset
  my_cols_Y1 = c("mean_Y1_target",
                 "mean_Y1_target_pred_rand",
                 "mean_Y1_target_pred_obs_rand",
                 "mean_Y1_target_pred_CCDS",
                 "mean_Y1_target_pred_CCDS_AIPW",
                 "mean_Y1_target_pred_CCDS_IPW",
                 #"mean_Y1_target_pred_2stageCCDS",
                 "mean_Y1_target_pred_weighted2stageCCDS",
                 #"mean_Y1_target_pred_2stageHCCDS",
                 "mean_Y1_target_pred_weighted2stageHCCDS",
                 #"mean_Y1_target_pred_2stageWD",
                 "mean_Y1_target_pred_weighted2stageWD",
                 "mean_Y1_target_pred_a1a2",
                 "mean_Y1_target_pred_a1a3",
                 "mean_Y1_target_pred_a1a2a3")
  my_cols_Y2 = c("mean_Y2_target",
                 "mean_Y2_target_pred_rand",
                 "mean_Y2_target_pred_obs_rand",
                 "mean_Y2_target_pred_CCDS",
                 "mean_Y2_target_pred_CCDS_AIPW",
                 "mean_Y2_target_pred_CCDS_IPW",
                 #"mean_Y2_target_pred_2stageCCDS",
                 "mean_Y2_target_pred_weighted2stageCCDS",
                 #"mean_Y2_target_pred_2stageHCCDS",
                 "mean_Y2_target_pred_weighted2stageHCCDS",
                 #"mean_Y2_target_pred_2stageWD",
                 "mean_Y2_target_pred_weighted2stageWD",
                 "mean_Y2_target_pred_a1a2",
                 "mean_Y2_target_pred_a1a3",
                 "mean_Y2_target_pred_a1a2a3")
  my_cols_Y3 = c("mean_Y3_target",
                 "mean_Y3_target_pred_rand",
                 "mean_Y3_target_pred_obs_rand",
                 "mean_Y3_target_pred_CCDS",
                 "mean_Y3_target_pred_CCDS_AIPW",
                 "mean_Y3_target_pred_CCDS_IPW",
                 #"mean_Y3_target_pred_2stageCCDS",
                 "mean_Y3_target_pred_weighted2stageCCDS",
                 #"mean_Y3_target_pred_2stageHCCDS",
                 "mean_Y3_target_pred_weighted2stageHCCDS",
                 #"mean_Y3_target_pred_2stageWD",
                 "mean_Y3_target_pred_weighted2stageWD")#,
                 # "mean_Y3_target_pred_a1a2",
                 # "mean_Y3_target_pred_a1a3",
                 # "mean_Y3_target_pred_a1a2a3")
  Estimator = factor(c("Oracle","RCT","Obs/RCT","CCDS","CCDS-AIPW","CCDS-IPW",#"2-stage CCDS",
                       "weighted 2-stage CCDS",
                       "weighted 2-stage hybrid CCDS",
                       #"2-stage WD", 
                       "weighted 2-stage WD",
                       "Lu-obs/rand1",
                       "Lu-RCT",
                       "Lu-obs/rand2"), #"Kallus"),
                     levels=c("Oracle","RCT","Lu-RCT",
                              "Obs/RCT",
                              "Lu-obs/rand1",
                              "Lu-obs/rand2",
                              "CCDS","CCDS-AIPW","CCDS-IPW",#"2-stage CCDS",
                              "weighted 2-stage CCDS",
                              "weighted 2-stage hybrid CCDS",
                              #"2-stage WD",
                              "weighted 2-stage WD"))#,"Kallus"))
  
  
  # Y1
  mean_Y1_means = sim_results %>% select(all_of(my_cols_Y1)) %>% summarize_all(.,mean,na.rm=T) %>% t()
  mean_Y1_vars = sim_results %>% select(all_of(my_cols_Y1)) %>% summarize_all(.,var,na.rm=T) %>% t() # estimator variance
  mean_Y1_CI_lower = sim_results %>% select(all_of(my_cols_Y1)) %>% summarize_all(.,quantile,probs=0.025,na.rm=T) %>% t()
  mean_Y1_CI_upper = sim_results %>% select(all_of(my_cols_Y1)) %>% summarize_all(.,quantile,probs=0.975,na.rm=T) %>% t()
  mean_Y1_bias = mean_Y1_means-rep(mean(sim_data$Y1),length(my_cols_Y1))
  mean_Y1_mse = mean_Y1_bias^2 + mean_Y1_vars
  
  Y1_table = data.frame(Estimator,"Mean"=mean_Y1_means,
                        "CI_lower"=mean_Y1_CI_lower,"CI_upper"=mean_Y1_CI_upper,
                        "Bias"=mean_Y1_bias,"Var"=mean_Y1_vars,"MSE"=mean_Y1_mse)
  rownames(Y1_table) = NULL
  
  # Y2
  mean_Y2_means = sim_results %>% select(all_of(my_cols_Y2)) %>% summarize_all(.,mean,na.rm=T) %>% t()
  mean_Y2_vars = sim_results %>% select(all_of(my_cols_Y2)) %>% summarize_all(.,var,na.rm=T) %>% t() # estimator variance
  mean_Y2_CI_lower = sim_results %>% select(all_of(my_cols_Y2)) %>% summarize_all(.,quantile,probs=0.025,na.rm=T) %>% t()
  mean_Y2_CI_upper = sim_results %>% select(all_of(my_cols_Y2)) %>% summarize_all(.,quantile,probs=0.975,na.rm=T) %>% t()
  mean_Y2_bias = mean_Y2_means-rep(mean(sim_data$Y2),length(my_cols_Y2))
  mean_Y2_mse = mean_Y2_bias^2 + mean_Y2_vars
  
  
  Y2_table = data.frame(Estimator,"Mean"=mean_Y2_means,
                        "CI_lower"=mean_Y2_CI_lower,"CI_upper"=mean_Y2_CI_upper,
                        "Bias"=mean_Y2_bias,"Var"=mean_Y2_vars,"MSE"=mean_Y2_mse)
  rownames(Y2_table) = NULL
  
  
  # Y1 - Y2
  Y1_means = sim_results %>% select(all_of(my_cols_Y1)) 
  Y2_means = sim_results %>% select(all_of(my_cols_Y2)) 
  diff_means = Y1_means - Y2_means
  mean_diff_means = diff_means %>% summarize_all(.,mean,na.rm=T) %>% t()
  mean_diff_vars = diff_means %>% summarize_all(.,var,na.rm=T) %>% t() # estimator variance
  mean_diff_CI_lower = diff_means %>% summarize_all(.,quantile,probs=0.025,na.rm=T) %>% t()
  mean_diff_CI_upper = diff_means %>% summarize_all(.,quantile,probs=0.975,na.rm=T) %>% t()
  mean_diff_bias = mean_diff_means-rep(mean(sim_data$Y1)-mean(sim_data$Y2),length(my_cols_Y1))
  mean_diff_mse = mean_diff_bias^2 + mean_diff_vars
  
  diff_table = data.frame(Estimator,"Mean"=mean_diff_means,
                          "CI_lower"=mean_diff_CI_lower,"CI_upper"=mean_diff_CI_upper,
                          "Bias"=mean_diff_bias,"Var"=mean_diff_vars,"MSE"=mean_diff_mse)
  rownames(diff_table) = NULL
  
  
  table = rbind(cbind("Estimand"="E(Y^1)",Y1_table),
                cbind("Estimand"="E(Y^2)",Y2_table),
                cbind("Estimand"="E(Y^1) - E(Y^2)",diff_table))
  
  return(table)
}



#####################################################################################################
#####################################################################################################
## Function to summarize Lu et al. 2019 comparison estimator confidence intervals from simulation results
#####################################################################################################
#####################################################################################################
get_lu2019_summary_table = function(sim_results, # output of run_simulation function
                                    sim_data, # target population data
                                    ensemble = T){ # whether ensemble was used to estimate 
 
  # Assess and eliminate iterations that failed
  print("Fraction of iterations that failed") # Note: most succeed when I rerun the simulation with that dataset again, but given how few iterations failed, I didn't think it worthwile. ensemble just happens to not converge for some seeds (https://stackoverflow.com/questions/15895897/line-search-fails-in-training-sim-prob-model)
  print(sum(sapply(sim_results,is.data.frame))/length(sim_results))
  cat("\n")
  
  if(ensemble %in% c(TRUE, "ksvm")){
    if(sum(sapply(sim_results,is.data.frame))/length(sim_results)==0){
      sim_results2 = data.frame(t(sim_results))
    } else {
      failed_indices = which(sapply(sim_results,is.data.frame))
      nonfailing_indices = setdiff(1:length(sim_results),failed_indices)
      
      sim_results2 = sim_results[-failed_indices]
      sim_results2 = do.call(rbind.data.frame, sim_results2) # get in right format
      colnames(sim_results2) = names(sim_results[[nonfailing_indices[1]]])
    }
  } else {
    sim_results2 = sim_results
  }
  
  ## Confidence interval width
  lu2019_results_lci = sim_results2[,grep("lci", colnames(sim_results2), value=TRUE)]
  lu2019_results_uci = sim_results2[,grep("uci", colnames(sim_results2), value=TRUE)]
  lu2019_results_uci = lu2019_results_uci[gsub("lci_","uci_",colnames(lu2019_results_lci))] # ensure columns align
  lu2019_ci_width = lu2019_results_uci - lu2019_results_lci
  
  
  ## Confidence interval coverage
  truth = apply(sim_data[,c("Y1","Y2")],2,mean)
  truth["Y1minusY2"] = truth["Y1"] - truth["Y2"]
  truth_vector = c(rep(c(truth["Y1"],truth["Y2"]),3),rep(truth["Y1minusY2"],3))
  lu2019_ci_coverage = sapply(1:length(truth_vector), function(i){
    lu2019_results_lci[,i] < truth_vector[i] & truth_vector[i] < lu2019_results_uci[,i]
  }) %>% 
    set_colnames(gsub("lci_","",colnames(lu2019_results_lci)))
 
  ## Standard error
  lu2019_results_se = sim_results2[,grep("se_", colnames(sim_results2), value=TRUE)]
  
  # Summarize simulation results
  mean_boot_se = apply(lu2019_results_se,2,mean)
  mean_ci_width = apply(lu2019_ci_width,2,mean)
  mean_coverage = apply(lu2019_ci_coverage,2,mean)
  mean_bounds_lower = apply(lu2019_results_lci,2,mean)
  mean_bounds_upper = apply(lu2019_results_uci,2,mean)
  
  boot_results = cbind.data.frame("mean_coverage"=mean_coverage,
                                  "mean_boot_se"=mean_boot_se,
                                  "mean_ci_width"=mean_ci_width,
                                  "mean_bounds_lower"=mean_bounds_lower,
                                  "mean_bounds_upper"=mean_bounds_upper)
  
  my_cols_Y1 = c("Y1_target_pred_a1a2",
                 "Y1_target_pred_a1a3",
                 "Y1_target_pred_a1a2a3")
  my_cols_Y2 = c("Y2_target_pred_a1a2",
                 "Y2_target_pred_a1a3",
                 "Y2_target_pred_a1a2a3")
  
  my_cols_Y1minusY2 = c("Y1minusY2_target_pred_a1a2",
                        "Y1minusY2_target_pred_a1a3",
                        "Y1minusY2_target_pred_a1a2a3")
  
  Estimator = factor(c("Lu-obs/rand1",
                       "Lu-RCT",
                       "Lu-obs/rand2"), 
                     levels=c("Lu-RCT",
                              "Lu-obs/rand1",
                              "Lu-obs/rand2"))
  
  # Pull out columns for each PTSM estimate
  Y1_boot_results = boot_results[my_cols_Y1,]  
  Y2_boot_results = boot_results[my_cols_Y2,]
  Y1minusY2_boot_results = boot_results[my_cols_Y1minusY2,]
  
  # Reorganize
  Y1_long = cbind(Estimand="E(Y^1)",
                  Estimator=Estimator,
                  Y1_boot_results)
  Y2_long = cbind(Estimand="E(Y^2)",
                  Estimator=Estimator,
                  Y2_boot_results)
  Y1minusY2_long = cbind(Estimand="E(Y^1) - E(Y^2)",
                         Estimator=Estimator,
                         Y1minusY2_boot_results)
  
  out = rbind(Y1_long, Y2_long, Y1minusY2_long) %>% 
    mutate(Estimator = fct_relevel(Estimator,c("Lu-RCT",
                                               "Lu-obs/rand1",
                                               "Lu-obs/rand2")),
           Estimand = fct_relevel(Estimand,c("E(Y^1)",
                                               "E(Y^2)",
                                               "E(Y^1) - E(Y^2)"))) %>% 
    arrange(Estimand, Estimator)
  
  return(out)
}




#####################################################################################################
#####################################################################################################
## Function to create plot of absolute bias and RMSE across estimators for only weighted estimators
#####################################################################################################
#####################################################################################################
plot_bias_rmse_results = function(my_tables, # stack of outputs of get_summary_table() function with "Type" column corresponding to names
                                  my_estimands = c("E(Y^1)","E(Y^2)","E(Y^1) - E(Y^2)"), # estimands to summarize in plots, from c("E(Y^1)","E(Y^2)","E(Y^3)","E(Y^1) - E(Y^2)")
                                  exclude_estimators = c("Oracle", "2-stage CCDS", "2-stage hybrid CCDS", "weighted 2-stage hybrid CCDS", 
                                                         "2-stage WD", "weighted 2-stage WD",
                                                         "Lu-RCT", "Lu-obs/rand1", "Lu-obs/rand2"), # Estimators in my_tables to exclude from plot
                                  my_ylim=c(0,max(sqrt(my_tables$MSE), na.rm=T))){ # y-axis limits
  
  # Format table
  my_tables2 = my_tables[my_tables$Estimand %in% my_estimands & 
                         my_tables$Estimator %notin%  exclude_estimators,] %>% 
    mutate(Estimator = gsub("weighted","",Estimator),
           Estimator = recode_factor(Estimator, "CCDS"="CCDS-OR"),
           Estimator = factor(Estimator, levels = c("RCT","Lu-RCT",
                                                    "Obs/RCT",
                                                    "Lu-obs/rand1",
                                                    "Lu-obs/rand2",
                                                    " 2-stage WD",
                                                    "CCDS-OR","CCDS-AIPW"," 2-stage CCDS","CCDS-IPW")),
           Abs_Bias = abs(Bias),
           RMSE = sqrt(MSE),
           RMSE_minus_Abs_Bias = RMSE - Abs_Bias)
  
  # Ensure facet_grid labels include exponents and line breaks
  my_estimands = factor(my_estimands, levels = my_estimands, labels = my_estimands)
  my_tables2$Type = factor(my_tables2$Type, labels = c(~atop("Correctly specified", "true overlap"),
                                                       ~atop("Main terms", "true overlap"),
                                                       ~atop("Ensemble", "true overlap"),
                                                       ~atop("Ensemble", "estimated overlap")))
  
  ############################################################################################
  ### Plots of predictions and biases (population average treatment mean, PATM, perspective)
  ### Plot of population mean estimates
  ############################################################################################
  my_tables2 %>% 
    select(Type, Estimand, Estimator, Abs_Bias, RMSE_minus_Abs_Bias) %>% 
    gather(Statistic, Amount, Abs_Bias:RMSE_minus_Abs_Bias) %>% 
    ggplot(., aes(x=Estimator, y=Amount, fill=forcats::fct_rev(Statistic))) +
    geom_bar(stat="identity", show.legend = FALSE) +
    theme_bw() + #base_size = 35
    xlab("") + ylab("Bias/RMSE") + ylim(my_ylim) +
    facet_grid(factor(Estimand,levels=my_estimands) ~ factor(Type), scales = "free",labeller = "label_parsed") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.caption = element_text(hjust = 0, face= "italic"))+
    scale_fill_manual(values=(c("#b2df8a","#1f78b4"))) 
  
}


#####################################################################################################
#####################################################################################################
## Function to create plot of coverage across estimators for only weighted estimators
#####################################################################################################
#####################################################################################################
plot_coverage_results = function(my_tables, # stack of outputs of get_summary_table() function with "Type" column corresponding to names
                                  my_estimands = c("E(Y^1)","E(Y^2)","E(Y^1) - E(Y^2)"), # estimands to summarize in plots, from c("E(Y^1)","E(Y^2)","E(Y^3)","E(Y^1) - E(Y^2)")
                                 exclude_estimators = c("Oracle", "2-stage CCDS", "2-stage hybrid CCDS", "weighted 2-stage hybrid CCDS", 
                                                        "2-stage WD", "weighted 2-stage WD"), # Estimators in my_tables to exclude from plot
                                  my_ylim=c(-0.05,1.05)){ # y-axis limits
  
  # Format table
  my_tables2 = my_tables[my_tables$Estimand %in% my_estimands & 
                           my_tables$Estimator %notin% exclude_estimators, ] %>% 
    mutate(Estimator = gsub("weighted","",Estimator),
           Estimator = recode_factor(Estimator, "CCDS"="CCDS-OR"),
           Estimator = factor(Estimator,levels = c("RCT","Obs/RCT",
                                                   "Lu-RCT",
                                                   "Lu-obs/rand1",
                                                   "Lu-obs/rand2"," 2-stage WD",
                                                   "CCDS-OR","CCDS-AIPW"," 2-stage CCDS","CCDS-IPW")))
  
  # Ensure facet_grid labels include exponents and line breaks
  my_estimands = factor(my_estimands, levels = my_estimands, labels = my_estimands)
  my_tables2$Type = factor(my_tables2$Type, labels = c(~atop("Correctly specified", "true overlap"),
                                                       ~atop("Main terms", "true overlap"),
                                                       #~atop("Ensemble", "true overlap"),
                                                       ~atop("Ensemble", "estimated overlap")))
  
  ############################################################################################
  ### Plots of coverage and CI width (for estimating population average treatment mean, PATM)
  ############################################################################################
  my_tables2 %>% 
    select(Type, Estimand, Estimator, mean_coverage, mean_ci_width) %>% 
    ggplot(., aes(x=mean_ci_width, y=mean_coverage, col=Estimator, shape=Estimator)) + 
    geom_point(position = position_jitter(width = 0.02, height = 0.02)) +
    geom_hline(yintercept=0.95, linetype='dashed') + 
    theme_bw() + theme(legend.position="right", 
                       legend.title = element_blank(), 
                       legend.text = element_text(size=8),
                       plot.caption = element_text(hjust = 0, face= "italic")) + 
    xlab("Mean CI Width") + ylab("Coverage Probability") + ylim(my_ylim) +
    scale_shape_manual(values = c(16,17,15,3,7,8,4)) +
    facet_grid(factor(Estimand,levels=my_estimands) ~ factor(Type), scales = "free",labeller = "label_parsed") #+
    #ggtitle("Coverage and CI width") +
    #labs(caption = "Dashed line at 95%")
  
}

#####################################################################################################
#####################################################################################################
## Function to create plot of bias and RMSE for EY1 across estimators with only weighted estimators
#####################################################################################################
#####################################################################################################
plot_bias_rmse_results_EY1 = function(my_tables, # stack of outputs of get_summary_table function with "Type" column corresponding names
                                      my_estimands = c("E(Y^1)"), # estimands to summarize in plots, from c("E(Y^1)","E(Y^2)","E(Y^3)","E(Y^1) - E(Y^2)")
                                      exclude_estimators = c("Oracle", "2-stage CCDS", "2-stage hybrid CCDS", "weighted 2-stage hybrid CCDS", 
                                                             "2-stage WD", "weighted 2-stage WD",
                                                             "Lu-RCT", "Lu-obs/RCT1", "Lu-obs/RCT2"), # Estimators in my_tables to exclude from plot
                                      my_ylim=c(0,max(sqrt(my_tables$MSE), na.rm=T))){ # target population on which simulation was run
  my_tables2 = my_tables[my_tables$Estimand %in% my_estimands & 
                           my_tables$Estimator %notin%  exclude_estimators,] %>% 
    mutate(Estimator = gsub("weighted","",Estimator),
           Estimator = recode_factor(Estimator, "CCDS"="CCDS-OR"),
           Estimator = factor(Estimator, levels = c("RCT","Lu-RCT",
                                                    "Obs/RCT",
                                                    "Lu-obs/RCT1",
                                                    "Lu-obs/RCT2",
                                                    " 2-stage WD",
                                                    "CCDS-OR","CCDS-AIPW"," 2-stage CCDS","CCDS-IPW")),
           Abs_Bias = abs(Bias),
           RMSE = sqrt(MSE),
           RMSE_minus_Abs_Bias = RMSE - Abs_Bias)
  
  # Ensure facet_grid labels include exponents and line breaks
  my_estimands = factor(my_estimands, levels = my_estimands, labels = my_estimands)
  my_tables2$Type = factor(my_tables2$Type, labels = c(~atop("Correctly specified", "true overlap"),
                                                       ~atop("Main terms", "true overlap"),
                                                       ~atop("Ensemble", "true overlap"),
                                                       ~atop("Ensemble", "estimated overlap")))
  ############################################################################################
  ### Plots of predictions and biases (population average treatment mean, PATM, perspective)
  ### Plot of population mean estimates
  ############################################################################################
  my_tables2 %>% 
    select(Type, Estimand, Estimator, Abs_Bias, RMSE_minus_Abs_Bias) %>% 
    gather(Statistic, Amount, Abs_Bias:RMSE_minus_Abs_Bias) %>% 
    ggplot(., aes(x=Estimator, y=Amount, fill=forcats::fct_rev(Statistic))) +
    geom_bar(stat="identity", show.legend = FALSE) +
    theme_bw() + #base_size = 35
    xlab("") + ylab("Bias/RMSE") + ylim(my_ylim) +
    facet_grid(. ~ factor(Type), scales = "free",labeller = "label_parsed") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.caption = element_text(hjust = 0, face= "italic"))+#,
    scale_fill_manual(values=(c("#b2df8a","#1f78b4"))) 
  
}


#####################################################################################################
#####################################################################################################
## Function to create plot of coverage for EY1 across estimators with only weighted estimators
#####################################################################################################
#####################################################################################################
plot_coverage_results_EY1 = function(my_tables, # stack of outputs of get_summary_table function with "Type" column corresponding names
                                     my_estimands = c("E(Y^1)"), # estimands to summarize in plots, from c("E(Y^1)","E(Y^2)","E(Y^3)","E(Y^1) - E(Y^2)")
                                     exclude_estimators = c("Oracle", "2-stage CCDS", "2-stage hybrid CCDS", "weighted 2-stage hybrid CCDS", 
                                                            "2-stage WD", "weighted 2-stage WD",
                                                            "Lu-RCT", "Lu-obs/RCT1", "Lu-obs/RCT2"), # Estimators in my_tables to exclude from plot
                                     my_ylim=c(-0.05,1.05)){ # target population on which simulation was run
  my_tables2 = my_tables[my_tables$Estimand %in% my_estimands & 
                           my_tables$Estimator %notin% exclude_estimators, ] %>% 
    mutate(Estimator = gsub("weighted","",Estimator),
           Estimator = recode_factor(Estimator, "CCDS"="CCDS-OR"),
           Estimator = factor(Estimator,levels = c("RCT","Obs/RCT",
                                                   "Lu-RCT",
                                                   "Lu-obs/RCT1",
                                                   "Lu-obs/RCT2"," 2-stage WD",
                                                   "CCDS-OR","CCDS-AIPW"," 2-stage CCDS","CCDS-IPW")))
  
  # Ensure facet_grid labels include exponents and line breaks
  my_estimands = factor(my_estimands, levels = my_estimands, labels = my_estimands)
  my_tables2$Type = factor(my_tables2$Type, labels = c(~atop("Correctly specified", "true overlap"),
                                                       ~atop("Main terms", "true overlap"),
                                                       ~atop("Ensemble", "estimated overlap")))
  
  ############################################################################################
  ### Plots of coverage and CI width (for estimating population average treatment mean, PATM)
  ############################################################################################
  my_tables2 %>% 
    select(Type, Estimand, Estimator, mean_coverage, mean_ci_width) %>% 
    #gather(Statistic, Amount, Abs_Bias:RMSE_minus_Abs_Bias) %>% 
    ggplot(., aes(x=mean_ci_width, y=mean_coverage, col=Estimator, shape=Estimator)) + #, fill=forcats::fct_rev(Statistic))) +
    geom_point(position = position_jitter(width = 0.02, height = 0.02)) +
    geom_hline(yintercept=0.95, linetype='dashed') + 
    theme_bw() + theme(legend.position="right", 
                       legend.title = element_blank(), 
                       legend.text = element_text(size=8),
                       plot.caption = element_text(hjust = 0, face= "italic")) + #base_size = 35
    xlab("Mean CI Width") + ylab("Coverage Probability") + ylim(my_ylim) +
    scale_shape_manual(values = c(16,17,15,3,7,8,4,0,1,9)) +
    facet_grid(. ~ factor(Type), scales = "free",labeller = "label_parsed") #+
  
  
}

#####################################################################################################
#####################################################################################################
## Function to create line plots of bias and RMSE results across estimators
#####################################################################################################
#####################################################################################################
plot_bias_rmse_separately = function(my_tables, # stack of outputs of get_summary_table() function with "Type" column corresponding to names
                                     my_estimands = c("E(Y^1) - E(Y^2)"), # estimands to summarize in plots, from c("E(Y^1)","E(Y^2)","E(Y^3)","E(Y^1) - E(Y^2)")
                                     estimators_remove = c("Oracle", "2-stage CCDS", "2-stage hybrid CCDS", 
                                                           "weighted 2-stage hybrid CCDS",
                                                           "2-stage WD", "weighted 2-stage WD",
                                                           "Lu-RCT", "Lu-obs/rand1", "Lu-obs/rand2")){ # Estimators in my_tables to exclude from plot
  # Format table
  my_tables2 = my_tables[my_tables$Estimand %in% my_estimands & 
                           my_tables$Estimator %notin% estimators_remove ,] %>% 
    mutate(Estimator = gsub("weighted","",Estimator),
           Estimator = recode_factor(Estimator, "CCDS"="CCDS-OR"),
           Estimator = factor(Estimator,levels = c("RCT","Obs/RCT",
                                                   "Lu-RCT",
                                                   "Lu-obs/rand1",
                                                   "Lu-obs/rand2"," 2-stage WD",
                                                   "CCDS-OR","CCDS-AIPW"," 2-stage CCDS","CCDS-IPW")), 
           Abs_Bias = abs(Bias),
           RMSE = sqrt(MSE),
           RMSE_minus_Abs_Bias = RMSE - Abs_Bias)
  
  # Plot
  my_ylim=c(0,max(sqrt(my_tables2$MSE), na.rm=T))
  
  my_tables2 %>% 
    select(Type, Setting, Estimand, Estimator, Abs_Bias, RMSE) %>% 
    rename(`Absolute bias` = Abs_Bias) %>% 
    gather(Statistic, Amount, `Absolute bias`:RMSE) %>% 
    ggplot(., aes(x=Setting, y=Amount, col=Estimator, group = Estimator)) +
    geom_line() +
    geom_point() +
    theme_bw() + theme(legend.position="right", 
                       legend.title = element_blank(), 
                       legend.text = element_text(size=8)) + #base_size = 35
    xlab("") + 
    ylab("") + ylim(my_ylim) + 
    facet_grid(factor(Statistic) ~ factor(Type), scales = "free") #+
  # ggtitle("Absolute bias and RMSE")
}

#####################################################################################################
#####################################################################################################
## Function to extract legend from ggplot, from http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
#####################################################################################################
#####################################################################################################
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#####################################################################################################
#####################################################################################################
## Function to summarize bootstrap simulation results and produce all output
#####################################################################################################
#####################################################################################################
summarize_simulation_bootstrap_results = function(sim_results){ # bootstrap simulation results from run_simulation_with_bootstrap() function
  # Transform lists to data.frames
  boot_se = data.frame(t(sapply(sim_results,function(x) x$boot_se)))
  ci_width = data.frame(t(sapply(sim_results,function(x) x$ci_width)))
  coverage = data.frame(t(sapply(sim_results,function(x) x$coverage)))
  bounds_lower = data.frame(t(sapply(sim_results,function(x) x$bounds[1,])))
  bounds_upper = data.frame(t(sapply(sim_results,function(x) x$bounds[2,])))
  
  # Summarize simulation results
  mean_boot_se = apply(boot_se,2,mean)
  mean_ci_width = apply(ci_width,2,mean)
  mean_coverage = apply(coverage,2,mean)
  mean_bounds_lower = apply(bounds_lower,2,mean)
  mean_bounds_upper = apply(bounds_upper,2,mean)
  
  boot_results = cbind.data.frame("mean_coverage"=mean_coverage,
                         "mean_boot_se"=mean_boot_se,
                         "mean_ci_width"=mean_ci_width,
                         "mean_bounds_lower"=mean_bounds_lower,
                         "mean_bounds_upper"=mean_bounds_upper)
  
  # Remove Y3 results
  boot_results = boot_results[!grepl('Y3', rownames(boot_results)), ]
  
  # Results of interest
  my_cols_Y1 = c("mean_Y1_target",
                 "mean_Y1_target_pred_rand",
                 "mean_Y1_target_pred_obs_rand",
                 "mean_Y1_target_pred_CCDS",
                 "mean_Y1_target_pred_CCDS_IPW",
                 "mean_Y1_target_pred_CCDS_AIPW",
                 #"mean_Y1_target_pred_2stageCCDS",
                 "mean_Y1_target_pred_weighted2stageCCDS",
                 #"mean_Y1_target_pred_2stageHCCDS",
                 "mean_Y1_target_pred_weighted2stageHCCDS",
                 #"mean_Y1_target_pred_2stageWD",
                 "mean_Y1_target_pred_weighted2stageWD",
                 "mean_Y1_target_pred_a1a2",
                 "mean_Y1_target_pred_a1a3",
                 "mean_Y1_target_pred_a1a2a3")
  my_cols_Y2 = c("mean_Y2_target",
                 "mean_Y2_target_pred_rand",
                 "mean_Y2_target_pred_obs_rand",
                 "mean_Y2_target_pred_CCDS",
                 "mean_Y2_target_pred_CCDS_IPW",
                 "mean_Y2_target_pred_CCDS_AIPW",
                 #"mean_Y2_target_pred_2stageCCDS",
                 "mean_Y2_target_pred_weighted2stageCCDS",
                 #"mean_Y2_target_pred_2stageHCCDS",
                 "mean_Y2_target_pred_weighted2stageHCCDS",
                 #"mean_Y2_target_pred_2stageWD",
                 "mean_Y2_target_pred_weighted2stageWD",
                 "mean_Y1_target_pred_a1a2",
                 "mean_Y1_target_pred_a1a3",
                 "mean_Y1_target_pred_a1a2a3")
  
  my_cols_Y1minusY2 = c("mean_Y1minusY2_target",
                 "mean_Y1minusY2_target_pred_rand",
                 "mean_Y1minusY2_target_pred_obs_rand",
                 "mean_Y1minusY2_target_pred_CCDS",
                 "mean_Y1minusY2_target_pred_CCDS_IPW",
                 "mean_Y1minusY2_target_pred_CCDS_AIPW",
                 #"mean_Y1minusY2_target_pred_2stageCCDS",
                 "mean_Y1minusY2_target_pred_weighted2stageCCDS",
                 #"mean_Y1minusY2_target_pred_2stageHCCDS",
                 "mean_Y1minusY2_target_pred_weighted2stageHCCDS",
                 #"mean_Y1minusY2_target_pred_2stageWD",
                 "mean_Y1minusY2_target_pred_weighted2stageWD")
  
  Estimator = factor(c("Oracle","RCT","Obs/RCT","CCDS","CCDS-IPW","CCDS-AIPW",#"2-stage CCDS",
                       "weighted 2-stage CCDS",
                       "weighted 2-stage hybrid CCDS",
                       #"2-stage WD", 
                       "weighted 2-stage WD",
                       "Lu-obs/rand1",
                       "Lu-RCT",
                       "Lu-obs/rand2"), #"Kallus"),
                     levels=c("Oracle","RCT","Obs/RCT",
                              "Lu-RCT",
                              "Lu-obs/rand1",
                              "Lu-obs/rand2",
                              "CCDS","CCDS-IPW","CCDS-AIPW",#"2-stage CCDS",
                              "weighted 2-stage CCDS",
                              "weighted 2-stage hybrid CCDS",
                              #"2-stage WD",
                              "weighted 2-stage WD"))#,"Kallus"))
  
  # Pull out columns for each PTSM estimate
  Y1_boot_results = boot_results[my_cols_Y1,]  
  Y2_boot_results = boot_results[my_cols_Y2,]
  Y1minusY2_boot_results = boot_results[my_cols_Y1minusY2,]
  
  # Reorganize
  Y1_long = cbind(Estimand="E(Y^1)",
                  Estimator=Estimator,
                  Y1_boot_results)
  Y2_long = cbind(Estimand="E(Y^2)",
                  Estimator=Estimator,
                  Y2_boot_results)
  Y1minusY2_long = cbind(Estimand="E(Y^1) - E(Y^2)",
                  Estimator=Estimator,
                  Y1minusY2_boot_results)
  
  out = rbind(Y1_long, Y2_long, Y1minusY2_long) %>% 
    mutate(Estimator = fct_relevel(Estimator,c("Oracle","RCT","Obs/RCT",
                                               "Lu-RCT",
                                               "Lu-obs/rand1",
                                               "Lu-obs/rand2",
                                               "CCDS","CCDS-IPW","CCDS-AIPW",#"2-stage CCDS",
                                               "weighted 2-stage CCDS",
                                               #"2-stage WD",
                                               "weighted 2-stage WD"))) 
  
  return(out)
}


#####################################################################################################
#####################################################################################################
## Function to get table for simulation results to feed into other results summaries -- NOTE: CAN SUBSTITUTE INTO FUNCTIONS ABOVE TO REMOVE DUPLICITY
#####################################################################################################
#####################################################################################################
get_table_simulation_results = function(sim_results, sim_data = sim_data_complex_noX1cube, ensemble=T,
                                        n_estimator = n_estimators){
  # Assess and eliminate iterations that failed
  print("Fraction of iterations that failed") # Note: most succeed when I rerun the simulation with that dataset again, but given how few iterations failed, I didn't think it worthwile. ensemble just happens to not converge for some seeds (https://stackoverflow.com/questions/15895897/line-search-fails-in-training-sim-prob-model)
  print(sum(sapply(sim_results,is.data.frame))/length(sim_results))
  cat("\n")
  
  if(ensemble %in% c(TRUE, "ksvm")){
    if(sum(sapply(sim_results,is.data.frame))/length(sim_results)==0){
      sim_results2 = data.frame(t(sim_results))
    } else {
      failed_indices = which(sapply(sim_results,is.data.frame))
      nonfailing_indices = setdiff(1:length(sim_results),failed_indices)
      
      sim_results2 = sim_results[-failed_indices]
      sim_results2 = do.call(rbind.data.frame, sim_results2) # get in right format
      colnames(sim_results2) = names(sim_results[[nonfailing_indices[1]]])
    }
  } else {
    sim_results2 = sim_results
  }
  
  # Summarize simulation results
  table_sim = get_summary_table(sim_results=sim_results2, sim_data = sim_data)
  
  return(table_sim)
}

#####################################################################################################
#####################################################################################################
## Function to summarize simulated data created by gen_conf_sim_data() function
#####################################################################################################
#####################################################################################################
# Note that some parts are duplicative to the get_estimates() function
summarize_sim_data = function(N=my_N, # target population size
                              sim_data, # simulated data
                              p_S_pop, # propensity scores for selection
                              p_A_pop, # propensity scores for treatment
                              n_target = 10000, # target sample size
                              overlap_function = "ifelse(X_target[,1]>qnorm(0.9),0,
                                  ifelse(X_target[,1]<qnorm(0.5),0,1))", # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                              min_X1_overlap=qnorm(0.5), # true X1 lower limit for overlap region
                              max_X1_overlap=qnorm(0.9), # true X1 upper limit for overlap region
                              beta_YA = -3, # strength of relationship between A and Y
                              beta_YU = my_beta_YU, # strength of relationship between U and Y
                              beta_AX = my_beta_AX, # coefficients for relationship between X and A
                              beta_YX = my_beta_YX, # coefficients for relationship between X and Y
                              beta_YAX = my_beta_YAX, # coefficient for interaction between X and A in Y regression
                              complex_gen_model = T, # whether to include complex model terms (X1^3) in the function for S, A, and Y
                              data_gen_scenario = F, # whether to change the true model specification to account for assumption violations
                              my_seed = random_seed){ # random seed
  ##########################################################################################
  # Target population: extract data
  ##########################################################################################
  target_pop = sim_data # population from which we'll be drawing samples
  X_pop = subset(target_pop,select=c(-Y,-S,-A,-U,-Y1,-Y2,-Y3))#sim_data_unconf %>% dplyr::select(-Y,-S,-A) %>% 
  X_pop_matrix = model.matrix(~.,X_pop)[,-1]
  Y_pop = target_pop %>% pull(Y)
  S_pop = target_pop %>% pull(S)
  A_pop = target_pop %>% pull(A)
  U_pop = target_pop %>% pull(U)
  Y1_pop = target_pop %>% pull(Y1)
  Y2_pop = target_pop %>% pull(Y2)
  Y3_pop = target_pop %>% pull(Y3)
  Y1minusY2_pop = Y1_pop - Y2_pop 
  
  # Target population: number of randomized individuals; percent randomized
  N_rand = sum(as.numeric(as.character(S_pop)))
  N_obs = N - N_rand
  p_S = mean(as.numeric(as.character(S_pop))) %>% round(3) # P(S=1)
  
  ##########################################################################################
  # Target sample: extract data (one instance for data summary purposes)
  ##########################################################################################
  set.seed(my_seed)
  i = sample(1:N, n_target, replace=TRUE)
  target_sample = target_pop[i,]
  X_target = subset(target_sample,select=c(-Y,-S,-A,-U,-Y1,-Y2,-Y3))
  X_target_matrix = model.matrix(~.,X_target)[,-1]
  Y_target = target_sample %>% pull(Y)
  S_target = target_sample %>% pull(S)
  A_target = target_sample %>% pull(A)
  levels_A = unique(A_target) %>% sort()
  U_target = target_sample %>% pull(U)
  
  # Counterfactual outcomes
  Y1_target = target_sample %>% pull(Y1)
  Y2_target = target_sample %>% pull(Y2)
  Y3_target = target_sample %>% pull(Y3)
  Y1minusY2_target = Y1_target - Y2_target  
  
  # Subsetting randomized and observational data
  Y_rand = Y_target[S_target==1]
  Y1_rand = Y1_target[S_target==1]
  Y2_rand = Y2_target[S_target==1]
  Y3_rand = Y3_target[S_target==1]
  A_rand = A_target[S_target==1]
  X_rand = X_target[S_target==1,]
  X_rand_matrix = model.matrix(~.,X_rand)[,-1]
  U_rand = U_target[S_target==1]
  
  Y_obs = Y_target[S_target==0]
  Y1_obs = Y1_target[S_target==0]
  Y2_obs = Y2_target[S_target==0]
  Y3_obs = Y3_target[S_target==0]
  A_obs = A_target[S_target==0]
  X_obs = X_target[S_target==0,]
  X_obs_matrix = model.matrix(~.,X_obs)[,-1]
  U_obs = U_target[S_target==0]
  
  # Target sample: number of randomized individuals; percent randomized
  n_rand = sum(as.numeric(as.character(S_target)))
  n_obs = n_target - n_rand
  p_S_1 = mean(as.numeric(as.character(S_target))) %>% round(3)
  
  # Create summary table
  table_n = rbind("Pop"=c('N_rand'=N_rand, 'N_obs'=N_obs, 'p_S'=p_S), "Sample"=c('n_rand'=n_rand, 'n_obs'=n_obs, 'p_S_1'=p_S_1))
  
  ##########################################################################################
  # Create table of potential outcomes and observed outcomes
  ##########################################################################################
  Y1_row = c(mean(Y1_pop), # Truth
             mean(Y1_pop[S_pop==1]), # RCT truth
             mean(Y_pop[S_pop==1 & A_pop==1]), # RCT pop observed Y for A=1
             mean(Y_rand[A_rand==1]), # RCT sample observed Y for A=1
             mean(Y1_pop[S_pop==0]), # Obs truth
             mean(Y_pop[S_pop==0 & A_pop==1]), # Obs pop observed Y for A=1
             mean(Y_obs[A_obs==1])) # Obs sample observed Y for A=1
  Y2_row = c(mean(Y2_pop),
             mean(Y2_pop[S_pop==1]),
             mean(Y_pop[S_pop==1 & A_pop==2]),
             mean(Y_rand[A_rand==2]),
             mean(Y2_pop[S_pop==0]),
             mean(Y_pop[S_pop==0 & A_pop==2]),
             mean(Y_obs[A_obs==2]))
  Y3_row = c(mean(Y3_pop),
             mean(Y3_pop[S_pop==1]),
             mean(Y_pop[S_pop==1 & A_pop==3]),
             mean(Y_rand[A_rand==3]),
             mean(Y3_pop[S_pop==0]),
             mean(Y_pop[S_pop==0 & A_pop==3]),
             mean(Y_obs[A_obs==3]))
  Y1minusY2_row = Y1_row - Y2_row
  
  table_outcomes = rbind(Y1_row,Y2_row,Y3_row,Y1minusY2_row)
  rownames(table_outcomes) = c("E(Y^1)","E(Y^2)","E(Y^3)","E(Y^1)-E(Y^2)")
  colnames(table_outcomes) = c("Truth", # Mean potential outcomes in target population
                               "RCT truth", "RCT pop observed Y", "RCT sample observed Y", # Check RCT internal validity, external validity bias, and sampling variability
                               "Obs truth",  "Obs pop observed Y", "Obs sample observed Y") # Check obs internal validity bias and sampling variability
  
  # Mean outcome in target population
  mean_Y_pop = mean(Y_pop) 
  
  ##########################################################################################
  # Examine propensity score for selection
  ##########################################################################################
  p_A_target = p_A_pop[i,]
  p_S_target = p_S_pop[i]
  
  # Given skewedness in the data, will examine logit propensity scores
  logit = qlogis
  logit_p_S_target = logit(p_S_target)
  
  # Plot propensity score
  # Percent
  plot_propensity_S_population = cbind.data.frame(S_pop,p_S_pop) %>% 
    mutate(S_pop = factor(S_pop, levels=c(0,1), labels=c("Observational","Randomized"))) %>% 
    ggplot(aes(p_S_pop, fill=S_pop)) +
    geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                   position='dodge') +
    ylab("Percentage") + theme_bw() + ggtitle("True propensity for S in superpopulation") +
    labs(fill="")
  
  plot_propensity_S_sample = cbind.data.frame(S_target,p_S_target) %>% 
    mutate(S_target = factor(S_target, levels=c(0,1), labels=c("Observational","Randomized"))) %>% 
    ggplot(aes(p_S_target, fill=S_target)) +
    geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                   position='dodge') +
    ylab("Percentage") + theme_bw() + ggtitle("Sample propensity for S") +
    labs(fill="")
  
  plot_propensity_S_sample_X1 = cbind.data.frame(S_target,p_S_target, "X1"=X_target[,1]) %>% 
    mutate(S_target = factor(S_target, levels=c(0,1), labels=c("Observational","Randomized"))) %>% 
    ggplot(aes(X1,p_S_target, color=S_target)) +
    geom_point(alpha=0.2) +
    theme_bw() + ggtitle("Sample propensity for S across X1") +
    labs(fill="")
  
  plot_log_propensity_S_sample_zoom = cbind.data.frame(S_target,logit_p_S_target) %>% 
    mutate(S_target = factor(S_target, levels=c(0,1), labels=c("Observational","Randomized"))) %>%
    ggplot( aes(x=logit_p_S_target, fill=S_target)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    #scale_fill_manual(values=c("#69b3a2", "#404080")) +
    theme_bw() + ggtitle("Sample logit propensity for S, zoom in") +
    coord_cartesian(ylim=c(0, 100))+
    labs(fill="")
  
  ##########################################################################################
  # Determine region of overlap
  ##########################################################################################
  # Tuning parameters for determining region of overlap
  # Given the skewdness of the data, determining overlap on logit scale
  if(overlap_function == "overlap1"){
    a = (max(p_S_target)-min(p_S_target)) %>% abs() %>% as.numeric()*.01# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 1% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    R_overlap = pw_overlap(ps=p_S_target, # logit propensity score for S
                           E=S_target, # "Exposure", for us, S
                           a=a, # tuning parameter a
                           b=b) # tuning parameter b
  } else if(overlap_function == "overlap2"){
    a = (max(p_S_target)-min(p_S_target)) %>% abs() %>% as.numeric()*.02# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 2% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    R_overlap = pw_overlap(ps=p_S_target, # logit propensity score for S
                           E=S_target, # "Exposure", for us, S
                           a=a, # tuning parameter a
                           b=b) # tuning parameter b
  } else if(overlap_function == "overlap3"){
    a = (max(logit_p_S_target)-min(logit_p_S_target)) %>% abs() %>% as.numeric()*.1# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.04*min(n_obs,n_rand),0) # less than 4% of randomized units in 10% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    R_overlap = pw_overlap(ps=logit_p_S_target, # logit propensity score for S
                           E=S_target, # "Exposure", for us, S
                           a=a, # tuning parameter a
                           b=b) # tuning parameter b
  } else if(overlap_function == "overlap4"){
    a = (max(logit_p_S_target)-min(logit_p_S_target)) %>% abs() %>% as.numeric()*.02# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 2% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    R_overlap = pw_overlap(ps=logit_p_S_target, # logit propensity score for S
                           E=S_target, # "Exposure", for us, S
                           a=a, # tuning parameter a
                           b=b) # tuning parameter b
  } else if(overlap_function==F){
    # Tuning parameters for determining region of overlap
    # Given the skewdness of the data, determining overlap on logit scale
    a = (max(logit_p_S_target)-min(logit_p_S_target)) %>% abs() %>% as.numeric()*.01# a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 1% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    R_overlap = pw_overlap(ps=logit_p_S_target, # propensity score for S
                           E=S_target, # "Exposure", for us, S
                           a=a, # tuning parameter a
                           b=b) # tuning parameter b
  } else {
    a = NA
    b = NA
    R_overlap=eval(parse(text=overlap_function))
  }
  
  
  
  # Regions of overlap
  overlap_target = R_overlap
  overlap_rand = R_overlap[S_target==1]
  overlap_obs = R_overlap[S_target==0]
  X_rand_overlap = X_rand[overlap_rand==1,]
  X_obs_overlap = X_obs[overlap_obs==1,]
  U_rand_overlap = X_rand[overlap_rand==1,]
  U_obs_overlap = X_obs[overlap_obs==1,]
  mean_overlap_target = mean(overlap_target) # proportion of data in region of overlap
  mean_overlap_rand = mean(overlap_rand) # proportion of randomized data in region of overlap
  mean_overlap_obs = mean(overlap_obs) # proportion of observational data in region of overlap
  
  # Order by propensity score and examine results
  order_p_S_target = p_S_target[order(p_S_target)]
  order_R_overlap = R_overlap[order(p_S_target)]
  order_S_target = S_target[order(S_target)]
  min_p_overlap = order_p_S_target[min(which(order_R_overlap==1))] # smallest propensity score in the region of overlap
  max_p_overlap = order_p_S_target[max(which(order_R_overlap==1))] # largest propensity score in the region of overlap
  
  # Check if overlap region is contiguous
  order_R_overlap_i = seq(min(which(order_R_overlap==1)),max(which(order_R_overlap==1)),1)
  
  # Summary table
  table_overlap = c('a'=a,
                    'b'=b,
                    'mean_overlap_target'=mean_overlap_target,
                    'mean_overlap_rand'=mean_overlap_rand,
                    'mean_overlap_obs'=mean_overlap_obs,
                    'min_p_overlap'=min_p_overlap, # Smallest propensity score in region of overlap
                    'max_p_overlap'=max_p_overlap, # Largest propensity score in region of overlap
                    'min_p_S_target'=min(order_p_S_target), # Smallest propensity score 
                    'max_p_S_target'=max(order_p_S_target), # Largest propensity score 
                    'nonoverlap_obs_in_R_overlap_span'=sum(order_R_overlap[order_R_overlap_i]==0) # No non-overlap observations within span of overlap region
  )
  
  # Plot results
  plot_propensity_S_sample_overlap = cbind.data.frame(S_target,p_S_target) %>% 
    mutate(S_target = factor(S_target, levels=c(0,1), labels=c("Observational","Randomized"))) %>% 
    ggplot(aes(p_S_target, fill=S_target)) +
    annotate("rect", xmin = min_p_overlap, xmax = max_p_overlap, ymin = 0, ymax = Inf, 
             fill="pink", color=NA, size=0.5) + 
    geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                   position='dodge') +
    ylab("Percentage") + theme_bw() + ggtitle("Sample propensity for S, overlap region in light pink") +
    labs(fill="") 
  
  plot_log_propensity_S_sample_overlap = cbind.data.frame(R_overlap,logit_p_S_target) %>%
    mutate(R_overlap = factor(R_overlap, levels=c(0,1), labels=c("Non-overlap","Overlap"))) %>%
    ggplot( aes(x=logit_p_S_target, fill=R_overlap)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', 
                    binwidth = ifelse(is.na(a),(max(logit_p_S_target)-min(logit_p_S_target)) %>% abs() %>% as.numeric()*.05,a)) +
    theme_bw() + ggtitle("Sample logit propensity for S") +
    labs(fill="")
  plot_propensity_S_sample_overlap2 = cbind.data.frame(R_overlap,p_S_target) %>% 
    mutate(R_overlap = factor(R_overlap, levels=c(0,1), labels=c("Non-overlap","Overlap"))) %>%  
    ggplot( aes(x=p_S_target, fill=R_overlap)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins=53) +
    theme_bw() + ggtitle("Sample propensity for S") +
    labs(fill="")
  plot_propensity_S_sample_overlap_zoom = cbind.data.frame(R_overlap,p_S_target) %>% mutate(R_overlap = factor(R_overlap, levels=c(0,1), labels=c("Non-overlap","Overlap"))) %>%  
    ggplot( aes(x=p_S_target, fill=R_overlap)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins=45) +
    theme_bw() + ggtitle("Sample propensity for S, zoom in") +
    xlim(0,0.05)+
    labs(fill="")
  
  
  #################################################################################################################################################
  ### EDA: U and Y distributions
  #################################################################################################################################################
  # Data distributions
  library(scales)
  plot_dist_U_by_group = cbind.data.frame(U_target,S_target) %>% 
    ggplot( aes(x=U_target, fill=as.factor(S_target))) +
    geom_bar(aes(y=..prop..), position = position_dodge(), color="#e9ecef", alpha=0.4) +
    scale_y_continuous(labels=percent_format()) +
    scale_fill_hc(name="", breaks=c("0", "1"),
                  labels=c("Observational", "Randomized")) +
    theme_bw() + ggtitle("Distribution of U by group") +
    labs(fill="")
  
  plot_dist_U_S0_by_overlap = U_obs %>% as.data.frame() %>% 
    ggplot( aes(x=U_obs, fill=as.factor(overlap_obs))) +
    geom_bar(aes(y=..prop..), position = position_dodge(), color="#e9ecef", alpha=0.4) +
    scale_y_continuous(labels=percent_format()) +
    scale_fill_ptol(name="", breaks=c("0", "1"),
                    labels=c("Non-overlap", "Overlap")) +
    theme_bw() + ggtitle("Distribution of U for S=0 by overlap status") + 
    labs(fill="")
  
  plot_dist_Y_by_group = cbind.data.frame(Y_target,S_target) %>% 
    ggplot( aes(x=Y_target, fill=as.factor(S_target))) +
    geom_density(color="#e9ecef", alpha=0.4) +    
    scale_fill_hc(name="", breaks=c("0", "1"),
                  labels=c("Observational", "Randomized")) +
    theme_bw() + ggtitle("Distribution of Y by group") +
    labs(fill="")
  
  plot_dist_Y_by_p_S = cbind.data.frame(Y_target,p_S_target,S_target) %>% 
    ggplot( aes(x=p_S_target, y=Y_target, fill=as.factor(S_target))) +
    geom_smooth(alpha=0.4) +    
    scale_fill_hc(name="", breaks=c("0", "1"),
                  labels=c("Observational", "Randomized")) +
    theme_bw() + ggtitle("Distribution of Y by P(S|X) and group") +
    xlab("P(S|X)") + ylab("Y") +
    labs(fill="")
  
  plot_dist_Ya_by_p_S = rbind(cbind("Y1",Y1_target,p_S_target,S_target),
                              cbind("Y2",Y2_target,p_S_target,S_target),
                              cbind("Y3",Y3_target,p_S_target,S_target)) %>% 
    set_colnames(.,c("Outcome","Ya","p_S","S")) %>% as.data.frame() %>% 
    mutate(Ya=as.numeric(as.character(Ya)),
           p_S=as.numeric(as.character(p_S)),
           S=as.numeric(as.character(S))) %>%  
    ggplot( aes(x=p_S, y=Ya))+#, fill=as.factor(S))) +
    geom_smooth(alpha=0.4) +    
    scale_fill_hc(name="", breaks=c("0", "1"),
                  labels=c("Observational", "Randomized")) +
    theme_bw() + ggtitle("Distribution of Y^a by P(S|X)") +
    xlab("P(S|X)") + ylab("Y^a") +
    labs(fill="") + facet_grid(Outcome~.)  
  
  plot_dist_Ya_by_p_S_group = rbind(cbind("Y1",Y1_target,p_S_target,S_target),
                                    cbind("Y2",Y2_target,p_S_target,S_target),
                                    cbind("Y3",Y3_target,p_S_target,S_target)) %>% 
    set_colnames(.,c("Outcome","Ya","p_S","S")) %>% as.data.frame() %>% 
    mutate(Ya=as.numeric(as.character(Ya)),
           p_S=as.numeric(as.character(p_S)),
           S=as.numeric(as.character(S))) %>%  
    ggplot( aes(x=p_S, y=Ya, fill=as.factor(S))) +
    geom_smooth(alpha=0.4) +    
    scale_fill_hc(name="", breaks=c("0", "1"),
                  labels=c("Observational", "Randomized")) +
    theme_bw() + ggtitle("Distribution of Y^a by P(S|X) and group") +
    xlab("P(S|X)") + ylab("Y^a") +
    labs(fill="") + facet_grid(Outcome~.) 
  
  my_sample = ifelse(n_target>1000,sample(1:n_target,1000,replace=T),seq(1:n_target))
  
  plot_dist_Ya_by_p_A_group = data.frame(rbind(cbind("Y1",Y1_target[my_sample],p_A_target[my_sample,1],S_target[my_sample]),
                                               cbind("Y2",Y2_target[my_sample],p_A_target[my_sample,2],S_target[my_sample]),
                                               cbind("Y3",Y3_target[my_sample],p_A_target[my_sample,3],S_target[my_sample]))) %>% 
    set_colnames(.,c("Outcome","Ya","p_A","S")) %>% as.data.frame() %>% 
    mutate(Ya=as.numeric(as.character(Ya)),
           p_A=as.numeric(as.character(p_A)),
           S=as.numeric(as.character(S))) %>%  
    ggplot( aes(x=p_A, y=Ya, fill=as.factor(S))) +
    geom_smooth(alpha=0.4) +    
    scale_fill_hc(name="", breaks=c("0", "1"),
                  labels=c("Observational", "Randomized")) +
    theme_bw() + ggtitle("Distribution of Y^a by P(A|X,U) and group") +
    xlab("P(A|X,U)") + ylab("Y^a") +
    labs(fill="") + facet_grid(Outcome~.) 
  
  plot_dist_Ya_by_p_A_group = data.frame(rbind(cbind("Y1",Y1_target,p_A_target[,1],S_target),
                                               cbind("Y2",Y2_target,p_A_target[,2],S_target),
                                               cbind("Y3",Y3_target,p_A_target[,3],S_target))) %>% 
    set_colnames(.,c("Outcome","Ya","p_A","S")) %>% as.data.frame() %>% 
    mutate(Ya=as.numeric(as.character(Ya)),
           p_A=as.numeric(as.character(p_A)),
           S=as.numeric(as.character(S))) %>%  
    ggplot( aes(x=p_A, y=Ya, fill=as.factor(S))) +
    geom_smooth(alpha=0.4,method = 'loess') +    
    scale_fill_hc(name="", breaks=c("0", "1"),
                  labels=c("Observational", "Randomized")) +
    theme_bw() + ggtitle("Distribution of Y^a by P(A|X,U) and group") +
    xlab("P(A|X,U)") + ylab("Y^a") +
    labs(fill="") + facet_grid(Outcome~.) 
  
  
  ##########################################################################################
  ### Assess CATE estimates from RCT and obs data
  ##########################################################################################
  ### Data in overlap region
  X_rand_overlap_matrix = model.matrix(~.,X_rand_overlap)[,-1]
  X_obs_overlap_matrix = model.matrix(~.,X_obs_overlap)[,-1]
  
  ##########################################################################################
  ### Generate fits
  # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
  fit_rand = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
    lm(Y ~ . + .:A,data=.) 
  # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
  fit_obs = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
    lm(Y ~ . + .:A,data=.) 
  
  # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
  fit_rand_overlap = cbind.data.frame("Y"=Y_rand,"A"=A_rand,X_rand_matrix) %>% 
    lm(Y ~ . + .:A,data=., subset=(overlap_rand==1)) 
  
  # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
  fit_obs_overlap = cbind.data.frame("Y"=Y_obs,"A"=A_obs,X_obs_matrix) %>% 
    lm(Y ~ . + .:A,data=., subset=(overlap_obs==1)) 
  
  
  ##########################################################################################
  ### Set up target data to estimate potential outcomes
  # Set up target data with the same variable names for prediction
  #    A2    A3     X     
  data_A1_target = cbind.data.frame("A"=levels_A[1],X_target) # target data with A=1 for all observations
  data_A2_target = cbind.data.frame("A"=levels_A[2],X_target) # target data with A=2 for all observations
  data_A3_target = cbind.data.frame("A"=levels_A[3],X_target) # target data with A=3 for all observations
  
  data_A1_target_matrix = cbind.data.frame("A"=levels_A[1],X_target_matrix) # target data with A=1 for all observations
  data_A2_target_matrix = cbind.data.frame("A"=levels_A[2],X_target_matrix) # target data with A=2 for all observations
  data_A3_target_matrix = cbind.data.frame("A"=levels_A[3],X_target_matrix) # target data with A=3 for all observations
  
  #               X      
  data_X_obs = as.data.frame(cbind(X_obs))
  data_X_obs_matrix = model.matrix(~.,data_X_obs)[,-1]
  
  
  ## A1
  # Needed for CCDS
  # Estimator name indicates data being used (e.g., obs_overlap) and model from which potential outcomes will be estimated (e.g., for_fit_rand). The latter is needed so variable names match for predicting from the model
  data_A1_rand = data_A1_target[S_target==1,]
  data_A1_obs = data_A1_target[S_target==0,]
  
  data_A1_rand_matrix = data_A1_target_matrix[S_target==1,]
  data_A1_obs_matrix = data_A1_target_matrix[S_target==0,]
  
  # Needed for 2-stage WD
  data_A1_rand_overlap = data_A1_target[S_target==1 & overlap_target==1,]
  data_A1_rand_matrix_overlap = data_A1_target_matrix[S_target==1 & overlap_target==1,]
  
  ## A2
  # Needed for CCDS
  data_A2_rand = data_A2_target[S_target==1,]
  data_A2_obs = data_A2_target[S_target==0,]
  
  data_A2_rand_matrix = data_A2_target_matrix[S_target==1,]
  data_A2_obs_matrix = data_A2_target_matrix[S_target==0,]
  
  # Needed for 2-stage WD
  data_A2_rand_overlap = data_A2_target[S_target==1 & overlap_target==1,]
  data_A2_rand_matrix_overlap = data_A2_target_matrix[S_target==1 & overlap_target==1,]
  
  ## A3
  # Needed for CCDS
  data_A3_rand = data_A3_target[S_target==1,]
  data_A3_obs = data_A3_target[S_target==0,]
  
  data_A3_rand_matrix = data_A3_target_matrix[S_target==1,]
  data_A3_obs_matrix = data_A3_target_matrix[S_target==0,]
  
  # Needed for 2-stage WD
  data_A3_rand_overlap = data_A3_target[S_target==1 & overlap_target==1,]
  data_A3_rand_matrix_overlap = data_A3_target_matrix[S_target==1 & overlap_target==1,]
  
  
  ##########################################################################################
  ### Estimate potential outcomes: novel estimators
  ### CCDS: predict counterfactuals for target sample
  # Warnings arise because of collinearity due to data generating mechanism
  # Y1
  Y1_rand_pred = predict(fit_rand,data_A1_rand_matrix)
  Y1_obs_pred = predict(fit_obs,data_A1_obs_matrix)
  Y1_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A1_obs_matrix)
  Y1_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A1_obs_matrix) 
  Y1_obs_pred_bias_CCDS = Y1_obs_pred_bias_CCDSa - Y1_obs_pred_bias_CCDSb  # bias estimate for obs data
  Y1_obs_pred_debiased_CCDS = Y1_obs_pred - # obs estimate
    Y1_obs_pred_bias_CCDS # bias estimate for obs data
  
  # Y2
  Y2_rand_pred = predict(fit_rand,data_A2_rand_matrix)
  Y2_obs_pred = predict(fit_obs,data_A2_obs_matrix)
  Y2_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A2_obs_matrix)
  Y2_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A2_obs_matrix) 
  Y2_obs_pred_bias_CCDS = Y2_obs_pred_bias_CCDSa - Y2_obs_pred_bias_CCDSb  # bias estimate for obs data
  Y2_obs_pred_debiased_CCDS = Y2_obs_pred - # obs estimate
    Y2_obs_pred_bias_CCDS # bias estimate for obs data
  
  # Y3
  if(!is.na(levels_A[3])){
    Y3_rand_pred = predict(fit_rand,data_A3_rand_matrix)
    Y3_obs_pred = predict(fit_obs,data_A3_obs_matrix)
    Y3_obs_pred_bias_CCDSa = predict(fit_obs_overlap,data_A3_obs_matrix)
    Y3_obs_pred_bias_CCDSb = predict(fit_rand_overlap,data_A3_obs_matrix) 
    Y3_obs_pred_bias_CCDS = Y3_obs_pred_bias_CCDSa - Y3_obs_pred_bias_CCDSb  # bias estimate for obs data
    Y3_obs_pred_debiased_CCDS = Y3_obs_pred - # obs estimate
      Y3_obs_pred_bias_CCDS # bias estimate for obs data
  }
  
  
  # Target population estimates
  Y1_target_pred_CCDS = rep(NA,n_target)
  Y1_target_pred_CCDS[S_target==1] = Y1_rand_pred
  Y1_target_pred_CCDS[S_target==0] = Y1_obs_pred_debiased_CCDS
  
  Y2_target_pred_CCDS = rep(NA,n_target)
  Y2_target_pred_CCDS[S_target==1] = Y2_rand_pred
  Y2_target_pred_CCDS[S_target==0] = Y2_obs_pred_debiased_CCDS
  
  Y3_target_pred_CCDS = rep(NA,n_target)
  if(!is.na(levels_A[3])){
    Y3_target_pred_CCDS[S_target==1] = Y3_rand_pred
    Y3_target_pred_CCDS[S_target==0] = Y3_obs_pred_debiased_CCDS
  }
  
  ##########################################################################################
  ### Estimate potential outcomes: comparison estimators
  ### Observational/randomized model predictions for observational/randomized data, respectively (ignoring unmeasured confounding)
  # Y1
  Y1_rand_pred_obs_rand = Y1_rand_pred # RCT estimate same as CCDS-OR
  Y1_obs_pred_obs_rand = Y1_obs_pred  # obs estimate same as biased estimate in CCDS-OR
  
  # Y2
  Y2_rand_pred_obs_rand = Y2_rand_pred # RCT estimate same as CCDS-OR
  Y2_obs_pred_obs_rand = Y2_obs_pred  # obs estimate same as biased estimate in CCDS-OR
  
  # Y3
  if(!is.na(levels_A[3])){
    Y3_rand_pred_obs_rand = Y3_rand_pred # RCT estimate same as CCDS-OR
    Y3_obs_pred_obs_rand = Y3_obs_pred  # obs estimate same as biased estimate in CCDS-OR
  }
  
  # Target population estimates
  Y1_target_pred_obs_rand = rep(NA,n_target)
  Y1_target_pred_obs_rand[S_target==1] = Y1_rand_pred_obs_rand
  Y1_target_pred_obs_rand[S_target==0] = Y1_obs_pred_obs_rand
  
  Y2_target_pred_obs_rand = rep(NA,n_target)
  Y2_target_pred_obs_rand[S_target==1] = Y2_rand_pred_obs_rand
  Y2_target_pred_obs_rand[S_target==0] = Y2_obs_pred_obs_rand
  
  Y3_target_pred_obs_rand = rep(NA,n_target)
  if(!is.na(levels_A[3])) {
    Y3_target_pred_obs_rand[S_target==1] = Y3_rand_pred_obs_rand
    Y3_target_pred_obs_rand[S_target==0] = Y3_obs_pred_obs_rand
  }
  
  ### Randomized model predictions for all data (ignoring positivity violations)
  # Target population estimates
  Y1_target_pred_rand = predict(fit_rand,data_A1_target_matrix)
  Y2_target_pred_rand = predict(fit_rand,data_A2_target_matrix)
  if(!is.na(levels_A[3])){
    Y3_target_pred_rand = predict(fit_rand,data_A3_target_matrix)
  }
  
  
  ##########################################################################################
  # Check relationship between RCT and obs predictions
  # Sample from sample for plotting feasibility
  sample_target = if(n_target>1000){sample(1:n_target,1000,replace=F)}else{seq(1:n_target)}
  sample_rand = if(n_rand>500){sample(1:n_rand,500,replace=F)}else{seq(1:n_rand)}
  sample_obs = if(n_obs>500){sample(1:n_obs,500,replace=F)}else{seq(1:n_obs)}
  
  # Y1
  plot_CATE_Y1 = cbind.data.frame("X1"=X_target[sample_target,1],
                                  "RCT"=Y1_target_pred_rand[sample_target],
                                  "Obs"=ifelse(S_target==0,Y1_target_pred_obs_rand,NA)[sample_target],
                                  "S"=as.factor(S_target[sample_target])) %>% 
    ggplot(aes(X1,RCT)) +
    annotate("rect", xmin = min_X1_overlap, xmax = max_X1_overlap, ymin = -Inf, ymax = Inf, 
             fill="lightyellow", color=NA, size=0.5) + 
    geom_point(aes(col="RCT"),alpha=0.2)+#color="blue") + 
    geom_point(aes(X1,Obs,col="Obs"),alpha=0.2) +
    geom_smooth(aes(col="RCT"))+
    geom_smooth(aes(X1,Obs,col="Obs"),alpha=0.2) +
    
    ylab("E(Y^a|X1) predictions") + 
    ggtitle("Obs and RCT Y1 model predictions across X1")
  
  
  
  # Y2
  plot_CATE_Y2 = cbind.data.frame("X1"=X_target[sample_target,1],
                                  "RCT"=Y2_target_pred_rand[sample_target],
                                  "Obs"=ifelse(S_target==0,Y2_target_pred_obs_rand,NA)[sample_target],
                                  "S"=as.factor(S_target[sample_target])) %>% 
    ggplot(aes(X1,RCT)) +
    annotate("rect", xmin = min_X1_overlap, xmax = max_X1_overlap, ymin = -Inf, ymax = Inf, 
             fill="lightyellow", color=NA, size=0.5) + 
    geom_point(aes(col="RCT"),alpha=0.2)+#color="blue") + 
    geom_point(aes(X1,Obs,col="Obs"),alpha=0.2) +
    geom_smooth(aes(col="RCT"))+
    geom_smooth(aes(X1,Obs,col="Obs"),alpha=0.2) +
    ylab("E(Y^a|X1) predictions") +  
    ggtitle("Obs and RCT Y2 model predictions across X1")
  
  # Add line for truth depending on DGM/which assumption is being violated
  if(complex_gen_model == F) {
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2,col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2,col="Truth"))
  } else if(data_gen_scenario=="knot"){
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+beta_YX[length(beta_YX)]*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 - 
                                                  15*(X1 < -1)*(X1 + 1),col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+(beta_YX[length(beta_YX)]+beta_YAX[2])*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 - 
                                                  30*(X1 < -1)*(X1 + 1),col="Truth"))
  } else if(data_gen_scenario=="knot2"){
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+beta_YX[length(beta_YX)]*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 - 
                                                  15*(X1 < -0.5)*(X1 + 0.5),col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+(beta_YX[length(beta_YX)]+beta_YAX[2])*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 - 
                                                  30*(X1 < -0.5)*(X1 + 0.5),col="Truth"))
  } else if(data_gen_scenario=="knot3"){
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+beta_YX[length(beta_YX)]*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 - 
                                                  2*(X1 < 0.5)*(X1 - 0.5),col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+(beta_YX[length(beta_YX)]+beta_YAX[2])*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 - 
                                                  4*(X1 < 0.5)*(X1 - 0.5),col="Truth"))
  } else if(data_gen_scenario=="UX1"){
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+beta_YX[length(beta_YX)]*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 + 1*X1,col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+(beta_YX[length(beta_YX)]+beta_YAX[2])*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 + 1.5*X1,col="Truth"))
  } else if(data_gen_scenario=="cbv"){
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+beta_YX[length(beta_YX)]*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2- 
                                                  7.5*(X1 < -0.5)*(X1 + 0.5),col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+(beta_YX[length(beta_YX)]+beta_YAX[2])*(X1+1.0)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2 - 
                                                  15*(X1 < -0.5)*(X1 + 0.5),col="Truth"))
  } else if (complex_gen_model == T) {
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+beta_YX[length(beta_YX)]*(X1+1)^3+#*(X1+1)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2,col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+(beta_YX[length(beta_YX)]+beta_YAX[2])*(X1+1)^3+#*(X1+1)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2,col="Truth"))
    
  } else {
    plot_CATE_Y1 = plot_CATE_Y1 + geom_line(aes(X1,-1.5+beta_YX[1]*X1+beta_YX[length(beta_YX)]*(X1-0.5)^3+#*(X1+1)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2,col="Truth"))
    plot_CATE_Y2 = plot_CATE_Y2 + geom_line(aes(X1,-1.5+beta_YA+(beta_YX[1]+beta_YAX[1])*X1+(beta_YX[length(beta_YX)]+beta_YAX[2])*(X1-0.5)^3+#*(X1+1)^3+
                                                  beta_YU[1]/2+beta_YU[2]*X1/2,col="Truth"))
    
  }
  
  # Plot biases in CCDS estimates
  Y1_rand_pred_overlap = predict(fit_rand_overlap,data_A1_rand_matrix)
  Y2_rand_pred_overlap = predict(fit_rand_overlap,data_A2_rand_matrix)
  
  plot_bias_obs_Y1 = cbind.data.frame(X1 = X_obs[,1], "Obs"= Y1_obs_pred - Y1_obs,
                                      "Obs overlap"=Y1_obs_pred_bias_CCDSa- Y1_obs, 
                                      "RCT overlap"=Y1_obs_pred_bias_CCDSb- Y1_obs, 
                                      "CCDS"=Y1_obs_pred_debiased_CCDS- Y1_obs) %>%  
    mutate(`Obs overlap` = ifelse(X1<min_X1_overlap,NA,`Obs overlap`),
           `RCT overlap` = ifelse(X1<min_X1_overlap,NA,`RCT overlap`)) %>% 
    gather(Model,Bias,-X1) %>% 
    ggplot(aes(X1,Bias,col=Model))+
    annotate("rect", xmin = min_X1_overlap, xmax = max_X1_overlap, ymin = -Inf, ymax = Inf, 
             fill="lightyellow", color=NA, size=0.5) +
    geom_smooth()+geom_hline(aes(yintercept=0)) +  
    ggtitle("Bias in E(Y1|X1) predictions") + scale_color_manual(values=c("purple", "orange", "red", "blue"))
  
  plot_bias_rand_Y1 = cbind.data.frame(X1=X_rand[,1], "RCT overlap" = Y1_rand_pred_overlap-Y1_rand, 
                                       "RCT"=Y1_rand_pred-Y1_rand) %>%  
    mutate(`RCT overlap` = ifelse(X1>max_X1_overlap,NA,`RCT overlap`)) %>%  
    gather(Model,Bias,-X1) %>% 
    ggplot(aes(X1,Bias,col=Model))+
    annotate("rect", xmin = min_X1_overlap, xmax = max_X1_overlap, ymin = -Inf, ymax = Inf, 
             fill="lightyellow", color=NA, size=0.5) +
    geom_smooth()+geom_hline(aes(yintercept=0)) +  
    ggtitle("Bias in E(Y1|X1) predictions") 
  
  plot_bias_obs_Y2 = cbind.data.frame(X1 = X_obs[,1], "Obs"= Y2_obs_pred - Y2_obs,
                                      "Obs overlap"=Y2_obs_pred_bias_CCDSa- Y2_obs, 
                                      "RCT overlap"=Y2_obs_pred_bias_CCDSb- Y2_obs, 
                                      "CCDS"=Y2_obs_pred_debiased_CCDS- Y2_obs) %>%  
    mutate(`Obs overlap` = ifelse(X1<min_X1_overlap,NA,`Obs overlap`),
           `RCT overlap` = ifelse(X1<min_X1_overlap,NA,`RCT overlap`)) %>% 
    gather(Model,Bias,-X1) %>% 
    ggplot(aes(X1,Bias,col=Model))+
    annotate("rect", xmin = min_X1_overlap, xmax = max_X1_overlap, ymin = -Inf, ymax = Inf, 
             fill="lightyellow", color=NA, size=0.5) +
    geom_smooth()+geom_hline(aes(yintercept=0)) +  
    ggtitle("Bias in E(Y2|X1) predictions") + scale_color_manual(values=c("purple", "orange", "red", "blue"))
  
  plot_bias_rand_Y2 = cbind.data.frame(X1=X_rand[,1], "RCT overlap" = Y2_rand_pred_overlap-Y2_rand, 
                                       "RCT"=Y2_rand_pred-Y2_rand) %>%  
    mutate(`RCT overlap` = ifelse(X1>max_X1_overlap,NA,`RCT overlap`)) %>%  
    gather(Model,Bias,-X1) %>% 
    ggplot(aes(X1,Bias,col=Model))+
    annotate("rect", xmin = min_X1_overlap, xmax = max_X1_overlap, ymin = -Inf, ymax = Inf, 
             fill="lightyellow", color=NA, size=0.5) +
    geom_smooth()+geom_hline(aes(yintercept=0)) +  
    ggtitle("Bias in E(Y2|X1) predictions") 
  
  return(list(table_n=table_n,
              table_outcomes=table_outcomes, 
              mean_Y_pop=mean_Y_pop,
              table_overlap=table_overlap,
              plots=list(
                plot_propensity_S_population=plot_propensity_S_population,
                plot_propensity_S_sample=plot_propensity_S_sample,
                plot_log_propensity_S_sample_zoom=plot_log_propensity_S_sample_zoom,
                plot_propensity_S_sample_overlap=plot_propensity_S_sample_overlap,
                plot_log_propensity_S_sample_overlap=plot_log_propensity_S_sample_overlap,
                plot_propensity_S_sample_overlap2=plot_propensity_S_sample_overlap2,
                plot_propensity_S_sample_overlap_zoom=plot_propensity_S_sample_overlap_zoom,
                plot_dist_U_by_group=plot_dist_U_by_group,
                plot_dist_U_S0_by_overlap=plot_dist_U_S0_by_overlap,
                plot_dist_Y_by_group=plot_dist_Y_by_group,
                plot_dist_Y_by_p_S=plot_dist_Y_by_p_S,
                plot_dist_Ya_by_p_S=plot_dist_Ya_by_p_S,
                plot_dist_Ya_by_p_S_group=plot_dist_Ya_by_p_S_group,
                plot_dist_Ya_by_p_A_group=plot_dist_Ya_by_p_A_group,
                plot_CATE_Y1=plot_CATE_Y1,
                plot_CATE_Y2=plot_CATE_Y2,
                plot_bias_obs_Y1=plot_bias_obs_Y1,
                plot_bias_rand_Y1=plot_bias_rand_Y1,
                plot_bias_obs_Y2=plot_bias_obs_Y2,
                plot_bias_rand_Y2=plot_bias_rand_Y2)))
}


#####################################################################################################
#####################################################################################################
## Function to summarize simulation results and produce all output
#####################################################################################################
#####################################################################################################
summarize_simulation_results = function(sim_results, sim_data = sim_data_complex_noX1cube, ensemble=T,
                                        n_estimator = n_estimators){
  # Assess and eliminate iterations that failed
  print("Fraction of iterations that failed") # Note: most succeed when I rerun the simulation with that dataset again, but given how few iterations failed, I didn't think it worthwile. ensemble just happens to not converge for some seeds (https://stackoverflow.com/questions/15895897/line-search-fails-in-training-sim-prob-model)
  print(sum(sapply(sim_results,is.data.frame))/length(sim_results))
  cat("\n")
  
  if(ensemble %in% c(TRUE, "ksvm")){
    if(sum(sapply(sim_results,is.data.frame))/length(sim_results)==0){
      sim_results2 = data.frame(t(sim_results))
    } else {
      failed_indices = which(sapply(sim_results,is.data.frame))
      nonfailing_indices = setdiff(1:length(sim_results),failed_indices)
      
      sim_results2 = sim_results[-failed_indices]
      sim_results2 = do.call(rbind.data.frame, sim_results2) # get in right format
      colnames(sim_results2) = names(sim_results[[nonfailing_indices[1]]])
    }
  } else {
    sim_results2 = sim_results
  }
  
  # Summarize simulation results
  table_sim = get_summary_table(sim_results=sim_results2, sim_data = sim_data) %>% 
    mutate(Estimand = factor(Estimand, levels = c("E(Y^1)","E(Y^2)","E(Y^1) - E(Y^2)"))) %>% 
    arrange(Estimand, Estimator)
  table_plot = get_plotable_results(sim_results=sim_results2)
  boxplot_summary_results_y1y2_subset(table_plot,sim_data) %>% print()
  true_potential_outcomes = c(mean(sim_data$Y1),
                              mean(sim_data$Y2),
                              NA,
                              mean(sim_data$Y1)-mean(sim_data$Y2))
  summary_table_sim = table_sim %>% mutate(Truth = c(rep(true_potential_outcomes[1],n_estimators),
                                                     rep(true_potential_outcomes[2],n_estimators),
                                                     #rep(true_potential_outcomes[3],n_estimators),
                                                     rep(true_potential_outcomes[4],n_estimators)),
                                           Bias = Bias/Truth*100) %>% rename(Percent_bias = Bias)
  cat("\n")
  kable(cbind(summary_table_sim[,c(1,2)], 
              round(summary_table_sim[,c(3:ncol(summary_table_sim))],2)), "html") %>% 
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover")) %>% 
    row_spec(n_estimator, hline_after=T) %>% print()
  cat("\n")
  
  # Summarize algorithm weights
  if(ensemble == TRUE){
    mean_alg_weights_sim = apply(sim_results2[,80:ncol(sim_results2)],2,mean)
    n_table_rows = length(mean_alg_weights_sim)/11
    mean_alg_weights_sim %>% 
      matrix(.,nrow= n_table_rows, ncol=11) %>% 
      set_rownames(sapply(names(mean_alg_weights_sim[1:n_table_rows]), 
                          function(x) str_match(x,"_(.*?)_"))[2,]) %>% 
      set_colnames(c("S", "Roverlap_rand","Roverlap_obs","A_rand","A_obs","A_rand_overlap",
                     "A_obs_overlap","RCT","obs","rand_overlap_","obs_overlap")) %>% 
      round(2) %>% kable("html") %>% 
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover")) %>% 
      row_spec(n_table_rows, hline_after=T) %>% print()
    cat("\n")
  }
  # Overlap and biases
  mean_overlap_sim = c("Overall"=mean(sim_results2$mean_overlap_target),
                       "Obs"=mean(sim_results2$mean_overlap_obs),
                       "RCT"=mean(sim_results2$mean_overlap_rand))
  print("Mean overlap"); cat("\n"); mean_overlap_sim %>% round(2) %>% kable() %>% print()
  cat("\n")
  
  mean_STSM_biases_CCDS_sim = apply(sim_results2[,(56+6):(73+6)],2,mean)
  print("STSM percent biases for CCDS estimator (=RCT estimator for RCT STSM)")
  cat("\n")
  mean_STSM_biases_CCDS_sim %>% matrix(.,nrow=9,ncol=2) %>% 
    set_colnames(.,c("Y1","Y2")) %>% set_rownames(.,c("Obs","Obs observed","Obs overlap","Obs nonoverlap",
                                                      "RCT","RCT observed","RCT overlap", "RCT overlap observed","RCT nonoverlap"))%>% 
    round(2)  %>% kable() %>%  print()
  cat("\n")
  
  # Focus on RCT and RCT overlap 
  print("STSM percent biases for Y1")
  cat("\n")
  summary_STSM_biases_CCDS_complex_linear = apply(sim_results2[,(60+6):(63+6)],2,
                                                  function(x)c("Mean"=mean(x, na.rm=T), 
                                                               quantile(x,0.025, na.rm=T), 
                                                               quantile(x,0.975, na.rm=T)))
  summary_STSM_biases_CCDS_complex_linear %>% set_colnames(.,c("RCT","RCT observed","RCT overlap", "RCT overlap observed"))%>%
    round(2) %>% kable() %>% print()
  cat("\n")
  
  ## Weights
  # Proportion trimmed
  mean_trim = cbind.data.frame("A1"=rbind("w1"=mean(sim_results2$n_trim_w1_A1),
                                          "w2"=mean(sim_results2$n_trim_w2_A1),
                                          "w3"=mean(sim_results2$n_trim_w3_A1),
                                          "w4"=mean(sim_results2$n_trim_w4_A1),
                                          "2-stage"=mean(sim_results2$n_trim_w_bias),
                                          "2-stage WD"=mean(sim_results2$n_trim_w_bias_wd)),
                               "A2"=rbind("w1"=mean(sim_results2$n_trim_w1_A2),
                                          "w2"=mean(sim_results2$n_trim_w2_A2),
                                          "w3"=mean(sim_results2$n_trim_w3_A2),
                                          "w4"=mean(sim_results2$n_trim_w4_A2),
                                          "2-stage"=NA,
                                          "2-stage WD"=NA))
  print("Proportion of weights trimmed");cat("\n"); mean_trim %>% mutate_if(is.numeric, round, 3) %>% kable() %>% print()
  cat("\n")
  
  # Average weight
  mean_weight = cbind.data.frame("A1"=rbind("w1"=mean(sim_results2$w1_A1) %>% round(2),
                                            "w2"=mean(sim_results2$w2_A1) %>% round(2),
                                            "w3"=mean(sim_results2$w3_A1) %>% round(2),
                                            "w4"=mean(sim_results2$w4_A1) %>% round(2),
                                            "2-stage"=mean(sim_results2$w_bias) %>% round(2),
                                            "2-stage WD"=mean(sim_results2$w_bias_wd)) %>% round(2),
                                 "A2"=rbind("w1"=mean(sim_results2$w1_A2) %>% round(2),
                                            "w2"=mean(sim_results2$w2_A2) %>% round(2),
                                            "w3"=mean(sim_results2$w3_A2) %>% round(2),
                                            "w4"=mean(sim_results2$w4_A2) %>% round(2),
                                            "2-stage"=NA,
                                            "2-stage WD"=NA))
  print("Mean weight"); cat("\n"); mean_weight %>% kable() %>% print();cat("\n")
}

