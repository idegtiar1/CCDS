#####################################################################################################
#####################################################################################################
### Helper files for CCDS real data analysis; similar to those for simulation but generalized beyond 3 treatments
#####################################################################################################
#####################################################################################################
# Note: the weighted2stageCCDS corresponds to the 2-stage CCDS presented in the manuscript. 
# The 2-stage HCCDS (hybrid CCDS) and weighted 2-stage HCCDS were not presented in the manuscript as their performance did not justify inclusion.
#####################################################################################################
#####################################################################################################
## Function to estimate PATSMs with bootstrap for confidence intervals
#####################################################################################################
#####################################################################################################
get_estimates_with_bootstrap = function(Y_rand, # outcome vector for randomized data
                                        A_rand, # treatment vector for randomized data. Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                        X_rand, # covariate dataframe for randomized data
                                        Y_obs, # outcome vector for observational data
                                        A_obs, # treatment vector for observational data
                                        X_obs, # covariate dataframe for observational data
                                        X_interaction = NULL, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression for 2-stage estimators
                                        weight_obs = w_obs, # scalar survey weight for obs data to reflect that we're using random subset of observational data but all randomized data
                                        complex_fit_models=F, # whether to use ensemble ("ensemble"), random forest ("ranger"), ksvm ("ksvm"), or linear (F) models to fit models for Y, A, and S
                                        SL_library=my_SL_library, # specifies ensemble SuperLearner library; only used when complex_fit_models = "ensemble"
                                        S_ridge=F, # whether to use ridge regression to fit S model if complex_fit_models=F; otherwise overriden by complex_fit_models argument
                                        overlap_function=F, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4) -- see code below for corresponding alpha/beta
                                        propensity_trim_threshold = 0.001, # threshold at which to trim propensity scores and their products used in weight denominators
                                        M=500, # number of bootstrap resamples
                                        parallel = T, # whether to run parSapply rather than sapply
                                        n_estimators=17, # number of estimators summarized in output of get_estimates() function
                                        include_plot=F, # whether to save plots of propensity score distributions
                                        ...){
  
  ##########################################################################################
  ### Run bootstraps
  ##########################################################################################
  A_target = as.factor(as.character(rbind(A_rand, A_obs)$A))
  levels_A = as.character(unique(A_target)) %>% sort()
  if(parallel){
    boot_results = data.frame(t(parSapply(cl,1:M,
                                       function(i,...){
                                         print(i)
                                         
                                         # Bootstrap sample
                                         r = sample(1:length(Y_rand), length(Y_rand), replace=TRUE)
                                         o = sample(1:length(Y_obs), length(Y_obs), replace=TRUE)
                                         
                                         Y_rand_boot = Y_rand[r]
                                         A_rand_boot = as.data.frame(A_rand[r,]) %>% 
                                           set_colnames("A")# Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                         X_rand_boot = X_rand[r,]
                                         Y_obs_boot = Y_obs[o]
                                         A_obs_boot = as.data.frame(A_obs[o,]) %>% 
                                           set_colnames("A")
                                         X_obs_boot = X_obs[o,]
                                         
                                         out = tryCatch({
                                           get_estimates(Y_rand = Y_rand_boot,
                                                         A_rand = A_rand_boot, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                                         X_rand = X_rand_boot,
                                                         Y_obs = Y_obs_boot,
                                                         A_obs = A_obs_boot,
                                                         X_obs = X_obs_boot,
                                                         X_interaction = X_interaction,
                                                         weight_obs = weight_obs, 
                                                         complex_fit_models=complex_fit_models, 
                                                         SL_library=SL_library, 
                                                         S_ridge=S_ridge, 
                                                         overlap_function=overlap_function, 
                                                         propensity_trim_threshold=propensity_trim_threshold,
                                                         include_plot=include_plot,...)
                                         },
                                         error=function(c) { # Row of NA's if get_estimates() function errors out
                                           if(complex_fit_models=="ensemble"){ # number of output columns for ensemble output
                                             {rep(NA,length(levels_A)*(n_estimators+11)+9+ 
                                                    7*length(SL_library)+4*length(SL_library)*length(levels_A))} # algorithm weight columns
                                           } else { # number of output columns for non-ensemble output
                                             rep(NA,length(levels_A)*(n_estimators+11)+9)
                                           }
                                         })
                                       })))
  } else {
    boot_results = data.frame(t(sapply(1:M,
                                       function(i,...){
                                         print(i)
                                         
                                         # Bootstrap sample
                                         r = sample(1:length(Y_rand), length(Y_rand), replace=TRUE)
                                         o = sample(1:length(Y_obs), length(Y_obs), replace=TRUE)
                                         
                                         Y_rand_boot = Y_rand[r]
                                         A_rand_boot = as.data.frame(A_rand[r,]) %>% 
                                           set_colnames("A") # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                         X_rand_boot = X_rand[r,]
                                         Y_obs_boot = Y_obs[o]
                                         A_obs_boot = as.data.frame(A_obs[o,]) %>% 
                                           set_colnames("A")
                                         X_obs_boot = X_obs[o,]
                                         
                                         out = tryCatch({
                                           get_estimates(Y_rand = Y_rand_boot,
                                                         A_rand = A_rand_boot, # Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                                                         X_rand = X_rand_boot,
                                                         Y_obs = Y_obs_boot,
                                                         A_obs = A_obs_boot,
                                                         X_obs = X_obs_boot,
                                                         weight_obs = weight_obs, 
                                                         X_interaction = X_interaction,
                                                         complex_fit_models=complex_fit_models, 
                                                         SL_library=SL_library, 
                                                         S_ridge=S_ridge, 
                                                         overlap_function=overlap_function, 
                                                         propensity_trim_threshold=propensity_trim_threshold,
                                                         include_plot=include_plot,...)
                                         },
                                         error=function(c) {
                                           if(complex_fit_models=="ensemble"){ # number of output columns for ensemble output
                                             {rep(NA,length(levels_A)*(n_estimators+11)+9+ 
                                                    7*length(SL_library)+4*length(SL_library)*length(levels_A))} # algorithm weight columns
                                           } else { # number of output columns for non-ensemble output
                                             rep(NA,length(levels_A)*(n_estimators+11)+9)
                                           }
                                         })
                                       })))
  }
  
  
  # Re-add names
  names_n_trim = c(paste0("n_trim_w1.", levels_A),
                    paste0("n_trim_w2.", levels_A),
                    paste0("n_trim_w3.", levels_A),
                    paste0("n_trim_w4.", levels_A),
                    "n_trim_w_bias",
                    "n_trim_w_bias_wd")
  
  names_weight_mean = c(paste0("w1.", levels_A),
                         paste0("w2.", levels_A),
                         paste0("w3.", levels_A),
                         paste0("w4.", levels_A),
                         "w_bias",
                         "w_bias_wd")
  
  names_mean_percent_biases = c(paste0("mean_bias_obs_observed.", levels_A),
                                 paste0("mean_bias_rand_observed.", levels_A),
                                 paste0("mean_bias_rand_overlap_observed.", levels_A))
  
  # Re-add names for algorithm weights if using ensemble
  if(complex_fit_models=="ensemble"){
    names_alg_weights_S = paste0("S_", SL_library)
    names_alg_weights_Roverlap_rand = paste0("Roverlap_rand_", SL_library)
    names_alg_weights_Roverlap_obs = paste0("Roverlap_obs_", SL_library)
    names_alg_weights_A_rand = apply(expand.grid(paste0("A_rand_", SL_library), levels_A), 1, paste, collapse=".")
    names_alg_weights_A_obs = apply(expand.grid(paste0("A_obs_", SL_library), levels_A), 1, paste, collapse=".")
    names_alg_weights_A_rand_overlap = apply(expand.grid(paste0("A_rand_overlap_", SL_library), levels_A), 1, paste, collapse=".")
    names_alg_weights_A_obs_overlap = apply(expand.grid(paste0("A_obs_overlap_", SL_library), levels_A), 1, paste, collapse=".")
    names_alg_weights_rand = paste0("rand_", SL_library)
    names_alg_weights_obs = paste0("obs_", SL_library)
    names_alg_weights_rand_overlap = paste0("rand_overlap_", SL_library)
    names_alg_weights_obs_overlap = paste0("obs_overlap_", SL_library)
  } 
    
  
  
  output_names = c(paste0("mean_Ya_rand_pred_rand.", levels_A),
                   paste0("mean_Ya_obs_pred_obs.", levels_A),
                   paste0("mean_Ya_rand_pred_TMLE.", levels_A),
                   paste0("mean_Ya_obs_pred_TMLE.", levels_A),
                   paste0("mean_Ya_rand_pred_AIPW.", levels_A),
                   paste0("mean_Ya_obs_pred_AIPW.", levels_A),
                   paste0("mean_target_pred_CCDS.", levels_A),
                   paste0("mean_target_pred_2stageCCDS.", levels_A),
                   paste0("mean_target_pred_weighted2stageCCDS.", levels_A),
                   paste0("mean_target_pred_2stageHCCDS.", levels_A),
                   paste0("mean_target_pred_weighted2stageHCCDS.", levels_A),
                   paste0("mean_target_pred_CCDS_IPW.", levels_A),
                   paste0("mean_target_pred_CCDS_AIPW.", levels_A),
                   paste0("mean_target_pred_obs_rand.", levels_A),
                   paste0("mean_target_pred_rand.", levels_A),
                   paste0("mean_target_pred_2stageWD.", levels_A),
                   paste0("mean_target_pred_weighted2stageWD.", levels_A),
               
                   "mean_overlap_target",
                   "mean_overlap_obs",
                   "mean_overlap_rand",
                   "min_p_overlap", 
                   "max_p_overlap",
                   
                   names_n_trim,
                   names_weight_mean,
                   names_mean_percent_biases)
  
  if(complex_fit_models=="ensemble"){
    output_names = c(output_names,
               names_alg_weights_S,
               names_alg_weights_Roverlap_rand,
               names_alg_weights_Roverlap_obs,
               names_alg_weights_A_rand,
               names_alg_weights_A_obs,
               names_alg_weights_A_rand_overlap,
               names_alg_weights_A_obs_overlap,
               names_alg_weights_rand,
               names_alg_weights_obs,
               names_alg_weights_rand_overlap,
               names_alg_weights_obs_overlap)
  }
  
  colnames(boot_results) = output_names
  
  return(boot_results)
}


#####################################################################################################
#####################################################################################################
## Function to get estimates from CCDS and comparison estimators
#####################################################################################################
#####################################################################################################
get_estimates = function(Y_rand, # outcome vector for randomized data
                         A_rand, # treatment vector for randomized data. Note: the first treatment should be common (not rare) for more stable estimates, since it'll be used as the reference category
                         X_rand, # covariate dataframe for randomized data
                         Y_obs, # outcome vector for observational data
                         A_obs, # treatment vector for observational data
                         X_obs, # covariate dataframe for observational data
                         X_interaction, # non-main terms in X matrix; i.e., terms in X matrix you don't want to use in second stage regression
                         weight_obs, # scalar; proportion of the target sample that the observational data should constitute; needed to account for us using random subset of observational data but all randomized data
                         complex_fit_models, # whether to use ensemble ("ensemble"), random forest ("ranger") or ksvm ("ksvm") to fit models for Y, A, and S. Note: ksvm will give error if after subsetting by treatment, some X_matrix columns just have 1 value
                         SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                         S_ridge, # whether to use ridge regression to fit S model; overriden by complex_fit_models argument
                         overlap_function, # true overlap region or way to determine overlap (F=default, overlap1, overlap2, overlap3, overlap4)
                         propensity_trim_threshold, # threshold at which to trim propensity scores and their products used in weight denominators
                         include_plot, # whether to save plots of propensity score distributions
                         ...){
  
  # Extract data
  X_target = rbind(X_rand, X_obs)
  if(is.null(X_interaction)){
    Xmain_target = X_target # distinguishes data to use for second stage regressions
  } else {
    Xmain_target = subset(X_target,select=-X_interaction)
  }
  
  X_target_matrix = model.matrix(~.,X_target)[,-1]
  Xmain_target_matrix = model.matrix(~.,Xmain_target)[,-1]
  
  Y_target = c(Y_rand, Y_obs)
  S_target = c(rep(1,length(Y_rand)),rep(0,length(Y_obs)))
  A_target = as.factor(as.character(rbind(A_rand, A_obs)$A))
  levels_A = as.character(unique(A_target)) %>% sort()
  
  # Subset randomized and observational data
  Y_rand = Y_target[S_target==1]
  A_rand = A_target[S_target==1] 
  A_rand_all = model.matrix(~ A_rand, 
                            contrasts.arg=list(A_rand=contrasts(A_rand, contrasts=F)))[,-1] %>% set_colnames(.,levels_A)
  A_rand_matrix = model.matrix(~ ., data.frame(A_rand))[,-1] %>% set_colnames(.,levels_A[-1])
  if(sum(A_rand_matrix[,2])==0){A_rand_matrix=A_rand_matrix[,1]}
  X_rand = X_target[S_target==1,]
  X_rand_matrix = model.matrix(~.,X_rand)[,-1]
  Xmain_rand = Xmain_target[S_target==1,]
  Xmain_rand_matrix = model.matrix(~.,Xmain_rand)[,-1]

  Y_obs = Y_target[S_target==0]
  A_obs = A_target[S_target==0]
  A_obs_all = model.matrix(~ A_obs, 
                            contrasts.arg=list(A_obs=contrasts(A_obs, contrasts=F)))[,-1] %>% set_colnames(.,levels_A)
  A_obs_matrix = model.matrix(~ ., data.frame(A_obs))[,-1] %>% set_colnames(.,levels_A[-1])
  if(sum(A_obs_matrix[,2])==0){A_obs_matrix=A_obs_matrix[,1]}
  X_obs = X_target[S_target==0,]
  X_obs_matrix = model.matrix(~.,X_obs)[,-1]
  Xmain_obs = Xmain_target[S_target==0,]
  Xmain_obs_matrix = model.matrix(~.,Xmain_obs)[,-1]

  # Target sample size and number of randomized and observational individuals
  n_target = length(S_target)
  n_rand = sum(as.numeric(as.character(S_target)))
  n_obs = n_target - n_rand
  
  # Weights to account for using a subset of observational data
  weight_rand = 1 - weight_obs
  weight_sampling = ifelse(S_target==1, weight_rand, weight_obs)
  
  # Remove unnecessary data
  rm(Xmain_target, Xmain_target_matrix)
  gc()
  
  ##########################################################################################
  ### Determine overlap and non-overlap regions
  ### Use method by Nethery et al. 2019
  ##########################################################################################
  # Propensity score for selection: P(S|X)
  if(complex_fit_models=="ensemble"){
    fit_S = SuperLearner(Y=S_target, X=as.data.frame(X_target_matrix), family="binomial", SL.library=SL_library, verbose=F) 
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = predict(fit_S, as.data.frame(X_target_matrix), onlySL = TRUE)$pred #as.data.frame(X_target_matrix)
    
  } else if(complex_fit_models=="ksvm"){
    fit_S = cbind.data.frame(S_target,X_target_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = predict(fit_S, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
  } else if(complex_fit_models=="ranger"){ # Ranger: random forest
    fit_S = cbind.data.frame(S_target,X_target_matrix) %>% ranger(formula(.), data=.,
                                                                  num.trees = 300, min.node.size=floor(nrow(X_target_matrix)*0.05), # these hyperparameters work reasonably in tests
                                                                  num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = cbind.data.frame(S_target,X_target_matrix) %>%
      predict(fit_S, data=.,type = "response") %>% .$predictions %>% .[,2]
    
  } else if(complex_fit_models=="bart"){
    fit_S = bart(X_target_matrix, S_target, X_target_matrix, verbose = F)
    pi_S = apply(pnorm(fit_S$yhat.test),2,mean)
    
  } else if(S_ridge==T){
    SX_target <- cbind.data.frame(S_target,X_target_matrix)
    fit_S <- cv.glmnet(X_target_matrix,S_target, data = SX_target, alpha=0,family="binomial")
    
    pi_S <- predict(fit_S,X_target_matrix, type="response") # predicted probability that y=1
    
  } else {
    fit_S = cbind.data.frame(S_target,X_target_matrix) %>% glm(formula(.), data=.,family="binomial")
    # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
    pi_S = predict(fit_S,type = "response")
    
  }
  
  # Trim pi_S 
  pi_S_trim = ifelse(pi_S<propensity_trim_threshold,propensity_trim_threshold,
                     ifelse(pi_S>(1-propensity_trim_threshold),(1-propensity_trim_threshold),pi_S))
  
  # logit propensity score for selection
  logit = qlogis
  logit_pi_S = logit(pi_S_trim)
  
  # Estimate overlap
  # a, b = tuning parameters for determining region of overlap
  if(overlap_function == "overlap1"){ # uses nominal propensity score (if changed, change plot that "a" feeds into)
    a = (max(pi_S)-min(pi_S)) %>% abs() %>% as.numeric()*.01# # a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 1% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=pi_S, # propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function == "overlap2"){ # uses nominal propensity score (if changed, change plot that "a" feeds into)
    a = 0.05
    b = 10
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=pi_S, # propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function == "overlap3"){ # uses logit propensity score (if changed, change plot that "a" feeds into)
    a = (max(logit_pi_S)-min(logit_pi_S)) %>% abs() %>% as.numeric()*.1
    b = round(0.04*min(n_obs,n_rand),0) # less than 4% of randomized units in 10% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=logit_pi_S, # logit propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function == "overlap4"){ # uses logit propensity score (if changed, change plot that "a" feeds into)
    a = (max(logit_pi_S)-min(logit_pi_S)) %>% abs() %>% as.numeric()*.02
    b = round(0.01*min(n_obs,n_rand),0) # less than 1% of randomized units in 2% range
    
    # Determine region of overlap
    overlap_target = pw_overlap(ps=logit_pi_S, # logit propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else if(overlap_function==F){ # uses logit propensity score (if changed, change plot that "a" feeds into)
    # Default overlap estimation approach
    # Given the skewdness of the data, determining overlap on logit scale
    a = (max(logit_pi_S)-min(logit_pi_S)) %>% abs() %>% as.numeric()*.01# # a=0.1*range(ps) suggested by Nethery et al. 2019; given that p_S are all small because a small proportion of the observations are randomized (vs. their simulations had an even split), I will choose a as the 1% quantile of propensity scores
    b = max(round(0.01*min(n_obs,n_rand),0),1) # less than 1% of randomized units in 1% range. Nethery et al. 2019 recommended 10 and used simulations with 250 in each treatment arm; 10/250=4%
    
    # Determine region of overlap -- slow
    overlap_target = pw_overlap(ps=logit_pi_S, # logit propensity score for S
                                E=S_target, # "Exposure", for us, S
                                a=a, # tuning parameter a
                                b=b) # tuning parameter b
  } else { # known overlap region
    overlap_target=eval(parse(text=overlap_function))
  }
  
  # Regions of overlap
  overlap_rand = overlap_target[S_target==1]
  overlap_obs = overlap_target[S_target==0]
  X_rand_overlap = X_rand[overlap_rand==1,]
  X_obs_overlap = X_obs[overlap_obs==1,]
  X_rand_overlap_matrix = model.matrix(~.,X_rand_overlap)[,-1]
  X_obs_overlap_matrix = model.matrix(~.,X_obs_overlap)[,-1]
  
  Xmain_rand_overlap = Xmain_rand[overlap_rand==1,]
  Xmain_obs_overlap = Xmain_obs[overlap_obs==1,]
  Xmain_rand_overlap_matrix = model.matrix(~.,Xmain_rand_overlap)[,-1]
  Xmain_obs_overlap_matrix = model.matrix(~.,Xmain_obs_overlap)[,-1]
  
  # Proportions of data in region of overlap
  mean_overlap_target = mean(overlap_target) # proportion of target sample data in region of overlap
  mean_overlap_rand = mean(overlap_rand) # proportion of randomized data in region of overlap
  mean_overlap_obs = mean(overlap_obs) # proportion of observational data in region of overlap
  
  # Bounds of overlap region (ends up being contiguous)
  order_pi_S = pi_S[order(pi_S)]
  order_overlap_target = overlap_target[order(pi_S)]
  order_S_target = S_target[order(S_target)]
  min_p_overlap = order_pi_S[min(which(order_overlap_target==1))] # smallest propensity score in the region of overlap
  max_p_overlap = order_pi_S[max(which(order_overlap_target==1))] # largest propensity score in the region of overlap
  
  # Remove unnecessary data
  rm(X_rand_overlap, X_obs_overlap, Xmain_rand, Xmain_rand_overlap, Xmain_obs_overlap)
  gc()
  
  # Propensity plots
  if(include_plot == T){
    # Label ending for saving plot
    ending = paste0(ifelse(complex_fit_models == F, "_linear", paste0("_",complex_fit_models)), 
                    ifelse(overlap_function == F, "", paste0("_",overlap_function)))
    
    # P(S) plot
    plot_propensity_S_sample_overlap = cbind.data.frame(S_target,pi_S) %>% 
      mutate(S_target = factor(S_target, levels=c(0,1), labels=c("Observational","Randomized"))) %>% 
      ggplot(aes(pi_S, fill=S_target)) +
      annotate("rect", xmin = min_p_overlap, xmax = max_p_overlap, ymin = 0, ymax = Inf, 
               fill="pink", color=NA, size=0.5) + 
      ylab("Percentage") + theme_bw() + ggtitle("Propensity for S, overlap region in light pink") +
      labs(fill="") 
    
    # logit(P(S)) plot
    plot_log_propensity_S_sample_overlap = cbind.data.frame(S_target,logit_pi_S) %>% 
      mutate(S_target = factor(S_target, levels=c(0,1), labels=c("Observational","Randomized"))) %>% 
      ggplot(aes(logit_pi_S, fill=S_target)) +
      annotate("rect", xmin = logit(min_p_overlap), xmax = logit(max_p_overlap), ymin = 0, ymax = Inf, 
               fill="pink", color=NA, size=0.5) + 
      ylab("Percentage") + theme_bw() + ggtitle("Log propensity for S, overlap region in light pink") +
      labs(fill="") 
    
    # Add binwidth = a to histogram specifications depending on whether a determined on propensity or log propensity scale
    if(overlap_function %in% c("overlap1", "overlap2")){ # a and b determined on nominal scale
      plot_propensity_S_sample_overlap = plot_propensity_S_sample_overlap + 
        geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                               ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                       position='dodge', binwidth = a)
      ggsave(paste0("/Results/Generalizability_results/plot_overlap",ending,".png"))
      
      plot_log_propensity_S_sample_overlap = plot_log_propensity_S_sample_overlap +
        geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                               ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                       position='dodge')
      ggsave(paste0("/Results/Generalizability_results/plot_log_overlap",ending,".png"))
      
    } else if (overlap_function %in% c(F, "overlap3", "overlap4")){ # a and b determined on log scale
      plot_propensity_S_sample_overlap = plot_propensity_S_sample_overlap + 
        geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                               ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                       position='dodge')
      ggsave(paste0("/Results/Generalizability_results/plot_overlap",ending,".png"))
      
      plot_log_propensity_S_sample_overlap = plot_log_propensity_S_sample_overlap +
        geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                               ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                       position='dodge', binwidth = a)
      ggsave(paste0("/Results/Generalizability_results/plot_log_overlap",ending,".png"))
      
    } else { # overlap region prespecified, a and b not used
      plot_propensity_S_sample_overlap = plot_propensity_S_sample_overlap + 
        geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                               ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                       position='dodge')
      ggsave(paste0("/Results/Generalizability_results/plot_overlap",ending,".png"))
      
      plot_log_propensity_S_sample_overlap = plot_log_propensity_S_sample_overlap +
        geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                               ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                       position='dodge')
      ggsave(paste0("/Results/Generalizability_results/plot_log_overlap",ending,".png"))
      
    }
    
  }
  
  ##########################################################################################
  ### Determine propensities and weights used in estimators 
  ##########################################################################################
  if(complex_fit_models=="ensemble"){
    # P(R_overlap|S=1,X) 
    fit_Roverlap_rand = SuperLearner(Y=overlap_rand, as.data.frame(X_rand_matrix),  
                                     family="binomial", SL.library=SL_library, verbose=F)
    pi_Roverlap_rand = predict(fit_Roverlap_rand, as.data.frame(X_target_matrix), onlySL = TRUE)$pred 
    
    # P(R_overlap|S=0,X)
    fit_Roverlap_obs = SuperLearner(Y=overlap_obs, X=as.data.frame(X_obs_matrix), 
                                    family="binomial", SL.library=SL_library, verbose=F)
    pi_Roverlap_obs = predict(fit_Roverlap_obs, as.data.frame(X_target_matrix), onlySL = TRUE)$pred #
    
    # P(A|S=1,X)
    # fit separately for each A then normalize because multinomial SuperLearner currently doesn't exist
    fit_A_rand = apply(A_rand_all, 2,
                           function(a) SuperLearner(Y=a, X=as.data.frame(X_rand_matrix), 
                                                    SL.library=SL_library, 
                                                    family="binomial",method="method.NNLS",
                                                    verbose=F))
    pi_A_rand_not_normalized = sapply(fit_A_rand, function(a) predict(a, as.data.frame(X_target_matrix), onlySL = TRUE)$pred)
    
    
    # P(A|S=0,X)
    fit_A_obs = apply(A_obs_all, 2,
                       function(a) SuperLearner(Y=a, X=as.data.frame(X_obs_matrix), 
                                                SL.library=SL_library, 
                                                family="binomial",method="method.NNLS",
                                                verbose=F))
    pi_A_obs_not_normalized = sapply(fit_A_obs, function(a) predict(a, as.data.frame(X_target_matrix), onlySL = TRUE)$pred)
    

    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap = apply(A_rand_all[overlap_rand==1,], 2,
                               function(a) SuperLearner(Y=a, X=as.data.frame(X_rand_matrix)[overlap_rand==1,], 
                                      family="binomial", SL.library=SL_library, verbose=F))
    pi_A_rand_overlap_not_normalized = sapply(fit_A_rand_overlap, function(a) predict(a, as.data.frame(X_target_matrix), onlySL = TRUE)$pred)
    
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap = apply(A_obs_all[overlap_obs==1,], 2, 
                              function(a) SuperLearner(Y=a, X=as.data.frame(X_obs_matrix)[overlap_obs==1,], 
                                     family="binomial", SL.library=SL_library, verbose=F))
    pi_A_obs_overlap_not_normalized = sapply(fit_A_obs_overlap, function(a) predict(a, as.data.frame(X_target_matrix), onlySL = TRUE)$pred)
    
    
    ### Obtain normalized probabilities
    pi_A_rand = pi_A_rand_not_normalized/apply(pi_A_rand_not_normalized,1,sum)
    pi_A_obs = pi_A_obs_not_normalized/apply(pi_A_obs_not_normalized,1,sum)
    pi_A_rand_overlap = pi_A_rand_overlap_not_normalized/apply(pi_A_rand_overlap_not_normalized,1,sum)
    pi_A_obs_overlap = pi_A_obs_overlap_not_normalized/apply(pi_A_obs_overlap_not_normalized,1,sum)
  
  } else if(complex_fit_models=="ksvm"){
    # P(R_overlap|S=1,X)
    fit_Roverlap_rand = cbind.data.frame("overlap"=overlap_rand,X_rand_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_Roverlap_rand = predict(fit_Roverlap_rand, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    # P(R_overlap|S=0,X)
    fit_Roverlap_obs = cbind.data.frame("overlap"=overlap_obs,X_obs_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_Roverlap_obs = predict(fit_Roverlap_obs, as.data.frame(X_target_matrix),type="probabilities")[,2]
    
    # P(A|S=1,X)
    fit_A_rand = cbind.data.frame("A"=droplevels(A_rand),X_rand_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_A_rand = predict(fit_A_rand, as.data.frame(X_target_matrix),type="probabilities")
    
    # P(A|S=0,X)
    fit_A_obs = cbind.data.frame("A"=droplevels(A_obs),X_obs_matrix) %>% ksvm(formula(.), data=.,type="C-svc",prob.model=T)
    pi_A_obs = predict(fit_A_obs, as.data.frame(X_target_matrix),type="probabilities")
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap = cbind.data.frame("A"=droplevels(A_rand)[overlap_rand==1],
                                          X_rand[overlap_rand==1,]) %>% ksvm(A ~ ., data=.,type="C-svc",prob.model=T)
    pi_A_rand_overlap = predict(fit_A_rand_overlap, as.data.frame(X_target_matrix),type="probabilities")
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap = cbind.data.frame("A"=droplevels(A_obs)[overlap_obs==1],
                                         X_obs[overlap_obs==1,]) %>% ksvm(A ~ ., data=.,type="C-svc",prob.model=T)
    pi_A_obs_overlap = predict(fit_A_obs_overlap, as.data.frame(X_target_matrix),type="probabilities")
    
  } else if(complex_fit_models=="ranger"){
    # P(R_overlap|S=1,X)
    fit_Roverlap_rand = cbind.data.frame("overlap"=overlap_rand,X_rand_matrix) %>% 
      ranger(formula(.), data=.,
             num.trees = 300, min.node.size=floor(nrow(X_rand_matrix)*0.05), # these hyperparameters work reasonably in tests
             num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_Roverlap_rand = cbind.data.frame(X_target_matrix) %>%
      predict(fit_Roverlap_rand, data=.,type = "response") %>% .$predictions %>% .[,2]
    
    # P(R_overlap|S=0,X)
    fit_Roverlap_obs = cbind.data.frame("overlap"=overlap_obs,X_obs_matrix) %>% 
      ranger(formula(.), data=.,
             num.trees = 300, min.node.size=floor(nrow(X_obs_matrix)*0.05), # these hyperparameters work reasonably in tests
             num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_Roverlap_obs = cbind.data.frame(X_target_matrix) %>%
      predict(fit_Roverlap_obs, data=.,type = "response") %>% .$predictions %>% .[,2]
    
    # P(A|S=1,X)
    fit_A_rand = cbind.data.frame("A"=A_rand,X_rand_matrix) %>% 
      ranger(formula(.), data=.,
             num.trees = 300, min.node.size=floor(nrow(X_rand_matrix)*0.05), # these hyperparameters work reasonably in tests
             num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_rand = cbind.data.frame(X_target_matrix) %>%
      predict(fit_A_rand, data=.,type = "response") %>% .$predictions 
    
    # P(A|S=0,X)
    fit_A_obs = cbind.data.frame("A"=A_obs,X_obs_matrix) %>% 
      ranger(formula(.), data=.,
             num.trees = 300, min.node.size=floor(nrow(X_obs_matrix)*0.05), # these hyperparameters work reasonably in tests
             num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_obs = cbind.data.frame(X_target_matrix) %>%
      predict(fit_A_obs, data=.,type = "response") %>% .$predictions 
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap = cbind.data.frame("A"=A_rand[overlap_rand==1],
                                          X_rand[overlap_rand==1,]) %>% 
      ranger(formula(.), data=.,
             num.trees = 300, min.node.size=floor(nrow(X_rand[overlap_rand==1,])*0.05), # these hyperparameters work reasonably in tests
             num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_rand_overlap = as.data.frame(X_target_matrix) %>%
      predict(fit_A_rand_overlap, data=.,type = "response") %>% .$predictions 
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap = cbind.data.frame("A"=A_obs[overlap_obs==1],
                                         X_obs[overlap_obs==1,]) %>% 
      ranger(formula(.), data=.,
             num.trees = 300, min.node.size=floor(nrow(X_obs[overlap_obs==1,])*0.05), # these hyperparameters work reasonably in tests
             num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
    pi_A_obs_overlap = as.data.frame(X_target_matrix) %>%
      predict(fit_A_obs_overlap, data=.,type = "response") %>% .$predictions 
    
  } else if(complex_fit_models=="bart"){
    # P(R_overlap|S=1,X)
    fit_Roverlap_rand = cbind.data.frame("overlap"=overlap_rand,X_rand_matrix) %>% 
      bart2(formula(.), data=., test=as.data.frame(X_target_matrix), verbose = F)
    pi_Roverlap_rand = apply(pnorm(fit_Roverlap_rand$yhat.test),3,mean)
    
    # P(R_overlap|S=0,X)
    fit_Roverlap_obs = cbind.data.frame("overlap"=overlap_obs,X_obs_matrix) %>% 
      bart2(formula(.), data=.,test=as.data.frame(X_target_matrix), verbose = F)
    pi_Roverlap_obs = apply(pnorm(fit_Roverlap_obs$yhat.test),3,mean)
    
    # P(A|S=1,X)
    fit_A_rand = mbart2(x.train = X_rand_matrix, 
                        y.train = droplevels(A_rand), 
                        x.test = X_target_matrix[,colnames(X_rand_matrix)],
                        offset = rep(0, length(levels_A)), ndpost=100L, printevery=2000L) #
    
    # Separate out predictions for each A=a: split prob.test into every ath element, then take posterior means
    pi_A_rand = sapply((1:length(levels_A))-1, 
                       function(a) apply(fit_A_rand$prob.test[,seq(1+a, ncol(fit_A_rand$prob.test),
                                                            length(levels_A))],2,mean)) %>% 
      set_colnames(levels_A)
    

    # P(A|S=0,X)
    fit_A_obs = mbart2(x.train = X_obs_matrix, 
                       y.train = droplevels(A_obs), 
                       x.test = X_target_matrix[,colnames(X_obs_matrix)],
                        offset = rep(0, length(levels_A)), ndpost=100L, printevery=2000L) #
    
    # Separate out predictions for each A=a: split prob.test into every ath element, then take posterior means
    pi_A_obs = sapply((1:length(levels_A))-1, 
                       function(a) apply(fit_A_obs$prob.test[,seq(1+a, ncol(fit_A_obs$prob.test),
                                                                   length(levels_A))],2,mean)) %>% 
      set_colnames(levels_A)
    
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap = mbart2(x.train = X_rand_matrix[overlap_rand==1,], 
                                y.train = droplevels(A_rand)[overlap_rand==1], 
                                x.test = X_target_matrix[,colnames(X_rand_matrix)],
                        offset = rep(0, length(levels_A)), ndpost=100L, printevery=2000L) #
    
    # Separate out predictions for each A=a: split prob.test into every ath element, then take posterior means
    pi_A_rand_overlap = sapply((1:length(levels_A))-1, 
                       function(a) apply(fit_A_rand_overlap$prob.test[,seq(1+a, ncol(fit_A_rand_overlap$prob.test),
                                                                   length(levels_A))],2,mean)) %>% 
      set_colnames(levels_A)
    
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap = mbart2(x.train = X_obs_matrix[overlap_obs==1,], 
                                y.train = droplevels(A_obs)[overlap_obs==1], 
                                x.test = X_target_matrix[,colnames(X_obs_matrix)],
                                offset = rep(0, length(levels_A)), ndpost=100L, printevery=2000L) #
    
    # Separate out predictions for each A=a: split prob.test into every ath element, then take posterior means
    pi_A_obs_overlap = sapply((1:length(levels_A))-1, 
                               function(a) apply(fit_A_obs_overlap$prob.test[,seq(1+a, ncol(fit_A_obs_overlap$prob.test),
                                                                                   length(levels_A))],2,mean)) %>% 
      set_colnames(levels_A)
    
  } else {
    # P(R_overlap|S=1,X)
    fit_Roverlap_rand =  cbind.data.frame("overlap"=overlap_rand,X_rand_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_Roverlap_rand = predict(fit_Roverlap_rand,as.data.frame(X_target_matrix),type = "response")
    
    # P(R_overlap|S=0,X)
    fit_Roverlap_obs =  cbind.data.frame("overlap"=overlap_obs,X_obs_matrix) %>% glm(formula(.), data=.,family="binomial")
    pi_Roverlap_obs = predict(fit_Roverlap_obs,as.data.frame(X_target_matrix),type = "response")
    
    # P(A|S=1,X) - eliminate residual noise through estimating although propensity is known
    fit_A_rand <- cbind.data.frame(A_rand,X_rand_matrix) %>% nnet::multinom(formula(.), data=., trace=F)
    pi_A_rand = predict(fit_A_rand,as.data.frame(X_target_matrix),type = "probs")
    
    # P(A|S=0,X)
    fit_A_obs =  cbind.data.frame("A"=A_obs,X_obs_matrix) %>% nnet::multinom(formula(.), data=., trace=F)
    pi_A_obs = predict(fit_A_obs,as.data.frame(X_target_matrix),type = "probs")
    
    # P(A|S=1,R_overlap=1,X)
    fit_A_rand_overlap =  cbind.data.frame("A"=A_rand[overlap_rand==1],
                                           X_rand_overlap_matrix) %>% nnet::multinom(formula(.), data=., trace=F)
    pi_A_rand_overlap = predict(fit_A_rand_overlap,as.data.frame(X_target_matrix),type = "probs")
    
    # P(A|S=0,R_overlap=1,X)
    fit_A_obs_overlap =  cbind.data.frame("A"=A_obs[overlap_obs==1],
                                          X_obs_overlap_matrix) %>% nnet::multinom(formula(.), data=., trace=F)
    pi_A_obs_overlap = predict(fit_A_obs_overlap,as.data.frame(X_target_matrix),type = "probs")
    
  } 
  
  
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
  w3_denom = sapply(1:length(levels_A), function(i) pi_Roverlap_obs_trim*(1-pi_A_obs_overlap_trim[,i])) %>% 
    set_colnames(levels_A)
  
  w4_denom = sapply(1:length(levels_A), function(i) pi_S_trim*pi_Roverlap_rand_trim*(1-pi_A_rand_overlap_trim[,i])) %>% 
    set_colnames(levels_A)
  
  # Trim
  w3_denom_trim = ifelse(w3_denom<propensity_trim_threshold,propensity_trim_threshold,
                                ifelse(w3_denom>(1-propensity_trim_threshold),(1-propensity_trim_threshold),w3_denom))
  w4_denom_trim = ifelse(w4_denom<propensity_trim_threshold,propensity_trim_threshold,
                         ifelse(w4_denom>(1-propensity_trim_threshold),(1-propensity_trim_threshold),w4_denom))
  
  ### Calculate weights
  # Untrimmed weights
  w1 = sapply(1:length(levels_A), function(i) ifelse(S_target==1 & A_target==levels_A[i],
                                                     1/(1-pi_A_rand[,i]),0)) %>% 
    set_colnames(levels_A)
  w2 = sapply(1:length(levels_A), function(i) ifelse(S_target==0 & A_target==levels_A[i],
                                                     1/(1-pi_A_obs[,i]),0)) %>% 
    set_colnames(levels_A)
  w3 = sapply(1:length(levels_A), function(i) ifelse(S_target==0 & A_target==levels_A[i] & overlap_target==1,
                                                     1/(w3_denom[,i]),0)) %>% 
    set_colnames(levels_A)
  w4 = sapply(1:length(levels_A), function(i) ifelse(S_target==1 & A_target==levels_A[i] & overlap_target==1,
                                                     (1-pi_S)/(w4_denom[,i]),0)) %>% 
    set_colnames(levels_A)
    
  # Trimmed weights
  w1_Aa_trim = sapply(1:length(levels_A), function(i) ifelse(S_target==1 & A_target==levels_A[i],
                                                     1/(1-pi_A_rand_trim[,i]),0)) %>% 
    set_colnames(levels_A)
  w2_Aa_trim = sapply(1:length(levels_A), function(i) ifelse(S_target==0 & A_target==levels_A[i],
                                                     1/(1-pi_A_obs_trim[,i]),0)) %>% 
    set_colnames(levels_A)
  w3_Aa_trim = sapply(1:length(levels_A), function(i) ifelse(S_target==0 & A_target==levels_A[i] & overlap_target==1,
                                                     1/(w3_denom_trim[,i]),0)) %>% 
    set_colnames(levels_A)
  w4_Aa_trim = sapply(1:length(levels_A), function(i) ifelse(S_target==1 & A_target==levels_A[i] & overlap_target==1,
                                                     (1-pi_S_trim)/(w4_denom_trim[,i]),0)) %>% 
    set_colnames(levels_A)
  
  
  # Summarize proportion of observations with trimmed weights
  n_trim_w1 = sapply(1:length(levels_A), function(i) mean(w1[,i] != w1_Aa_trim[,i]))
  n_trim_w2 = sapply(1:length(levels_A), function(i) mean(w2[,i] != w2_Aa_trim[,i]))
  n_trim_w3 = sapply(1:length(levels_A), function(i) mean(w3[,i] != w3_Aa_trim[,i] |
                                                            pi_Roverlap_obs != pi_Roverlap_obs_trim |
                                                            pi_A_obs_overlap[,i] != pi_A_obs_overlap_trim[,i]))
  n_trim_w4 = sapply(1:length(levels_A), function(i) mean(w4[,i] != w4_Aa_trim[,i] |
                                                            pi_S != pi_S_trim |
                                                            pi_Roverlap_rand != pi_Roverlap_rand_trim |
                                                            pi_A_rand_overlap[,i] != pi_A_rand_overlap_trim[,i]))

  
  ##########################################################################################
  ### Model conditional distributions used in estimators 
  ##########################################################################################
  if(complex_fit_models=="ensemble"){
    # Estimate restricted mean model E(Y|S=1,A=a,X) within the randomized data
    fit_rand <- SuperLearner(Y=Y_rand, X=cbind.data.frame(A_rand_matrix,as.data.frame(X_rand_matrix)), 
                             family="gaussian", SL.library=SL_library, verbose=F)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X) within the observational data
    fit_obs = SuperLearner(Y=Y_obs, X=cbind.data.frame(A_obs_matrix,as.data.frame(X_obs_matrix)), 
                           family="gaussian", SL.library=SL_library, verbose=F)
    
    # Estimate restricted mean model E(Y|S=1,A=a,X_overlap,X) within the randomized data in the overlap region
    fit_rand_overlap = SuperLearner(Y=Y_rand[overlap_rand==1], 
                                    X=cbind.data.frame(A_rand_matrix,as.data.frame(X_rand_matrix))[overlap_rand==1,], 
                                    family="gaussian", SL.library=SL_library, verbose=F)
    
    # Estimate restricted mean model E(Y|S=0,A=a,X_overlap,X) within the observational data in the overlap region
    fit_obs_overlap = SuperLearner(Y=Y_obs[overlap_obs==1], 
                                   X=cbind.data.frame(A_obs_matrix,as.data.frame(X_obs_matrix))[overlap_obs==1,], 
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
             num.threads=1, respect.unordered.factors=T)#, seed=my_seed)
    
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
    
    
  } else if(complex_fit_models=="bart"){
    ## Using BART package because dbarts version of bart subsets to non-collinear columns in a non-reproducible way (uncertain which columns it subsets to), while I can get column names in the BART package version
    fit_rand <- gbart(x.train=cbind.data.frame(A_rand_all,X_rand_matrix), 
                      y.train=Y_rand, type = "wbart", printevery=2000L)
    
    fit_obs <- gbart(x.train=cbind.data.frame(A_obs_all,X_obs_matrix), 
                      y.train=Y_obs, type = "wbart", printevery=2000L)
    
    fit_rand_overlap <- gbart(x.train=cbind.data.frame(A_rand_all,X_rand_matrix)[overlap_rand==1,],
                              y.train=Y_rand[overlap_rand==1], type = "wbart", printevery=2000L)
    
    fit_obs_overlap <- gbart(x.train=cbind.data.frame(A_obs_all,X_obs_matrix)[overlap_obs==1,], 
                              y.train=Y_obs[overlap_obs==1], type = "wbart", printevery=2000L)
    
  } else { # Include interactions between all covariates and A, changing names for common names across data.frames so predict() function works
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
  # Set up potential outcome target data with the same variable names for prediction
  # Matrix format depends on algorithm(s) used to fit regressions
  data_Aa_target = lapply(1:length(levels_A), function(i) cbind.data.frame("A"=levels_A[i],X_target)) %>% 
    set_names(levels_A)# target data with A=a for all observations
  
  if(complex_fit_models %in% c("ensemble")){
    data_A0_target_matrix = matrix(data=0,nrow=n_target,ncol=length(levels_A)-1) %>% set_colnames(levels_A[-1]) %>% 
      cbind.data.frame(.,X_target_matrix) # target data with A=0 for all observations
    data_Aa_target_matrix = lapply(1:(length(levels_A)-1), function(i){ data_A0_target_matrix[,i]=1; return(data_A0_target_matrix)}) %>% 
      prepend(.,list(data_A0_target_matrix)) %>% 
      set_names(levels_A)
    
  } else if(complex_fit_models %in% c("bart")){
    data_Aa_target_matrix = lapply(1:(length(levels_A)), function(i){ 
      dat = cbind.data.frame(matrix(rep(0, length(levels_A)), ncol=length(levels_A)),
                 X_target_matrix) %>% 
        set_colnames(c(levels_A, colnames(X_target_matrix))) 
      dat[,i] = 1
      return(dat)}) %>% 
      set_names(levels_A)
  } else {
    data_Aa_target_matrix = lapply(1:length(levels_A), function(i) cbind.data.frame("A"=levels_A[i],X_target_matrix)) %>% 
      set_names(levels_A) # target data with A=a for all observations
  }
  
  # Subset data needed for estimators
  data_X_obs = as.data.frame(cbind(X_obs))
  data_X_obs_matrix = model.matrix(~.,data_X_obs)[,-1]
  
  data_Xmain_obs = as.data.frame(cbind(Xmain_obs))
  data_Xmain_obs_matrix = model.matrix(~.,data_Xmain_obs)[,-1]
  
  
  # Needed for CCDS estimators
  # Estimator name indicates data being used to fit regression (e.g., obs_overlap) 
  data_Aa_rand = lapply(1:(length(levels_A)), function(i) data_Aa_target[[i]][S_target==1,])
  data_Aa_obs = lapply(1:(length(levels_A)), function(i) data_Aa_target[[i]][S_target==0,])
  
  data_Aa_rand_matrix = lapply(1:(length(levels_A)), function(i) data_Aa_target_matrix[[i]][S_target==1,])
  data_Aa_obs_matrix = lapply(1:(length(levels_A)), function(i) data_Aa_target_matrix[[i]][S_target==0,])
  
  # Needed for 2-stage WD
  data_Aa_rand_overlap = lapply(1:(length(levels_A)), function(i) data_Aa_target[[i]][S_target==1 & overlap_target==1,])
  data_Aa_rand_matrix_overlap = lapply(1:(length(levels_A)), function(i) data_Aa_target_matrix[[i]][S_target==1 & overlap_target==1,])
  
  
  ##########################################################################################
  ##########################################################################################
  ### Estimate target population means (PTSMs)
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  ### Estimate potential outcomes: CCDS 
  ##########################################################################################
  # Estimate terms in CCDS estimators
  if(complex_fit_models=="ensemble"){
    Ya_rand_pred = sapply(1:(length(levels_A)), function(i) predict(fit_rand,data_Aa_rand_matrix[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A) 
    Ya_obs_pred = sapply(1:(length(levels_A)), function(i) predict(fit_obs,data_Aa_obs_matrix[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSa = sapply(1:(length(levels_A)), function(i) predict(fit_obs_overlap,data_Aa_obs_matrix[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i) predict(fit_rand_overlap,data_Aa_obs_matrix[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDS = Ya_obs_pred_bias_CCDSa - Ya_obs_pred_bias_CCDSb  # bias estimate for obs data
    Ya_obs_pred_debiased_CCDS = Ya_obs_pred - # obs estimate
      Ya_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Ya_rand_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i) predict(fit_rand_overlap,data_Aa_rand_matrix[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A) 
    
  } else if(complex_fit_models=="ranger"){
    Ya_rand_pred = sapply(1:(length(levels_A)), function(i) predict(fit_rand,data_Aa_rand[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_obs_pred = sapply(1:(length(levels_A)), function(i) predict(fit_obs,data_Aa_obs[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSa = sapply(1:(length(levels_A)), function(i) predict(fit_obs_overlap,data_Aa_obs[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i) predict(fit_rand_overlap,data_Aa_obs[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDS = Ya_obs_pred_bias_CCDSa - Ya_obs_pred_bias_CCDSb  # bias estimate for obs data
    Ya_obs_pred_debiased_CCDS = Ya_obs_pred - # obs estimate
      Ya_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Ya_rand_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i) predict(fit_rand_overlap,data_Aa_rand[[i]])$predictions) %>% 
      set_colnames(levels_A) 
    
  } else if(complex_fit_models=="bart"){
    Ya_rand_pred = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_rand,data_Aa_rand_matrix[[i]][,colnames(fit_rand$varcount)]) %>%
        apply(.,2,mean)}) %>% 
        set_colnames(levels_A)
    Ya_obs_pred = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_obs,data_Aa_obs_matrix[[i]][,colnames(fit_obs$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSa = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_obs_overlap,data_Aa_obs_matrix[[i]][,colnames(fit_obs_overlap$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_rand_overlap,data_Aa_obs_matrix[[i]][,colnames(fit_rand_overlap$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDS = Ya_obs_pred_bias_CCDSa - Ya_obs_pred_bias_CCDSb  # bias estimate for obs data
    Ya_obs_pred_debiased_CCDS = Ya_obs_pred - # obs estimate
      Ya_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Ya_rand_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_rand_overlap,data_Aa_rand_matrix[[i]][,colnames(fit_rand_overlap$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    
  } else { # note: error for ksvm if some columns of a covariate all have the same values for a given treatment 
    Ya_rand_pred = sapply(1:(length(levels_A)), function(i) predict(fit_rand,data_Aa_rand_matrix[[i]])) %>% 
      set_colnames(levels_A)
    Ya_obs_pred = sapply(1:(length(levels_A)), function(i) predict(fit_obs,data_Aa_obs_matrix[[i]])) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSa = sapply(1:(length(levels_A)), function(i) predict(fit_obs_overlap,data_Aa_obs_matrix[[i]])) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i) predict(fit_rand_overlap,data_Aa_obs_matrix[[i]])) %>% 
      set_colnames(levels_A)
    Ya_obs_pred_bias_CCDS = Ya_obs_pred_bias_CCDSa - Ya_obs_pred_bias_CCDSb  # bias estimate for obs data
    Ya_obs_pred_debiased_CCDS = Ya_obs_pred - # obs estimate
      Ya_obs_pred_bias_CCDS # bias estimate for obs data
    
    # Used in CCDS AIPW
    Ya_rand_pred_bias_CCDSb = sapply(1:(length(levels_A)), function(i) predict(fit_rand_overlap,data_Aa_rand_matrix[[i]])) %>% 
      set_colnames(levels_A) 
  }
  
  # Target population estimates
  Ya_target_pred_CCDS = matrix(NA,n_target, length(levels_A)) %>% 
    set_colnames(levels_A)
  Ya_target_pred_CCDS[S_target==1,] = Ya_rand_pred
  Ya_target_pred_CCDS[S_target==0,] = Ya_obs_pred_debiased_CCDS
  
  # Mean potential outcome estimates
  mean_Ya_target_pred_CCDS = apply(Ya_target_pred_CCDS, 2, weighted.mean, w = weight_sampling)
  
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
  w_bias = ifelse(S_target==1 & overlap_target==1,(1-pi_S_trim)/w_bias_denom,0)
  w_bias_trim = ifelse(S_target==1 & overlap_target==1,(1-pi_S_trim)/w_bias_denom_trim,0)
  
  # Amount of trimming
  n_trim_w_bias = mean(w_bias != w_bias_trim)
  
  # Normalize weights
  w_bias_norm = w_bias_trim/sum(w_bias_trim)
  w_bias_norm_rand_overlap = w_bias_norm[S_target==1 & overlap_target==1]
  
  ### Second-stage projection of bias estimates onto X: predict bias using randomized data
  if(complex_fit_models=="ensemble"){
    Ya_rand_pred_bias_overlap_2stageCCDS_a = sapply(1:(length(levels_A)), function(i) 
      predict(fit_obs_overlap,data_Aa_rand_matrix_overlap[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A) 
    Ya_rand_pred_bias_overlap_2stageCCDS_b = sapply(1:(length(levels_A)), function(i) 
      predict(fit_rand_overlap,data_Aa_rand_matrix_overlap[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A) 
    Ya_rand_pred_bias_overlap_2stageCCDS = Ya_rand_pred_bias_overlap_2stageCCDS_a-
      Ya_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
  } else if(complex_fit_models=="ranger"){
    Ya_rand_pred_bias_overlap_2stageCCDS_a = sapply(1:(length(levels_A)), function(i)
      predict(fit_obs_overlap,data_Aa_rand_overlap[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageCCDS_b = sapply(1:(length(levels_A)), function(i)
      predict(fit_rand_overlap,data_Aa_rand_overlap[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageCCDS = Ya_rand_pred_bias_overlap_2stageCCDS_a-
      Ya_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
  } else if(complex_fit_models=="bart"){
    Ya_rand_pred_bias_overlap_2stageCCDS_a = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_obs_overlap,data_Aa_rand_matrix_overlap[[i]][,colnames(fit_obs_overlap$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageCCDS_b = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_rand_overlap,data_Aa_rand_matrix_overlap[[i]][,colnames(fit_rand_overlap$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
   Ya_rand_pred_bias_overlap_2stageCCDS = Ya_rand_pred_bias_overlap_2stageCCDS_a-
      Ya_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
  } else {
    Ya_rand_pred_bias_overlap_2stageCCDS_a = sapply(1:(length(levels_A)), function(i) predict(fit_obs_overlap,data_Aa_rand_matrix_overlap[[i]])) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageCCDS_b = sapply(1:(length(levels_A)), function(i) predict(fit_rand_overlap,data_Aa_rand_matrix_overlap[[i]])) %>%
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageCCDS = Ya_rand_pred_bias_overlap_2stageCCDS_a-
      Ya_rand_pred_bias_overlap_2stageCCDS_b # bias estimate for obs data
  }

  ## Estimate restricted mean model for bias term E(Ya_rand_pred_bias_overlap_2stageCCDS|S=0,X_overlap,X) within the randomized data in the overlap region
  fit_Ya_rand_overlap_2stageCCDS = as.data.frame(Xmain_rand_overlap_matrix) %>%
    lm(Ya_rand_pred_bias_overlap_2stageCCDS ~ ., data=.)
  
  Ya_obs_pred_bias_2stageCCDS =predict(fit_Ya_rand_overlap_2stageCCDS,
                                       as.data.frame(data_Xmain_obs_matrix)) # bias estimate for obs data
  
  Ya_obs_pred_debiased_2stageCCDS = Ya_obs_pred - # obs estimate
    Ya_obs_pred_bias_2stageCCDS # bias estimate for obs data
  
  ## Weighted version
  fit_Ya_rand_overlap_weighted2stageCCDS = as.data.frame(Xmain_rand_overlap_matrix) %>%
    lm(Ya_rand_pred_bias_overlap_2stageCCDS ~ ., data=., weights=w_bias_norm_rand_overlap)
  
  Ya_obs_pred_bias_weighted2stageCCDS = predict(fit_Ya_rand_overlap_weighted2stageCCDS,
                                                as.data.frame(data_Xmain_obs_matrix)) # bias estimate for obs data
  
  Ya_obs_pred_debiased_weighted2stageCCDS = Ya_obs_pred - # obs estimate
    Ya_obs_pred_bias_weighted2stageCCDS # bias estimate for obs data

  # Target population estimates
  Ya_target_pred_2stageCCDS = Ya_target_pred_weighted2stageCCDS = matrix(NA,n_target, length(levels_A)) %>% 
    set_colnames(levels_A)
  Ya_target_pred_2stageCCDS[S_target==1,] = Ya_rand_pred
  Ya_target_pred_2stageCCDS[S_target==0,] = Ya_obs_pred_debiased_2stageCCDS
  
  # Mean potential outcome estimates
  mean_Ya_target_pred_2stageCCDS = apply(Ya_target_pred_2stageCCDS, 2, weighted.mean, w = weight_sampling)
  
  # Weighted target population estimates
  Ya_target_pred_weighted2stageCCDS[S_target==1,] = Ya_rand_pred
  Ya_target_pred_weighted2stageCCDS[S_target==0,] = Ya_obs_pred_debiased_weighted2stageCCDS
  
  # Mean potential outcome estimates
  mean_Ya_target_pred_weighted2stageCCDS = apply(Ya_target_pred_weighted2stageCCDS, 2, weighted.mean, w = weight_sampling)
  
  ##########################################################################################
  ### Estimate potential outcomes: 2-stage hybrid CCDS and weighted 2-stage hybrid CCDS
  ##########################################################################################
  ### Use model fits from entire rand and entire obs data but estimate bias term via predictions 
  ###   for randomized data in overlap region
  if(complex_fit_models=="ensemble"){
    Ya_rand_pred_bias_overlap_2stageHCCDS_a = sapply(1:(length(levels_A)), function(i)
      predict(fit_obs,data_Aa_rand_matrix_overlap[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS_b = sapply(1:(length(levels_A)), function(i)
      predict(fit_rand,data_Aa_rand_matrix_overlap[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS = Ya_rand_pred_bias_overlap_2stageHCCDS_a-
      Ya_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
    
  } else if(complex_fit_models=="ranger"){
    Ya_rand_pred_bias_overlap_2stageHCCDS_a = sapply(1:(length(levels_A)), function(i)
      predict(fit_obs,data_Aa_rand_overlap[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS_b = sapply(1:(length(levels_A)), function(i)
      predict(fit_rand,data_Aa_rand_overlap[[i]])$predictions) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS = Ya_rand_pred_bias_overlap_2stageHCCDS_a-
      Ya_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
  } else if(complex_fit_models=="bart"){
    Ya_rand_pred_bias_overlap_2stageHCCDS_a = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_obs,data_Aa_rand_matrix_overlap[[i]][,colnames(fit_obs$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS_b = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_rand,data_Aa_rand_matrix_overlap[[i]][,colnames(fit_rand$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS = Ya_rand_pred_bias_overlap_2stageHCCDS_a-
      Ya_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
  } else {
    Ya_rand_pred_bias_overlap_2stageHCCDS_a = sapply(1:(length(levels_A)), function(i) 
      predict(fit_obs,data_Aa_rand_matrix_overlap[[i]])) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS_b = sapply(1:(length(levels_A)), function(i) 
      predict(fit_rand,data_Aa_rand_matrix_overlap[[i]])) %>% 
      set_colnames(levels_A)
    Ya_rand_pred_bias_overlap_2stageHCCDS = Ya_rand_pred_bias_overlap_2stageHCCDS_a-
      Ya_rand_pred_bias_overlap_2stageHCCDS_b # bias estimate for obs data
  }
  
  ## Estimate restricted mean model E(Ya_rand_pred_bias_overlap_2stageHCCDS|S=0,X_overlap,X) within the observational data in the overlap region
  fit_Ya_rand_overlap_2stageHCCDS = as.data.frame(Xmain_rand_overlap_matrix) %>%
    lm(Ya_rand_pred_bias_overlap_2stageHCCDS ~ ., data=.)
  
  Ya_obs_pred_bias_2stageHCCDS =  predict(fit_Ya_rand_overlap_2stageHCCDS,
                                          as.data.frame(data_Xmain_obs_matrix)) # bias estimate for obs data
  
  Ya_obs_pred_debiased_2stageHCCDS = Ya_obs_pred - # obs estimate
    Ya_obs_pred_bias_2stageHCCDS # bias estimate for obs data
  
  # Weighted version
  fit_Ya_rand_overlap_weighted2stageHCCDS = as.data.frame(Xmain_rand_overlap_matrix) %>%
    lm(Ya_rand_pred_bias_overlap_2stageHCCDS ~ ., data=., weights=w_bias_norm_rand_overlap)
  
  Ya_obs_pred_bias_weighted2stageHCCDS = predict(fit_Ya_rand_overlap_weighted2stageHCCDS,
                                                 as.data.frame(data_Xmain_obs_matrix)) # bias estimate for obs data
  
  Ya_obs_pred_debiased_weighted2stageHCCDS = Ya_obs_pred - # obs estimate
    Ya_obs_pred_bias_weighted2stageHCCDS # bias estimate for obs data
  
  
  # Target population estimates
  Ya_target_pred_2stageHCCDS = Ya_target_pred_weighted2stageHCCDS = matrix(NA,n_target, length(levels_A)) %>% 
    set_colnames(levels_A)
  Ya_target_pred_2stageHCCDS[S_target==1,] = Ya_rand_pred
  Ya_target_pred_2stageHCCDS[S_target==0,] = Ya_obs_pred_debiased_2stageHCCDS
  
  # Mean potential outcome estimates
  mean_Ya_target_pred_2stageHCCDS = apply(Ya_target_pred_2stageHCCDS, 2, weighted.mean, w = weight_sampling)
  
  # Weighted target population estimates
  Ya_target_pred_weighted2stageHCCDS[S_target==1,] = Ya_rand_pred
  Ya_target_pred_weighted2stageHCCDS[S_target==0,] = Ya_obs_pred_debiased_weighted2stageHCCDS
  
  # Mean potential outcome estimates
  mean_Ya_target_pred_weighted2stageHCCDS = apply(Ya_target_pred_weighted2stageHCCDS, 2, weighted.mean, w = weight_sampling)
  
  
  
  ##########################################################################################
  ### Estimate potential outcomes: CCDS-IPW 
  ##########################################################################################
  mean_Ya_rand_pred_CCDS_IPW = apply(w1_Aa_trim, 2, function(x) sum(x*Y_target)/sum(x))
  mean_Ya_obs_pred_CCDS_IPW = apply(w2_Aa_trim, 2, function(x) sum(x*Y_target)/sum(x))
  mean_Ya_obs_pred_bias_CCDS_IPWa = apply(w3_Aa_trim, 2, function(x) sum(x*Y_target)/sum(x))
  mean_Ya_obs_pred_bias_CCDS_IPWb = apply(w4_Aa_trim, 2, function(x) sum(x*Y_target)/sum(x))
  mean_Ya_obs_pred_bias_CCDS_IPW = mean_Ya_obs_pred_bias_CCDS_IPWa - mean_Ya_obs_pred_bias_CCDS_IPWb  # bias estimate for obs data
  mean_Ya_obs_pred_CCDS_IPW_debiased = mean_Ya_obs_pred_CCDS_IPW - # obs estimate
    mean_Ya_obs_pred_bias_CCDS_IPW # bias estimate for obs data
  
  # Target population estimates
  mean_Ya_target_pred_CCDS_IPW = weight_rand*mean_Ya_rand_pred_CCDS_IPW+
    weight_obs*mean_Ya_obs_pred_CCDS_IPW_debiased
  
  
  ##########################################################################################
  ### Estimate potential outcomes: CCDS-AIPW 
  ##########################################################################################
  Ya_rand_pred_CCDS_AIPW = sapply(1:(length(levels_A)), function(i) 
    weight_rand*w1_Aa_trim[S_target==1,i]/sum(w1_Aa_trim[,i])*(Y_rand-Ya_rand_pred[,i]) + Ya_rand_pred[,i]) #pre-subsetting to rand data; weights then subset to rand A=1 data
  Ya_obs_pred_CCDS_AIPW =  sapply(1:(length(levels_A)), function(i)
    weight_obs*w2_Aa_trim[S_target==0,i]/sum(w2_Aa_trim[,i])*(Y_obs-Ya_obs_pred[,i]) + Ya_obs_pred[,i]) #pre-subsetting to obs data
  Ya_obs_pred_bias_CCDS_AIPWa = sapply(1:(length(levels_A)), function(i)
    weight_obs*w3_Aa_trim[S_target==0,i]/sum(w3_Aa_trim[,i])*(Y_obs-Ya_obs_pred_bias_CCDSa[,i]) + Ya_obs_pred_bias_CCDSa[,i]) #pre-subsetting to obs data
  Ya_pred_bias_CCDS_AIPWb = sapply(1:(length(levels_A)), function(i)
    ifelse(S_target==1,weight_obs*w4_Aa_trim[S_target==1,i]/sum(w4_Aa_trim[,i])*(Y_rand-Ya_rand_pred_bias_CCDSb[,i]),
           Ya_obs_pred_bias_CCDSb[,i])) #combo of rand (first component) and obs (second component) data
  
  
  # Target population estimates: mean potential outcome estimates
  mean_Ya_target_pred_CCDS_AIPW = sapply(1:(length(levels_A)), function(i) 
    1/n_target*(sum(Ya_rand_pred_CCDS_AIPW[,i])+sum(Ya_obs_pred_CCDS_AIPW[,i])-
                                                sum(Ya_obs_pred_bias_CCDS_AIPWa[,i])+sum(Ya_pred_bias_CCDS_AIPWb[,i]))) %>% 
    set_names(levels_A)
  
  ##########################################################################################
  ### Estimate potential outcomes: obs/rand
  ##########################################################################################
  ### Observational/randomized model predictions for observational/randomized data, respectively (ignoring unmeasured confounding)
  Ya_rand_pred_obs_rand = Ya_rand_pred # rand estimate same as CCDS-OR
  Ya_obs_pred_obs_rand = Ya_obs_pred  # obs estimate same as biased preliminary estimate in CCDS-OR
  
  # Target population estimates
  Ya_target_pred_obs_rand = matrix(NA,n_target, length(levels_A)) %>% 
    set_colnames(levels_A)
  Ya_target_pred_obs_rand[S_target==1,] = Ya_rand_pred_obs_rand
  Ya_target_pred_obs_rand[S_target==0,] = Ya_obs_pred_obs_rand
  
  # Mean potential outcome estimates
  mean_Ya_target_pred_obs_rand = apply(Ya_target_pred_obs_rand, 2, weighted.mean, w = weight_sampling)
  
  ##########################################################################################
  ### Estimate potential outcomes: rand
  ##########################################################################################
  ### Randomized model predictions for all data (ignoring positivity violations)
  # Target population estimates
  if(complex_fit_models=="ensemble"){
    Ya_target_pred_rand = sapply(1:(length(levels_A)), function(i)
      predict(fit_rand,data_Aa_target_matrix[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A) 
  } else if(complex_fit_models=="ranger"){
    Ya_target_pred_rand = sapply(1:(length(levels_A)), function(i)
      predict(fit_rand,data_Aa_target[[i]])$predictions) %>% 
      set_colnames(levels_A) 
  } else if(complex_fit_models=="bart"){
    Ya_target_pred_rand = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_rand,data_Aa_target_matrix[[i]][,colnames(fit_rand$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
  } else {
    Ya_target_pred_rand = sapply(1:(length(levels_A)), function(i)
      predict(fit_rand,data_Aa_target_matrix[[i]])) %>% 
      set_colnames(levels_A) 
  }
  
  
  # Mean potential outcome estimates
  mean_Ya_target_pred_rand = apply(Ya_target_pred_rand, 2, weighted.mean, w = weight_sampling)
  
  
  
  ##########################################################################################
  ### Estimate potential outcomes: 2-stage whole data 
  ##########################################################################################
  # (1) Use S=0 data to create CATE (Y_a) model: fit_obs
  # (2) Predict CATE (Y_a) for S=1 data using model from (1)
  if(complex_fit_models=="ensemble"){
    Ya_rand_pred_from_fit_obs_2stageWD = sapply(1:(length(levels_A)), function(i)
      predict(fit_obs,data_Aa_rand_matrix[[i]], onlySL = TRUE)$pred) %>% 
      set_colnames(levels_A)
    
  } else if(complex_fit_models=="ranger"){
    Ya_rand_pred_from_fit_obs_2stageWD = sapply(1:(length(levels_A)), function(i)
      predict(fit_obs,data_Aa_rand[[i]])$predictions) %>% 
      set_colnames(levels_A)
    
  } else if(complex_fit_models=="bart"){
    Ya_rand_pred_from_fit_obs_2stageWD = sapply(1:(length(levels_A)), function(i){ 
      predict(fit_obs,data_Aa_rand_matrix[[i]][,colnames(fit_obs$varcount)]) %>%
        apply(.,2,mean)}) %>% 
      set_colnames(levels_A)
    
  } else {
    Ya_rand_pred_from_fit_obs_2stageWD = sapply(1:(length(levels_A)), function(i)
      predict(fit_obs,data_Aa_rand_matrix[[i]])) %>% 
      set_colnames(levels_A)
    
  }
  # (3) Use S=1 data to estimate CATE (Y_a) for S=1 
  Ya_rand_pred_2stageWD = Ya_rand_pred
  
  # (4) Estimate bias: eta_pred = (3) - (2)
  Ya_eta_pred_2stageWD = Ya_rand_pred_2stageWD - Ya_rand_pred_from_fit_obs_2stageWD
  
  # (5) Create linear estimator for eta_hat
  fit_Ya_eta_pred_2stageWD = as.data.frame(Xmain_rand_matrix) %>% 
    lm(Ya_eta_pred_2stageWD ~ ., data=.)

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
  fit_Ya_eta_pred_weighted2stageWD = as.data.frame(Xmain_rand_matrix) %>% 
    lm(Ya_eta_pred_2stageWD ~ ., data=., weights=w_bias_wd_norm_rand) 
  
  # (6) Debias CATE (Y_a) prediction for S=0: use (1) to predict for X_obs + eta_hat(X_obs)
  Ya_obs_pred_2stageWD = Ya_obs_pred # biased observational data prediction
  Ya_obs_pred_bias_2stageWD = predict(fit_Ya_eta_pred_2stageWD,as.data.frame(data_Xmain_obs_matrix)) # estimate of bias from linear model in (5) 
  Ya_obs_pred_debiased_2stageWD = Ya_obs_pred_2stageWD + Ya_obs_pred_bias_2stageWD # debiased observational data estimate
  
  # Weighted version
  Ya_obs_pred_weighted2stageWD = Ya_obs_pred # biased observational data prediction
  Ya_obs_pred_bias_weighted2stageWD = predict(fit_Ya_eta_pred_weighted2stageWD,as.data.frame(data_Xmain_obs_matrix)) # estimate of bias from linear model in (5) 
  Ya_obs_pred_debiased_weighted2stageCCDSWD = Ya_obs_pred_weighted2stageWD + Ya_obs_pred_bias_weighted2stageWD # debiased observational data estimate
  
  
  # (7) Create target population estimate by averaging across predictions for all the target sample data
  Ya_target_pred_2stageWD = Ya_target_pred_weighted2stageWD = matrix(NA,n_target, length(levels_A)) %>% 
    set_colnames(levels_A)
  Ya_target_pred_2stageWD[S_target==1,] = Ya_rand_pred_2stageWD
  Ya_target_pred_2stageWD[S_target==0,] = Ya_obs_pred_debiased_2stageWD
  
  mean_Ya_target_pred_2stageWD = apply(Ya_target_pred_2stageWD, 2, weighted.mean, w = weight_sampling)
  
  # Weighted version
  Ya_target_pred_weighted2stageWD[S_target==1,] = Ya_rand_pred_2stageWD
  Ya_target_pred_weighted2stageWD[S_target==0,] = Ya_obs_pred_debiased_weighted2stageCCDSWD
  
  mean_Ya_target_pred_weighted2stageWD = apply(Ya_target_pred_weighted2stageWD, 2, weighted.mean, w = weight_sampling)

  ##########################################################################################
  ##########################################################################################
  ### Estimate study population means
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  ### Outcome regression
  ##########################################################################################
  mean_Ya_rand_pred_rand = apply(Ya_rand_pred, 2, mean)
  mean_Ya_obs_pred_obs = apply(Ya_obs_pred, 2, mean)
  
  ##########################################################################################
  ### TMLE
  ##########################################################################################
  ## Rand
  fit_tmle_rand <- lapply(seq(1:(length(levels_A))), function (i) tmle(Y=Y_rand, A=NULL, W=X_rand_matrix, 
                                                                       Delta=A_rand_all[,i], Q = cbind(0,Ya_rand_pred[,i]), gbound = propensity_trim_threshold, #Q=c(E(Y|A=0,X), E(Y|A=a,X)); First column of values doesn't matter. Symmetric bounds are assumed, pDelta1 <- .bound(pDelta1, c(1, min(gbounds)))
                                                                       pDelta1 = cbind(0,pi_A_rand_trim[S_target == 1,i]))) # First column of values doesn't matter
  names(fit_tmle_rand) = levels_A
  mean_Ya_rand_pred_TMLE = t(sapply(fit_tmle_rand,function(x) x$estimates$EY1$psi))
  
  ## Obs
  fit_tmle_obs <- lapply(seq(1:(length(levels_A))), function (i) tmle(Y=Y_obs, A=NULL, W=X_obs_matrix, 
                                                                       Delta=A_obs_all[,i], Q = cbind(0,Ya_obs_pred[,i]), gbound = propensity_trim_threshold, #Q=c(E(Y|A=0,X), E(Y|A=a,X)); First column of values doesn't matter. Symmetric bounds are assumed, pDelta1 <- .bound(pDelta1, c(1, min(gbounds)))
                                                                       pDelta1 = cbind(0,pi_A_obs_trim[S_target == 0,i]))) # First column of values doesn't matter
  names(fit_tmle_obs) = levels_A
  mean_Ya_obs_pred_TMLE = t(sapply(fit_tmle_obs,function(x) x$estimates$EY1$psi))
  
  ##########################################################################################
  ### AIPW
  ##########################################################################################
  ## Rand
  # Subset to rand study propensities 
  pi_A_rand_norm = pi_A_rand_trim[S_target == 1,] 
  
  # Obtain AIPW estimate
  mean_Ya_rand_pred_AIPW = sapply(1:length(levels_A), function(i) mean(Ya_rand_pred[,i] + (A_rand == levels_A[i])/pi_A_rand_norm[,i]*(Y_rand - Ya_rand_pred[,i]))) %>% 
    set_names(levels_A) # Equiv: mean((A_rand == levels_A[i])*Y_rand/pi_A_rand_norm[,i] - Ya_rand_pred*((A_rand == levels_A[i]) - pi_A_rand_norm[,i])/pi_A_rand_norm[,i])
  
  ## Obs
  # Subset to obs study propensities 
  pi_A_obs_norm = pi_A_obs_trim[S_target == 0,] 
  
  # Obtain AIPW estimate
  mean_Ya_obs_pred_AIPW = sapply(1:length(levels_A), function(i) mean(Ya_obs_pred[,i] + (A_obs == levels_A[i])/pi_A_obs_norm[,i]*(Y_obs - Ya_obs_pred[,i]))) %>% 
    set_names(levels_A) # Equiv: mean((A_obs == levels_A[i])*Y_obs/pi_A_obs_norm[,i] - Ya_obs_pred*((A_obs == levels_A[i]) - pi_A_obs_norm[,i])/pi_A_obs_norm[,i])
  
  
  ##########################################################################################
  ##########################################################################################
  ### Summarize results
  ##########################################################################################
  ##########################################################################################
  ## Evaluate bias in region of overlap and non-overlap
  mean_bias_Ya_obs_observed = sapply(1:(length(levels_A)), function(i) (mean(Ya_obs_pred[,i])-
                                  mean(Y_obs[A_obs==levels_A[i]]))/mean(Y_obs[A_obs==levels_A[i]])) %>% 
    set_names(levels_A)
  mean_bias_Ya_rand_observed = sapply(1:(length(levels_A)), function(i) (mean(Ya_rand_pred[,i])-
                                  mean(Y_rand[A_rand==levels_A[i]]))/mean(Y_rand[A_rand==levels_A[i]])) %>% 
    set_names(levels_A)
  mean_bias_Ya_rand_overlap_observed = sapply(1:(length(levels_A)), function(i) (mean(Ya_rand_pred[overlap_rand==1,i])-
                                  mean(Y_rand[A_rand==levels_A[i] & overlap_rand==1]))/mean(Y_rand[A_rand==levels_A[i] & overlap_rand==1])) %>% 
    set_names(levels_A)
  
  mean_percent_biases = c(mean_bias_Ya_obs_observed, mean_bias_Ya_rand_observed, mean_bias_Ya_rand_overlap_observed)
  names(mean_percent_biases) = c(paste0("mean_bias_obs_observed.", levels_A),
                                 paste0("mean_bias_rand_observed.", levels_A),
                                 paste0("mean_bias_rand_overlap_observed.", levels_A))
  
  # Summarize number of trimmed weights
  n_trim = c(n_trim_w1,
             n_trim_w2,
             n_trim_w3,
             n_trim_w4,
             n_trim_w_bias,
             n_trim_w_bias_wd)
  names(n_trim) = c(paste0("n_trim_w1.", levels_A),
                    paste0("n_trim_w2.", levels_A),
                    paste0("n_trim_w3.", levels_A),
                    paste0("n_trim_w4.", levels_A),
                    "n_trim_w_bias",
                    "n_trim_w_bias_wd")
  
  # Summarize mean trimmed weights
  weight_mean = c(apply(w1, 2, mean),
                  apply(w2, 2, mean),
                  apply(w3, 2, mean),
                  apply(w4, 2, mean),
                  mean(w_bias),
                  mean(w_bias_wd))
  names(weight_mean) = c(paste0("w1.", levels_A),
                         paste0("w2.", levels_A),
                         paste0("w3.", levels_A),
                         paste0("w4.", levels_A),
                         "w_bias",
                         "w_bias_wd")
  
  # Summarize algorithm weights if using ensemble 
  if(complex_fit_models=="ensemble"){
    alg_weights_S = fit_S$coef
    alg_weights_Roverlap_rand = fit_Roverlap_rand$coef
    alg_weights_Roverlap_obs = fit_Roverlap_obs$coef
    alg_weights_A_rand = sapply(fit_A_rand, function(x) x$coef) %>% t() %>% c() %>%
      set_names(apply(expand.grid(SL_library, levels_A), 1, paste, collapse="."))
    alg_weights_A_obs = sapply(fit_A_obs, function(x) x$coef) %>% t() %>% c() %>%
      set_names(apply(expand.grid(SL_library, levels_A), 1, paste, collapse="."))
    alg_weights_A_rand_overlap = sapply(fit_A_rand_overlap, function(x) x$coef) %>% t() %>% c() %>%
      set_names(apply(expand.grid(SL_library, levels_A), 1, paste, collapse="."))
    alg_weights_A_obs_overlap = sapply(fit_A_obs_overlap, function(x) x$coef) %>% t() %>% c() %>%
      set_names(apply(expand.grid(SL_library, levels_A), 1, paste, collapse="."))
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
  
  ## Remove model fits (can be large objects)
  rm(fit_S, fit_Roverlap_rand, fit_Roverlap_obs, fit_A_rand, fit_A_obs, fit_A_rand_overlap, fit_A_obs_overlap,
     fit_obs, fit_obs_overlap, fit_rand, fit_rand_overlap, fit_tmle_obs, fit_tmle_rand)
  
  ### Output
  output = c(mean_Ya_rand_pred_rand,
             mean_Ya_obs_pred_obs,
             mean_Ya_rand_pred_TMLE,
             mean_Ya_obs_pred_TMLE,
             mean_Ya_rand_pred_AIPW,
             mean_Ya_obs_pred_AIPW,
             mean_Ya_target_pred_CCDS,
             mean_Ya_target_pred_2stageCCDS,
             mean_Ya_target_pred_weighted2stageCCDS,
             mean_Ya_target_pred_2stageHCCDS,
             mean_Ya_target_pred_weighted2stageHCCDS,
             mean_Ya_target_pred_CCDS_IPW,
             mean_Ya_target_pred_CCDS_AIPW,
             mean_Ya_target_pred_obs_rand,
             mean_Ya_target_pred_rand,
             mean_Ya_target_pred_2stageWD,
             mean_Ya_target_pred_weighted2stageWD,
             
             mean_overlap_target,
             mean_overlap_obs,
             mean_overlap_rand,
             min_p_overlap, 
             max_p_overlap,
             
             n_trim,
             weight_mean,
             
             mean_percent_biases)
  
  names(output) = c(paste0("mean_Ya_rand_pred_rand.", levels_A),
                    paste0("mean_Ya_obs_pred_obs.", levels_A),
                    paste0("mean_Ya_rand_pred_TMLE.", levels_A),
                    paste0("mean_Ya_obs_pred_TMLE.", levels_A),
                    paste0("mean_Ya_rand_pred_AIPW.", levels_A),
                    paste0("mean_Ya_obs_pred_AIPW.", levels_A),
                    paste0("mean_target_pred_CCDS.", levels_A),
                    paste0("mean_target_pred_2stageCCDS.", levels_A),
                    paste0("mean_target_pred_weighted2stageCCDS.", levels_A),
                    paste0("mean_target_pred_2stageHCCDS.", levels_A),
                    paste0("mean_target_pred_weighted2stageHCCDS.", levels_A),
                    paste0("mean_target_pred_CCDS_IPW.", levels_A),
                    paste0("mean_target_pred_CCDS_AIPW.", levels_A),
                    paste0("mean_target_pred_obs_rand.", levels_A),
                    paste0("mean_target_pred_rand.", levels_A),
                    paste0("mean_target_pred_2stageWD.", levels_A),
                    paste0("mean_target_pred_weighted2stageWD.", levels_A),
                    
                    "mean_overlap_target",
                    "mean_overlap_obs",
                    "mean_overlap_rand",
                    "min_p_overlap", 
                    "max_p_overlap",
                    
                    names(n_trim),
                    names(weight_mean),
                    names(mean_percent_biases))
  
  if(complex_fit_models %in% c("ensemble")){
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
  
  return(t(output))
}

#####################################################################################################
#####################################################################################################
## Function to summarize results via table: point estimates from get_estimates and bootstrap CI from 
##    get_estimates_with_bootstrap function
#####################################################################################################
#####################################################################################################
get_summary_table = function(results, 
                             boot_results,
                             levels_A,
                             ensemble = F, 
                             n_estimators = 17){ # number of estimators in get_estimates() function output
  
  # Estimators that will be summarized
  my_cols_Ya = c("mean_Ya_rand_pred_rand", "mean_Ya_obs_pred_obs",
                 "mean_Ya_rand_pred_AIPW", "mean_Ya_obs_pred_AIPW",
                 "mean_Ya_rand_pred_TMLE", "mean_Ya_obs_pred_TMLE",
                 "mean_target_pred_rand",
                 "mean_target_pred_obs_rand",
                 "mean_target_pred_CCDS",
                 "mean_target_pred_CCDS_AIPW",
                 "mean_target_pred_CCDS_IPW",
                 #"mean_target_pred_2stageCCDS",
                 "mean_target_pred_weighted2stageCCDS")#,
                 #"mean_target_pred_2stageHCCDS",
                 #"mean_target_pred_weighted2stageHCCDS",
                 #"mean_target_pred2",
                 #"mean_target_pred3",
                 #"mean_target_pred_2stageWD",
                 #"mean_target_pred_weighted2stageWD")
  
  
  Estimators = factor(c("STSM rand OR", "STSM obs OR", # OR = outcome regression
                        "STSM rand AIPW", "STSM obs AIPW", 
                        "STSM rand TMLE", "STSM obs TMLE",
                        "Rand OR","Obs/rand OR",
                       "CCDS","CCDS-AIPW","CCDS-IPW",#"2-stage CCDS",
                       "2-stage CCDS"),#,
                       #"weighted 2-stage hybrid CCDS",
                       #"2-stage WD", 
                       #"weighted 2-stage WD"), 
                     levels=c("STSM rand OR", "STSM rand AIPW", "STSM rand TMLE",
                              "STSM obs OR", "STSM obs AIPW", "STSM obs TMLE",
                              "Rand OR","Obs/rand OR","CCDS","CCDS-AIPW","CCDS-IPW",#"2-stage CCDS",
                              "2-stage CCDS"))#,
                              #"weighted 2-stage hybrid CCDS",
                              #"2-stage WD",
                              #"weighted 2-stage WD"))
  
  ##########################################################################################
  ### Obtain confidence intervals and SE
  ##########################################################################################
  bounds = apply(boot_results[,1:(n_estimators*length(levels_A))],2, # number of estimators * number of levels of A
                 function(x) quantile(x,
                                      probs = c(0.025,0.975),
                                      na.rm = TRUE))
  
  # Adjust CI for multiple testing
  bounds_adj = apply(boot_results[,1:(n_estimators*length(levels_A))],2, # number of estimators * number of levels of A
                 function(x) quantile(x,
                                      probs = c(0.05/(2*length(levels_A)),1-0.05/(2*length(levels_A))),
                                      na.rm = TRUE))

  boot_se = apply(boot_results[,1:(n_estimators*length(levels_A))],2, # number of estimators * number of levels of A
                  function(x) sqrt(var(x, na.rm = TRUE)))

  ci_width = bounds[2,]-bounds[1,]
  
  ci_width_adj = bounds_adj[2,]-bounds_adj[1,]

  
  # Convert to table format
  estimates = results[,1:(n_estimators*length(levels_A))]
  estimate_table = cbind.data.frame(matrix(unlist(strsplit(names(estimates), "[.]")),ncol=2,byrow=T), estimates) %>% 
    set_colnames(c("Estimator","Intervention","Estimate"))
  
  table = cbind(estimate_table, 
                "2.5% Adj" = bounds_adj[1,],
                "97.5% Adj" = bounds_adj[2,],
                "CI Width Adj" = ci_width_adj,
                "2.5% Unadj" = bounds[1,],
                "97.5% Unadj" = bounds[2,],
                "CI Width Unadj" = ci_width,
                "SE Unadj" = boot_se)
  
  # Subset estimators
  table = table %>% filter(Estimator %in% my_cols_Ya) %>% mutate(Estimator = plyr::mapvalues(Estimator,my_cols_Ya, as.character(Estimators)) %>% 
                                                                   factor(.,levels=levels(Estimators)),
                                                                 Intervention = factor(Intervention, levels=levels_A)) %>% 
    arrange(Intervention, Estimator)
  
  
  return(table)
}

#####################################################################################################
#####################################################################################################
## Function to plot results, using as input the table created via get_summary_table()
#####################################################################################################
#####################################################################################################
plot_results = function(my_table, A_levels, my_ylim=NULL){
  # Crop error bars at my_ylim so they show up
  if(!is.null(my_ylim)){
    my_table[,"2.5% Adj"] = ifelse(my_table[,"2.5% Adj"]<my_ylim[1],my_ylim[1],my_table[,"2.5% Adj"])
    my_table[,"97.5% Adj"] = ifelse(my_table[,"97.5% Adj"]>my_ylim[2],my_ylim[2],my_table[,"97.5% Adj"])
  }
  # Line for STSM rand estimate
  STSM_mean <- my_table %>%
    filter(Estimator %in% "STSM rand AIPW") %>% 
    group_by(Intervention) %>%
    summarise(STSM = mean(Estimate))
  
  # Add color grouping for estimators; subset to one set of STSM estimates
  my_table = my_table %>%  filter(!(Estimator %in% c("STSM rand OR", "STSM obs OR",
                                                     "STSM rand TMLE", "STSM obs TMLE"))) %>% 
    mutate(Group = ifelse(Estimator %in% c("STSM rand AIPW", "STSM obs AIPW"), "STSM",
                       ifelse(Estimator %in% c("Rand OR", "Obs/rand OR"), "PTSM comparison",
                              "PTSM novel")),
           Estimator = recode(Estimator, "CCDS"="CCDS-OR", 
                              "STSM rand AIPW" = "Rand AIPW",
                              "STSM obs AIPW" = "Obs AIPW",
                              "Rand OR" = "Rand",
                              "Obs/rand OR" = "Obs/rand"),
           Group = factor(Group, levels = c("STSM", "PTSM comparison", "PTSM novel")))  

  
  # Plot results
  dodge <- position_dodge(width=0.25)  
  my_plot = my_table %>% ggplot(aes(x=Estimator, y = Estimate, col=Group)) +  
    geom_hline(data=STSM_mean, aes(yintercept = STSM)) +
    geom_point(size=1) +
    geom_pointrange(aes(x=Estimator, ymin=`2.5% Adj`, ymax=`97.5% Adj`))+
    theme_bw() + 
    xlab("") + ylab("log spending")+
    facet_wrap(~ factor(Intervention, A_levels)) +
    theme(axis.text.x = element_text(angle=77, hjust=1))+
    theme(legend.position="bottom", legend.title = element_blank()) +
    scale_color_manual(values=(c("#453009", "#80cdc1", "#048b77"))) 
  
  if(!is.null(my_ylim)){
    my_plot = my_plot + ylim(my_ylim)
  }
  return(my_plot)
}
