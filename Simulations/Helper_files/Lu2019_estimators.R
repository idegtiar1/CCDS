#####################################################################################################
#####################################################################################################
## Function with Lu et al. 2019 comparison estimators - Causal Inference for Comprehensive Cohort Studies (arxiv)
## Replaces gam regressions with linear or ensemble regressions to match fit choices for other estimators
## Only outputs PTSM estimates (mu) not STSM estimates (nu)
#####################################################################################################
#####################################################################################################
# Inputs are duplicative, for ease of use
gen_estimates_lu2019_match_fit = function(datc = target_sample, # data sample data
                                          X_target = X_target, # covariate dataframe for target sample
                                          X_target_matrix = X_target_matrix, # matrix version of covariate dataframe for target sample
                                          Y_target = Y_target, # outcome vector for target sample
                                          S_target = S_target, # randomized group status vector for target sample
                                          A_target = A_target, # treatment status vector for target sample
                                          X1_high_target = X1_high_target, # X1 threshold above which S=1 for all observations
                                          X1_low_target = X1_low_target, # X1 threshold below which S=0 for all observations
                                          complex_fit_models = complex_fit_models, # whether to use ensemble ("ensemble"), random forest ("ranger"), ksvm ("ksvm"), or correctly specified models ("correctly_specified") to fit models for Y, A, and S
                                          SL_library = SL_library, # specifies SuperLearner library for complex_fit_models="ensemble"
                                K=5) # number of folds (sample splits)
{
  A_target_01 = ifelse(A_target == 2, 1, 0) # re-label A from 1 & 2 to 0 & 1
  n = nrow(datc) # target sample size
  YA_target_01_pred_a1a2 = matrix(0,K,4) # 4 = Y1hat, Y0hat, se(Y1hat), se(Y0hat)
  YA_target_01_pred_a1a3 = matrix(0,K,4)
  YA_target_01_pred_a1a2a3 = matrix(0,K,4)
  Ya_study_pred_a1 = matrix(0,K,6)
  Ya_study_pred_a1a2a3 = matrix(0,K,6)
  sam = sample(1:K,n,replace=TRUE)
  
  for (k in 1: K) {
    datmk = datc[sam!=k,] # minus kth split
    datk = datc[sam==k,] # kth split
    
    # Minus kth split
    X_target_mk = X_target[sam!=k,]
    X_target_matrix_mk = as.data.frame(X_target_matrix[sam!=k,])
    Y_target_mk = Y_target[sam!=k]
    S_target_mk = S_target[sam!=k]
    A_target_01_mk = A_target_01[sam!=k]

    X1_high_target_mk = X1_high_target[sam!=k]
    X1_low_target_mk = X1_low_target[sam!=k]
    
    # kth split
    X_target_k = X_target[sam==k,]
    X_target_matrix_k = as.data.frame(X_target_matrix[sam==k,])
    Y_target_k = Y_target[sam==k]
    S_target_k = S_target[sam==k]
    A_target_01_k = A_target_01[sam==k]
    
    X1_high_target_k = X1_high_target[sam==k]
    X1_low_target_k = X1_low_target[sam==k]
    
    ## Fit regressions
    if(complex_fit_models=="ensemble"){
      # model for Y for A=1 patients
      # predicted probabilties on minus kth split and kth split
      fit_treated_mk <- SuperLearner(Y=Y_target_mk[A_target_01_mk==1], 
                                         X=X_target_matrix_mk[A_target_01_mk==1,], 
                               family="gaussian", SL.library=SL_library, verbose=F)
      pred_treated_k = predict(fit_treated_mk,X_target_matrix_k, onlySL = TRUE)$pred 

      
      # model for Y for A=0 patients
      # predicted probabilties on minus kth split and kth split
      fit_control_mk <- SuperLearner(Y=Y_target_mk[A_target_01_mk==0], 
                                         X=X_target_matrix_mk[A_target_01_mk==0,], 
                                         family="gaussian", SL.library=SL_library, verbose=F)
      pred_control_k = predict(fit_control_mk,X_target_matrix_k, onlySL = TRUE)$pred 
      
      
      # model for Y for A=1 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_treated_mk <- SuperLearner(Y=Y_target_mk[S_target_mk==1 & A_target_01_mk==1], 
                                         X=X_target_matrix_mk[S_target_mk==1 & A_target_01_mk==1,], 
                                         family="gaussian", SL.library=SL_library, verbose=F)
      pred_rand_treated_k = predict(fit_rand_treated_mk,X_target_matrix_k, onlySL = TRUE)$pred 
      
      # model for Y for A=0 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_control_mk <- SuperLearner(Y=Y_target_mk[S_target_mk==1 & A_target_01_mk==0], 
                                              X=X_target_matrix_mk[S_target_mk==1 & A_target_01_mk==0,], 
                                              family="gaussian", SL.library=SL_library, verbose=F)
      pred_rand_control_k = predict(fit_rand_control_mk,X_target_matrix_k, onlySL = TRUE)$pred 
      
      # model for enrollment into RCT
      # predicted probabilties on minus kth split and kth split
      fit_S_mk = SuperLearner(Y=S_target_mk, X=as.data.frame(X_target_matrix_mk), 
                                  family="binomial", SL.library=SL_library, verbose=F) 
      # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
      pred_S_k = predict(fit_S_mk, as.data.frame(X_target_matrix_k), onlySL = TRUE)$pred #as.data.frame(X_target_matrix)
      
      # model for A in OBS
      # predicted probabilties on minus kth split and kth split
      fit_A_obs_mk = SuperLearner(Y=A_target_01_mk[S_target_mk==0], 
                                      X=as.data.frame(X_target_matrix_mk[S_target_mk==0,]), # Note: only works for binary A; need to create wrappers for multinomial classification: https://github.com/ecpolley/SuperLearner/issues/16
                                family="binomial", SL.library=SL_library, verbose=F)
      pred_A_obs_k = predict(fit_A_obs_mk, as.data.frame(X_target_matrix_k), onlySL = TRUE)$pred #
      
    } else if(complex_fit_models=="ksvm"){
      # model for Y for A=1 patients
      # predicted probabilties on minus kth split and kth split
      fit_treated_mk <- cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==1], 
                                             X_target_mk[A_target_01_mk==1,]) %>% 
        ksvm(Y ~ .,data = .)
      pred_treated_k = predict(fit_treated_mk, X_target_k)
      
      # model for Y for A=0 patients
      # predicted probabilties on minus kth split and kth split
      fit_control_mk <- cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==0], 
                                             X_target_mk[A_target_01_mk==0,]) %>% 
        ksvm(Y ~ .,data = .)
      pred_control_k = predict(fit_control_mk, X_target_k)
      
      # model for Y for A=1 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_treated_mk <- cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==1], 
                                                  X_target_mk[S_target_mk==1 & A_target_01_mk==1,]) %>% 
        ksvm(Y ~ .,data = .)
      pred_rand_treated_k = predict(fit_rand_treated_mk, X_target_k)
      
      # model for Y for A=0 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_control_mk <- cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==0], 
                                                  X_target_mk[S_target_mk==1 & A_target_01_mk==0,]) %>% 
        ksvm(Y ~ .,data = .)
      pred_rand_control_k = predict(fit_rand_control_mk, X_target_k)
      
      # model for enrollment into RCT
      # predicted probabilties on minus kth split and kth split
      fit_S_mk = cbind.data.frame(S_target_mk,X_target_matrix_mk) %>% 
        ksvm(formula(.), data=.,type="C-svc",prob.model=T)
      # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
      pred_S_k = predict(fit_S_mk, as.data.frame(X_target_matrix_k),type="probabilities")[,2]
      
      # model for A in OBS
      # predicted probabilties on minus kth split and kth split
      fit_A_obs_mk = cbind.data.frame("A"=A_target_01_mk[S_target_mk==0],
                                          X_target_matrix_mk[S_target_mk==0,]) %>% 
        ksvm(formula(.), data=.,type="C-svc",prob.model=T)
      pred_A_obs_k = predict(fit_A_obs_mk, as.data.frame(X_target_matrix_k),type="probabilities")[,2]
      
    } else if(complex_fit_models=="ranger"){
      # model for Y for A=1 patients
      # predicted probabilties on minus kth split and kth split
      fit_treated_mk <- cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==1], 
                                             X_target_mk[A_target_01_mk==1,]) %>%
        ranger(Y ~ .,data = .,
               num.trees = 300, min.node.size=floor(nrow(X_target_mk[A_target_01_mk==1,])*0.05),
               num.threads=1, respect.unordered.factors=T)
      pred_treated_k = predict(fit_treated_mk,X_target_k)$predictions
      
      # model for Y for A=0 patients
      # predicted probabilties on minus kth split and kth split
      fit_control_mk <- cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==0], 
                                             X_target_mk[A_target_01_mk==0,]) %>%
        ranger(Y ~ .,data = .,
               num.trees = 300, min.node.size=floor(nrow(X_target_mk[A_target_01_mk==0,])*0.05),
               num.threads=1, respect.unordered.factors=T)
      pred_control_k = predict(fit_control_mk,X_target_k)$predictions
      
      # model for Y for A=1 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_treated_mk <- cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==1], 
                                             X_target_mk[S_target_mk==1 & A_target_01_mk==1,]) %>%
        ranger(Y ~ .,data = .,
               num.trees = 300, min.node.size=floor(nrow(X_target_mk[S_target_mk==1 & A_target_01_mk==1,])*0.05),
               num.threads=1, respect.unordered.factors=T)
      pred_rand_treated_k = predict(fit_rand_treated_mk,X_target_k)$predictions

      # model for Y for A=0 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_control_mk <- cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==0], 
                                                  X_target_mk[S_target_mk==1 & A_target_01_mk==0,]) %>%
        ranger(Y ~ .,data = .,
               num.trees = 300, min.node.size=floor(nrow(X_target_mk[S_target_mk==1 & A_target_01_mk==0,])*0.05),
               num.threads=1, respect.unordered.factors=T)
      pred_rand_control_k = predict(fit_rand_control_mk,X_target_k)$predictions
      
      
      # model for enrollment into RCT
      # predicted probabilties on minus kth split and kth split
      fit_S_mk = cbind.data.frame(S_target_mk,X_target_matrix_mk) %>% 
        ranger(formula(.), data=.,
               num.trees = 300, min.node.size=floor(nrow(X_target_matrix_mk)*0.05), # these hyperparameters work reasonably in tests
               num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
      # Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred - unstable because so many values close to zero
      pred_S_k = cbind.data.frame(S_target_k,X_target_matrix_k) %>%
        predict(fit_S_mk, .,type = "response") %>% .$predictions %>% .[,2]
      
      # model for A in OBS
      # predicted probabilties on minus kth split and kth split
      fit_A_obs_mk = cbind.data.frame("A"=A_target_01_mk[S_target_mk==0],X_target_matrix_mk[S_target_mk==0,]) %>% 
        ranger(formula(.), data=.,
               num.trees = 300, min.node.size=floor(nrow(X_target_matrix_mk[S_target_mk==0,])*0.05), # these hyperparameters work reasonably in tests
               num.threads=1, classification=T,probability=T, respect.unordered.factors=T)
      pred_A_obs_k = cbind.data.frame(X_target_matrix_k) %>%
        predict(fit_A_obs_mk, .,type = "response") %>% .$predictions %>% .[,2]

      
    } else if(complex_fit_models=="correctly_specified"){
      # model for Y for A=1 patients
      # predicted probabilties on minus kth split and kth split
      fit_treated_mk = cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==1], 
                                            X_target_matrix_mk[A_target_01_mk==1,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_treated_k = predict(fit_treated_mk,X_target_matrix_k)
      
      # model for Y for A=0 patients
      # predicted probabilties on minus kth split and kth split
      fit_control_mk = cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==0], 
                                            X_target_matrix_mk[A_target_01_mk==0,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_control_k = predict(fit_control_mk,X_target_matrix_k)
      
      # model for Y for A=1 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_treated_mk = cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==1], 
                                            X_target_matrix_mk[S_target_mk==1 & A_target_01_mk==1,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_rand_treated_k = predict(fit_rand_treated_mk,X_target_matrix_k)

      # model for Y for A=0 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_control_mk = cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==0], 
                                                 X_target_matrix_mk[S_target_mk==1 & A_target_01_mk==0,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_rand_control_k = predict(fit_rand_control_mk,X_target_matrix_k)
      
      
      # model for enrollment into RCT
      # predicted probabilties on minus kth split and kth split
      fit_S_mk = cbind.data.frame(S_target_mk,X_target_matrix_mk,
                                      "X1_high"=X1_high_target_mk,
                                      "X1_low"=X1_low_target_mk) %>% 
        glm(formula(.), data=.,family="binomial")
      pred_S_k = cbind.data.frame(X_target_matrix_k,
                                      "X1_high"=X1_high_target_k,
                                      "X1_low"=X1_low_target_k) %>%
        predict(fit_S_mk, newdata = ., type = "response")
      
      # model for A in OBS
      # predicted probabilties on minus kth split and kth split
      fit_A_obs_mk =  cbind.data.frame("A"=A_target_01_mk[S_target_mk==0],
                                           X_target_matrix_mk[S_target_mk==0,]) %>% 
        glm(formula(.), data=.,family="binomial")
      pred_A_obs_k = predict(fit_A_obs_mk,as.data.frame(X_target_matrix_k),type = "response")
      

    } else { # NOTE: not using S_ridge==T scenario since I didn't end up using it throughout simulations
      # model for Y for A=1 patients
      # predicted probabilties on minus kth split and kth split
      # model for Y for A=1 patients
      # predicted probabilties on minus kth split and kth split
      fit_treated_mk = cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==1], 
                                            X_target_matrix_mk[A_target_01_mk==1,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_treated_k = predict(fit_treated_mk,X_target_matrix_k)
      
      # model for Y for A=0 patients
      # predicted probabilties on minus kth split and kth split
      fit_control_mk = cbind.data.frame("Y"=Y_target_mk[A_target_01_mk==0], 
                                            X_target_matrix_mk[A_target_01_mk==0,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_control_k = predict(fit_control_mk,X_target_matrix_k)
      
      # model for Y for A=1 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_treated_mk = cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==1], 
                                                 X_target_matrix_mk[S_target_mk==1 & A_target_01_mk==1,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_rand_treated_k = predict(fit_rand_treated_mk,X_target_matrix_k)
      
      # model for Y for A=0 patients in RCT
      # predicted probabilties on minus kth split and kth split
      fit_rand_control_mk = cbind.data.frame("Y"=Y_target_mk[S_target_mk==1 & A_target_01_mk==0], 
                                                 X_target_matrix_mk[S_target_mk==1 & A_target_01_mk==0,]) %>% 
        lm(Y ~ ., data=.) # separate regressions already interact X1 and X1_cube with A
      pred_rand_control_k = predict(fit_rand_control_mk,X_target_matrix_k)
      
      
      # model for enrollment into RCT
      # predicted probabilties on minus kth split and kth split
      fit_S_mk = cbind.data.frame(S_target_mk,X_target_matrix_mk) %>% glm(formula(.), data=.,family="binomial")
      pred_S_k = predict(fit_S_mk, X_target_matrix_k, type = "response")

      # model for A in OBS
      # predicted probabilties on minus kth split and kth split
      fit_A_obs_mk =  cbind.data.frame("A"=A_target_01_mk[S_target_mk==0],
                                           X_target_matrix_mk[S_target_mk==0,]) %>% 
        glm(formula(.), data=.,family="binomial")
      pred_A_obs_k = predict(fit_A_obs_mk,as.data.frame(X_target_matrix_k),type = "response")
      
    }
    
    # estimated conditional probability of Y under A=1 (A=0) given covariates
    # estimated on minus kth split and kth split
    Y1_pred_k = pred_treated_k
    Y0_pred_k = pred_control_k
    
    # probability of A in RCT for minus kth split and kth split
    pi_A_S1_mk = 0.6
    pi_A_S1_k = 0.6
    
    
    # estimated conditional probability of A=1 (A=0) given covariates
    # estimated on minus kth split and kth split
    pi_A1_k = pred_S_k * pi_A_S1_mk + (1-pred_S_k) * pred_A_obs_k
    pi_A0_k = pred_S_k * (1-pi_A_S1_mk) + (1-pred_S_k) * (1-pred_A_obs_k)
    
    # kth contribution to estimate of mu under Assumptions treated, A2
    # kth contribution to standard error calculation
    YA_target_01_pred_a1a2[k,1] = mean((datk$A==1)*(datk$Y)/pi_A1_k +
                                          (1-(datk$A==1)/pi_A1_k)*Y1_pred_k) # Y1, resmua1a2
    YA_target_01_pred_a1a2[k,2] = mean((datk$A==0)*(datk$Y)/pi_A0_k +
                                          (1-(datk$A==0)/pi_A0_k)*Y0_pred_k) # Y0
    YA_target_01_pred_a1a2[k,3] = sum(((datk$A==1)*(datk$Y)/pi_A1_k +
                                          (1-(datk$A==1)/pi_A1_k)*Y1_pred_k-YA_target_01_pred_a1a2[k,1])^2) # var(Y1) 
    YA_target_01_pred_a1a2[k,4] = sum(((datk$A==0)*(datk$Y)/pi_A0_k +
                                          (1-(datk$A==0)/pi_A0_k)*Y0_pred_k-YA_target_01_pred_a1a2[k,2])^2) #var(Y0)
    
    # kth contribution to estimate of mu under Assumptions treated, A3
    # kth contribution to standard error calculation
    YA_target_01_pred_a1a3[k,1] = mean((datk$A==1)*(datk$S)*(datk$Y)/(pred_S_k*pi_A_S1_mk) +
                                          (1-(datk$A==1)*(datk$S)/(pred_S_k*pi_A_S1_mk))*pred_rand_treated_k) # YA_target_01_pred_a1a3
    YA_target_01_pred_a1a3[k,2] = mean((datk$A==0)*(datk$S)*(datk$Y)/(pred_S_k*(1-pi_A_S1_mk)) +
                                          (1-(datk$A==0)*(datk$S)/(pred_S_k*(1-pi_A_S1_mk)))*pred_rand_control_k)
    YA_target_01_pred_a1a3[k,3] = sum(((datk$A==1)*(datk$S)*(datk$Y)/(pred_S_k*pi_A_S1_mk) +
                                          (1-(datk$A==1)*(datk$S)/(pred_S_k*pi_A_S1_mk))*pred_rand_treated_k-YA_target_01_pred_a1a3[k,1])^2)
    YA_target_01_pred_a1a3[k,4] = sum(((datk$A==0)*(datk$S)*(datk$Y)/(pred_S_k*(1-pi_A_S1_mk)) +
                                          (1-(datk$A==0)*(datk$S)/(pred_S_k*(1-pi_A_S1_mk)))*pred_rand_control_k-YA_target_01_pred_a1a3[k,2])^2)
    
    # kth contribution to estimate of mu under Assumptions treated, A2, A3
    # kth contribution to standard error calculation
    YA_target_01_pred_a1a2a3[k,1] = mean((datk$A==1)*(datk$Y)/pi_A1_k +
                                            (1-(datk$A==1)/pi_A1_k)*Y1_pred_k)
    YA_target_01_pred_a1a2a3[k,2] = mean((datk$A==0)*(datk$Y)/pi_A0_k +
                                            (1-(datk$A==0)/pi_A0_k)*Y0_pred_k)
    YA_target_01_pred_a1a2a3[k,3] = sum(((datk$A==1)*(datk$Y)/pi_A1_k +
                                            (1-(datk$A==1)/pi_A1_k)*Y1_pred_k-YA_target_01_pred_a1a2a3[k,1])^2)
    YA_target_01_pred_a1a2a3[k,4] = sum(((datk$A==0)*(datk$Y)/pi_A0_k +
                                            (1-(datk$A==0)/pi_A0_k)*Y0_pred_k-YA_target_01_pred_a1a2a3[k,2])^2)
    
    
  }
  # Estimates, standard errors and confidence intervals of mu
  # mu = target population estimates
  mua1a2 = apply(YA_target_01_pred_a1a2,2,mean)[1:2] # mean across k folds
  semua1a2 = sqrt((apply(YA_target_01_pred_a1a2,2,sum)/n^2)[3:4]) # sum across k folds
  lcimua1a2 = mua1a2 - 1.96*semua1a2
  ucimua1a2 = mua1a2 + 1.96*semua1a2
  mua1a3 = apply(YA_target_01_pred_a1a3,2,mean)[1:2]
  semua1a3 = sqrt((apply(YA_target_01_pred_a1a3,2,sum)/n^2)[3:4])
  lcimua1a3 = mua1a3 - 1.96*semua1a3
  ucimua1a3 = mua1a3 + 1.96*semua1a3
  mua1a2a3 = apply(YA_target_01_pred_a1a2a3,2,mean)[1:2]
  semua1a2a3 = sqrt((apply(YA_target_01_pred_a1a2a3,2,sum)/n^2)[3:4])
  lcimua1a2a3 = mua1a2a3 - 1.96*semua1a2a3
  ucimua1a2a3 = mua1a2a3+ 1.96*semua1a2a3
  
  
  # Estimates, standard errors and confidence intervals of treatment effects
  deltamua1a2 = mua1a2[2]-mua1a2[1] # Y0 - Y1 (i.e., Y1 - Y2)
  sedeltamua1a2 = sqrt(sum(semua1a2^2))
  cideltamua1a2 =c(deltamua1a2-1.96*sedeltamua1a2,
                   deltamua1a2+1.96*sedeltamua1a2)
  deltamua1a3 = mua1a3[2]-mua1a3[1]
  sedeltamua1a3 = sqrt(sum(semua1a3^2))
  cideltamua1a3 =c(deltamua1a3-1.96*sedeltamua1a3,
                   deltamua1a3+1.96*sedeltamua1a3)
  deltamua1a2a3 = mua1a2a3[2]-mua1a2a3[1]
  sedeltamua1a2a3 = sqrt(sum(semua1a2a3^2))
  cideltamua1a2a3 =c(deltamua1a2a3-1.96*sedeltamua1a2a3,
                     deltamua1a2a3+1.96*sedeltamua1a2a3)
  
  out = c(
    # Target sample estimates
    "mean_Y1_target_pred_a1a2" = mua1a2[2],
    "se_Y1_target_pred_a1a2" = semua1a2[2],
    "lci_Y1_target_pred_a1a2" = lcimua1a2[2],
    "uci_Y1_target_pred_a1a2" = ucimua1a2[2],
    "mean_Y2_target_pred_a1a2" = mua1a2[1],
    "se_Y2_target_pred_a1a2" = semua1a2[1], 
    "lci_Y2_target_pred_a1a2" = lcimua1a2[1],
    "uci_Y2_target_pred_a1a2" = ucimua1a2[1],
    
    "mean_Y1_target_pred_a1a3" = mua1a3[2],
    "se_Y1_target_pred_a1a3" = semua1a3[2],
    "lci_Y1_target_pred_a1a3" = lcimua1a3[2],
    "uci_Y1_target_pred_a1a3" = ucimua1a3[2],
    "mean_Y2_target_pred_a1a3" = mua1a3[1],
    "se_Y2_target_pred_a1a3" = semua1a3[1], 
    "lci_Y2_target_pred_a1a3" = lcimua1a3[1],
    "uci_Y2_target_pred_a1a3" = ucimua1a3[1],
    
    "mean_Y1_target_pred_a1a2a3" = mua1a2a3[2],
    "se_Y1_target_pred_a1a2a3" = semua1a2a3[2],
    "lci_Y1_target_pred_a1a2a3" = lcimua1a2a3[2],
    "uci_Y2_target_pred_a1a2a3" = ucimua1a2a3[1],
    "mean_Y2_target_pred_a1a2a3" = mua1a2a3[1],
    "se_Y2_target_pred_a1a2a3" = semua1a2a3[1], 
    "lci_Y2_target_pred_a1a2a3" = lcimua1a2a3[1],
    "uci_Y1_target_pred_a1a2a3" = ucimua1a2a3[2],
    
    
    "mean_Y1minusY2_target_pred_a1a2" = deltamua1a2,
    "se_Y1minusY2_target_pred_a1a2" = sedeltamua1a2,
    "lci_Y1minusY2_target_pred_a1a2" = cideltamua1a2[1],
    "uci_Y1minusY2_target_pred_a1a2" = cideltamua1a2[2],
    
    "mean_Y1minusY2_target_pred_a1a3" = deltamua1a3,
    "se_Y1minusY2_target_pred_a1a3" = sedeltamua1a3,
    "lci_Y1minusY2_target_pred_a1a3" = cideltamua1a3[1],
    "uci_Y1minusY2_target_pred_a1a3" = cideltamua1a3[2],
    
    "mean_Y1minusY2_target_pred_a1a2a3" = deltamua1a2a3,
    "se_Y1minusY2_target_pred_a1a2a3" = sedeltamua1a2a3,
    "lci_Y1minusY2_target_pred_a1a2a3" = cideltamua1a2a3[1],
    "uci_Y1minusY2_target_pred_a1a2a3" = cideltamua1a2a3[2]#,
    
    
  )
  
}

