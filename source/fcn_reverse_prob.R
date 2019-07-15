# This function calculates the reverse probability (favoring H1 at interim analysis but fail to reject at decision analysis, or vice versa).


library(ltmle)
library(mvtnorm)

source("source/fcn_err_sp_bdry.R")
source("source/fcn_compcov.R")
source("source/fcn_trialresults.R")

reverse_prob <- function(v, parallel_filenames_H0, parallel_filenames_H1){
  # v needs to be a list of: nsim, K, alpha, beta, rho, f_err, g_err
  #   note: nsim should be the total number, i.e. add all parallel jobs up
  # caseH0, caseH1 are the file names of the stored result, e.g. "est_H0" and "est_H1"
  # idx is the index of parallel result files, currently are 1:50. (used in compcov() function)
  
  
  ## output: vector of length 12,
  ## each entry is the reversal probability under a scenario.
  ## "good" indicates wrong interim analysis but correct decision analysis.
  ## "bad" indicates correct interim analysis but incorrect decision analysis (want to avoid this!).
  output <- data.frame(ltmle_H0 = 999, ltmle_H0_good = 999, ltmle_H0_bad = 999,
                       ltmle_H1 = 999, ltmle_H1_good = 999, ltmle_H1_bad = 999,
                       unadj_H0 = 999, unadj_H0_good = 999, unadj_H0_bad = 999,
                       unadj_H1 = 999, unadj_H1_good = 999, unadj_H1_bad = 999)
  
  
  nsim = v$nsim  
  
  covH0 <- compcov(parallel_filenames_H0)
  covH1 <- compcov(parallel_filenames_H1)
  
  ### For ltmle ###
  
  # extract the mean and covariance structure of the simulated estimators
  mvmean0 <- covH0$mean_ltmle_std
  mvmean1 <- covH1$mean_ltmle_std
  mvcov0 <- covH0$cov_ltmle_std
  mvcov1 <- covH1$cov_ltmle_std
  
  # compute the error that will be spent at each interim analysis
  errs <- errs(covH0$cov_ltmle, v)
  
  # compute the error spending boundry, with the given multivariate normal distribution structure
  bdry_ltmle <- err_sp_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  
  ## H0 ##
  # the result (when to stop, and trial result) of each trial, under H0
  results_ltmle_H0 <- trialresults(covH0$mat_ltmle_std, bdry_ltmle)
  
  # compute the reversal probabilities, of H0 
  reversals <- subset(results_ltmle_H0, (whichstop != 5) & (interimresult != decisionresult))
  output$ltmle_H0 <- nrow(reversals) / nsim
  output$ltmle_H0_good <- nrow(subset(reversals, decisionresult == 0)) / nsim
  output$ltmle_H0_bad <- nrow(subset(reversals, decisionresult == 1)) / nsim

  ## H1 ##
  # the result (when to stop, and trial result) of each trial, under H1
  results_ltmle_H1 <- trialresults(covH1$mat_ltmle_std, bdry_ltmle)
  
  # compute the reversal probabilities, of H1   
  reversals <- subset(results_ltmle_H1, (whichstop != 5) & (interimresult != decisionresult))
  output$ltmle_H1 <- nrow(reversals) / nsim
  output$ltmle_H1_good <- nrow(subset(reversals, decisionresult == 1)) / nsim
  output$ltmle_H1_bad <- nrow(subset(reversals, decisionresult == 0)) / nsim
  
  
  ### For unadj ###
  
  # extract the mean and covariance structure of the simulated estimators
  mvmean0 <- covH0$mean_unadj_std
  mvmean1 <- covH1$mean_unadj_std
  mvcov0 <- covH0$cov_unadj_std
  mvcov1 <- covH1$cov_unadj_std
  
  # compute the error that will be spent at each interim analysis
  errs <- errs(covH0$cov_unadj, v)
  
  # compute the error spending boundry, with the given multivariate normal distribution structure
  bdry_unadj <- err_sp_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  
  ## H0 ##
  # the result (when to stop, and trial result) of each trial, under H0
  results_unadj_H0 <- trialresults(covH0$mat_unadj_std, bdry_unadj)
  
  # compute the reversal probabilities, of H0 
  reversals <- subset(results_unadj_H0, (whichstop != 5) & (interimresult != decisionresult))
  output$unadj_H0 <- nrow(reversals) / nsim
  output$unadj_H0_good <- nrow(subset(reversals, decisionresult == 0)) / nsim
  output$unadj_H0_bad <- nrow(subset(reversals, decisionresult == 1)) / nsim
  
  ## H1 ##
  # the result (when to stop, and trial result) of each trial, under H1
  results_unadj_H1 <- trialresults(covH1$mat_unadj_std, bdry_unadj)
  
  # compute the reversal probabilities, of H1   
  reversals <- subset(results_unadj_H1, (whichstop != 5) & (interimresult != decisionresult))
  output$unadj_H1 <- nrow(reversals) / nsim
  output$unadj_H1_good <- nrow(subset(reversals, decisionresult == 1)) / nsim
  output$unadj_H1_bad <- nrow(subset(reversals, decisionresult == 0)) / nsim
  
  return(output)
}








