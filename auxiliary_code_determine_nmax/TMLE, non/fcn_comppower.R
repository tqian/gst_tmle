library(ltmle)
library(mvtnorm)

source("fcn_err_spending_boundry.R")
source("fcn_compcov.R")
source("fcn_trialresults.R")

comppower <- function(v, caseH0, caseH1, idx = 0:99, nonbinding = FALSE){
  # v needs to be a list of: nsim, K, alpha, beta, rho, f_err, g_err
  #   note: nsim should be the total number, i.e. add all parallel jobs up
  # caseH0, caseH1 are the file names of the stored result, e.g. "est_H0" and "est_H1"
  # idx is the index of parallel result files, currently are 1:100. (used in compcov() function)
  
  
  ## output: vector of length 4,
  ## "ltmle_H0" means the total type 1 error of ltmle under H0
  ## "ltmle_H1" means the total power of ltmle under H1
  ## similar to unadj.
  output <- rep(NA, 4)
  names(output) <- c("ltmle_H0", "ltmle_H1", "unadj_H0", "unadj_H1")
  
  
  nsim = v$nsim
  
  covH0 <- compcov(caseH0, idx)
  covH1 <- compcov(caseH1, idx)
  
  ### For ltmle ###
  
  # extract the mean and covariance structure of the simulated estimators
  mvmean0 <- covH0$mean_ltmle_std
  mvmean1 <- covH1$mean_ltmle_std
  mvcov0 <- covH0$cov_ltmle_std
  mvcov1 <- covH1$cov_ltmle_std
  
  # compute the error that will be spent at each interim analysis
  errs <- compute.errs(covH0$cov_ltmle, v)
  
  # compute the error spending boundry, with the given multivariate normal distribution structure
  if (isTRUE(nonbinding)){
    bdry_ltmle <- err_spending_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  } else {
    bdry_ltmle <- err_spending_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  }
  print("ltmle: errors spent at each stage:")
  print(errs)
  print("ltmle: error spending boundary:")
  print(bdry_ltmle)
  
  ## ltmle_H0
  results <- trialresults(covH0$mat_ltmle_std, bdry_ltmle)
  output[1] <- mean(results$decisionresult == 1)
  
  ## ltmle_H1
  results <- trialresults(covH1$mat_ltmle_std, bdry_ltmle)
  output[2] <- mean(results$decisionresult == 1)
  
  ### For unadj ###
  
  # extract the mean and covariance structure of the simulated estimators
  mvmean0 <- covH0$mean_unadj_std
  mvmean1 <- covH1$mean_unadj_std
  mvcov0 <- covH0$cov_unadj_std
  mvcov1 <- covH1$cov_unadj_std
  
  # compute the error that will be spent at each interim analysis
  errs <- compute.errs(covH0$cov_unadj, v)
  
  # compute the error spending boundry, with the given multivariate normal distribution structure
  if (isTRUE(nonbinding)){
    bdry_unadj <- err_spending_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  } else {
    bdry_unadj <- err_spending_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  }

#   print("unadj: errors spent at each stage:")
#   print(errs)
#   print("unadj: error spending boundary:")
#   print(bdry_unadj)
  
  ## unadj_H0
  results <- trialresults(covH0$mat_unadj_std, bdry_unadj)
  output[3] <- mean(results$decisionresult == 1)
  
  ## unadj_H1
  results <- trialresults(covH1$mat_unadj_std, bdry_unadj)
  output[4] <- mean(results$decisionresult == 1)
    
  return(output)
}






