library(ltmle)
library(mvtnorm)

source("source/fcn_err_sp_bdry.R")
source("source/fcn_compcov.R")
source("source/fcn_trialresults.R")

comppower <- function(v, parallel_filenames_H0, parallel_filenames_H1, nonbinding = FALSE, rho = 2, only_unadj = FALSE){
  # v needs to be a list of: nsim, K, alpha, beta, rho, f_err, g_err
  #   note: nsim should be the total number, i.e. add all parallel jobs up
  # parallel_filenames_H0, parallel_filenames_H1 are the filenames of the stored result
  
  
  ## output: vector of length 4,
  ## "ltmle_H0" means the total type 1 error of ltmle under H0
  ## "ltmle_H1" means the total power of ltmle under H1
  ## similar to unadj.
  output <- rep(NA, 4)
  names(output) <- c("ltmle_H0", "ltmle_H1", "unadj_H0", "unadj_H1")
  
  covH0 <- compcov(parallel_filenames_H0)
  covH1 <- compcov(parallel_filenames_H1)
  
  nsim <- v$nsim
  
  if (!only_unadj){
    ### For ltmle ###
    
    # extract the mean and covariance structure of the simulated estimators
    mvmean0 <- covH0$mean_ltmle_std
    mvmean1 <- covH1$mean_ltmle_std
    mvcov0 <- covH0$cov_ltmle_std
    mvcov1 <- covH1$cov_ltmle_std
    
    # compute the error spending boundry, with the given multivariate normal distribution structure
    if (nonbinding){
      errs <- errs_nonbinding(covH0$cov_ltmle, v, rho = rho)
      bdry_ltmle <- err_sp_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    } else {
      errs <- errs(covH0$cov_ltmle, v, rho = rho)
      bdry_ltmle <- err_sp_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    }
    
    print_errs_bdry(errs, bdry_ltmle, "ltmle")
    
    ## ltmle_H0
    results <- trialresults(covH0$mat_ltmle_std, bdry_ltmle)
    output[1] <- mean(results$decisionresult == 1)
    
    ## ltmle_H1
    results <- trialresults(covH1$mat_ltmle_std, bdry_ltmle)
    output[2] <- mean(results$decisionresult == 1)
  }
  
  ### For unadj ###
  
  # extract the mean and covariance structure of the simulated estimators
  mvmean0 <- covH0$mean_unadj_std
  mvmean1 <- covH1$mean_unadj_std
  mvcov0 <- covH0$cov_unadj_std
  mvcov1 <- covH1$cov_unadj_std
  
  # compute the error spending boundry, with the given multivariate normal distribution structure
  if (nonbinding){
    errs <- errs_nonbinding(covH0$cov_unadj, v, rho = rho)
    bdry_unadj <- err_sp_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  } else {
    errs <- errs(covH0$cov_unadj, v, rho = rho)
    bdry_unadj <- err_sp_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
  }

  print_errs_bdry(errs, bdry_unadj, "unadj")
  
  ## unadj_H0
  results <- trialresults(covH0$mat_unadj_std, bdry_unadj)
  output[3] <- mean(results$decisionresult == 1)
  
  ## unadj_H1
  results <- trialresults(covH1$mat_unadj_std, bdry_unadj)
  output[4] <- mean(results$decisionresult == 1)
    
  return(output)
}


print_errs_bdry <- function(errs, bdry, name = c("ltmle", "unadj")){
  ##########
  # Print out computed errs and bdry.
  #
  # Keyword Arguments:
  #     errs -- errors spent at each stage
  #     bdry -- boundaries (u, l, c) at each stage
  #     name -- "ltmle" or "unadj"
  ##########
  alphas <- errs[1, ]
  betas <- errs[2, ]
  u <- bdry$u
  c <- bdry$c
  l <- bdry$l
  cat(sprintf("\n%s: errors and boundaries at each stage:\n", name))
#   cat(sprintf("stage   %6.0f   %6.0f   %6.0f   %6.0f   %6.0f\n",
#               1, 2, 3, 4, 5))
  cat(sprintf("alpha   %1.4f   %1.4f   %1.4f   %1.4f   %1.4f\n",
              alphas[1], alphas[2], alphas[3], alphas[4], alphas[5]))
  cat(sprintf("beta    %1.4f   %1.4f   %1.4f   %1.4f   %1.4f\n",
              betas[1],betas[2],betas[3],betas[4],betas[5]))
  cat(sprintf("u       %3.2f   %3.2f   %3.2f   %3.2f   %3.2f\n",
              u[1],u[2],u[3],u[4],u[5]))
  cat(sprintf("c       %3.2f   %3.2f   %3.2f   %3.2f   %3.2f\n",
              c[1],c[2],c[3],c[4],c[5]))
  cat(sprintf("l       %3.2f   %3.2f   %3.2f   %3.2f   %3.2f\n",
              l[1],l[2],l[3],l[4],l[5]))
}