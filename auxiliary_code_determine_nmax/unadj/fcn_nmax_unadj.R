
library(ltmle)

source("fcn_err_spending_boundry.R")
source("fcn_compcov.R")
source("fcn_trialresults.R")

simtrials.unadjonly <- function(dt, v){
  
  # parameters set-up
  nsim <- v$nsim
  nmax <- v$nmax
  K <- v$K
  enrollrate <- v$enrollrate
  randomA <- v$randomA
  A_to_L1 <- v$A_to_L1
  L1_to_Y <- v$L1_to_Y
  interim_times  <- (1:K)/K * nmax/enrollrate + A_to_L1 + L1_to_Y # so that we have (k/K)*nmax patients with final outcome at k-th inteirm analysis
  decision_times <- pmin(interim_times + A_to_L1 + L1_to_Y, interim_times[K])
  # at K-th interim analysis, we have everyone's final outcome observed.
  
  # output initializaiton
  output <- list(unadj_est = matrix(NA, nrow = nsim, ncol = 2*K))
  
  # Step 0: generate twins
  dtaug <- DGM.step0.twins(dt)  
  
  
  for( itrial in 1:nsim ){
    
    if( itrial %% 10 == 0){
      cat(paste0("\nTrial number ", itrial))
    }
    
    # Sample from the twin data with replacement, to have desired sample size
    ids <- sample(nrow(dtaug), nmax, replace = TRUE)
    dtc <- dtaug[ids, ]
    
    # Step 1: Calibrate Y with pY, to mimic Treatment Effect
    dtc <- DGM.step1.pY.mimicTE(dtc, pY = v$pY, direction = v$dir)
    
    # Step 2: Calibrate Y with e0 and e1, to mimic Relative Efficiency
    dtc <- DGM.step2.e.mimicRE(dtc, e = v$e, p0 = v$p0, p1 = v$p1)
    
    if (randomA == TRUE) { # Reset A as Bernoulli(0.5), i.e. assume no effect
      dtc$A <- rbinom(nrow(dtc), 1, 0.5)      
    }
    
    # Assign enrollment time for each patient
    
    dtc$enrolltime <- sample(1:nmax, nmax) / enrollrate
    dtc$L1time <- dtc$enrolltime + A_to_L1
    dtc$Ytime  <- dtc$L1time + L1_to_Y
    
    for( istage in 1:K ){ # Do each interim analysis and decision analysis, compute the estimators
      
      ## Interim analysis
      
      intrtime <- interim_times[istage]
      
      C0 <- as.numeric(dtc$enrolltime <= intrtime)
      C1 <- as.numeric(dtc$L1time     <= intrtime)
      C2 <- as.numeric(dtc$Ytime      <= intrtime)
      
      dttmle <- data.frame(C0 = C0, S = dtc$S, W1 = dtc$W1, W2 = dtc$W2, A = dtc$A,
                           C1 = C1, L1 = dtc$L1,
                           C2 = C2, Y = dtc$Y)
      
      output$unadj_est[itrial, 2*istage - 1] <- 
        mean(subset(dttmle, (C2 == 1) & (A == 1))$Y) - mean(subset(dttmle, (C2 == 1) & (A == 0))$Y)
      
      ## Decision analysis
      
      dcsntime <- decision_times[istage]
      
      C0 <- as.numeric(dtc$enrolltime <= dcsntime)
      C1 <- as.numeric(dtc$L1time     <= dcsntime)
      C2 <- as.numeric(dtc$Ytime      <= dcsntime)
      
      dttmle <- data.frame(C0 = C0, S = dtc$S, W1 = dtc$W1, W2 = dtc$W2, A = dtc$A,
                           C1 = C1, L1 = dtc$L1,
                           C2 = C2, Y = dtc$Y)
      
      output$unadj_est[itrial, 2*istage] <- 
        mean(subset(dttmle, (C2 == 1) & (A == 1))$Y) - mean(subset(dttmle, (C2 == 1) & (A == 0))$Y)
      
    } # end for loop: istage
  } # end for loop: itrial
  
  return(output)
}



compcov.unadjonly <- function(generic_filename, idx, list_name = "test", directory = ""){
  # generic_filename: the part of the filename that is shared across all the result files.
  # idx: the indices of the parallel jobs (for example, 1:50)
  # list_name: the name of the list that is saved in the result .rda file. Default: "test".
  # directory: the directory where the .rda results are stored.
  
  n <- length(idx)
  idx <- as.character(idx)
  
  # initialize overall matrix to collect all parallel simulated esitmators
  mat_unadj <- matrix(NA, nrow = 1, ncol = 10)
  
  # rbind each resulting matrix to the overall matrix
  for (ifile in 1:n){
    load(paste0(directory, generic_filename, "_", idx[ifile], ".rda"))
    eval(parse(text = paste0("mat_unadj <- rbind(mat_unadj, ", list_name, "$unadj_est)")))
  }
  
  # delete the first line (initializer)
  mat_unadj <- mat_unadj[-1, ]
  
  # Compute mean and covariance of ltmle and unadj
  mean_unadj <- apply(mat_unadj, 2, mean)
  cov_unadj <- cov(mat_unadj)
  
  # Standardize ltmle and unadj (divide by standard error, make it Z-statisitc)
  mat_unadj_std <- mat_unadj %*% diag(diag(cov_unadj) ^ (-1/2))
  
  # Compute mean and covariance of the standardized ones
  mean_unadj_std <- apply(mat_unadj_std, 2, mean)
  cov_unadj_std <- cov(mat_unadj_std)
  
  # gather the cov/corr results
  result <- list(mat_unadj = mat_unadj,
                 mean_unadj = mean_unadj,
                 cov_unadj = cov_unadj,
                 mat_unadj_std = mat_unadj_std,
                 mean_unadj_std = mean_unadj_std,
                 cov_unadj_std = cov_unadj_std)
  return(result)
}


comppower.unadjonly <- function(v, caseH0, caseH1, idx = 1:500, nonbinding = FALSE){
  # v needs to be a list of: nsim, K, alpha, beta, rho, f_err, g_err
  #   note: nsim should be the total number, i.e. add all parallel jobs up
  # caseH0, caseH1 are the file names of the stored result, e.g. "est_H0" and "est_H1"
  # idx is the index of parallel result files, currently are 1:100. (used in compcov() function)
  
  
  ## output: vector of length 4,
  ## "ltmle_H0" means the total type 1 error of ltmle under H0
  ## "ltmle_H1" means the total power of ltmle under H1
  ## similar to unadj.
  output <- rep(NA, 2)
  names(output) <- c("unadj_H0", "unadj_H1")
  
  
  nsim = v$nsim
  
  covH0 <- compcov.unadjonly(caseH0, idx)
  covH1 <- compcov.unadjonly(caseH1, idx)
  
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
  
    print("unadj: errors spent at each stage:")
    print(errs)
    print("unadj: error spending boundary:")
    print(bdry_unadj)
  
  ## unadj_H0
  results <- trialresults(covH0$mat_unadj_std, bdry_unadj)
  output[1] <- mean(results$decisionresult == 1)
  
  ## unadj_H1
  results <- trialresults(covH1$mat_unadj_std, bdry_unadj)
  output[2] <- mean(results$decisionresult == 1)
  
  return(output)
}
