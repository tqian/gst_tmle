source("source/fcn_compcov.R")

get_jobname <- function(directory = "Result", include_dir = TRUE){
  
  if (substr(directory, nchar(directory), nchar(directory)) != "/"){
    directory <- paste0(directory, "/")
  }
  
  files <- list.files(directory)
  unique_jobs <- unique(sub("_[0-9]{1,2}.rda", "", files))
  
  jobs_H0 <- unique_jobs[which(grepl("H0", unique_jobs))]
  jobs_H1 <- unique_jobs[which(grepl("H1", unique_jobs))]
  
  if (length(jobs_H0) + length(jobs_H1) != length(unique_jobs)) {
    stop("length(jobs_H0) + length(jobs_H1) != length(unique_jobs)")
  }
  if (length(jobs_H0) != length(jobs_H1)) {
    stop("length(jobs_H0) != length(jobs_H1)")
  }
  
  if (include_dir) {
    jobs_H0 <- paste0(directory, jobs_H0)
    jobs_H1<- paste0(directory, jobs_H1)
  }  
  
  return(list(H0 = sort(jobs_H0), H1 = sort(jobs_H1)))
}


percentage_stopping <- function(v, parallel_filenames_H0, parallel_filenames_H1, nonbinding = FALSE, rho = 2,
                                estimator = c("ltmle", "unadj"),
                                digits = 1, print_more = FALSE, scn){
  ##########
  # Return a list of: a matrix of percentage stopping & rejection under H0 and H1;
  #                   a table of nmax, ESS_H0, ESS_H1 (expected sample size.)
  #
  ##########

  nmax <- v$nmax
  erate <- v$enrollrate
  nsim <- v$nsim
  
  t_AL1 <- v$A_to_L1
  t_L1Y <- v$L1_to_Y
  num_full_obs_at_interim <- floor(seq(0.2, 1, 0.2)*nmax)
  t_at_interim <- num_full_obs_at_interim/erate + t_AL1 + t_L1Y
  num_W_obs_at_interim <- pmin(t_at_interim * erate, nmax)
  
  covH0 <- compcov(parallel_filenames_H0)
  covH1 <- compcov(parallel_filenames_H1)
  
  # compute errors, boundaries, and trial results
  if (identical(estimator, c("ltmle", "unadj")) | identical(estimator, "ltmle")){    
    mvmean0 <- covH0$mean_ltmle_std
    mvmean1 <- covH1$mean_ltmle_std
    mvcov0 <- covH0$cov_ltmle_std
    mvcov1 <- covH1$cov_ltmle_std    
    if (nonbinding){
      errs <- errs_nonbinding(covH0$cov_ltmle, v, rho = rho)
      bdry_ltmle <- err_sp_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    } else {
      errs <- errs(covH0$cov_ltmle, v, rho = rho)
      bdry_ltmle <- err_sp_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    }
    resultsH0 <- trialresults(covH0$mat_ltmle_std, bdry_ltmle)
    resultsH1 <- trialresults(covH1$mat_ltmle_std, bdry_ltmle)
  } else if (estimator == "unadj"){
    mvmean0 <- covH0$mean_unadj_std
    mvmean1 <- covH1$mean_unadj_std
    mvcov0 <- covH0$cov_unadj_std
    mvcov1 <- covH1$cov_unadj_std    
    if (nonbinding){
      errs <- errs_nonbinding(covH0$cov_unadj, v, rho = rho)
      bdry_unadj <- err_sp_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    } else {
      errs <- errs(covH0$cov_unadj, v, rho = rho)
      bdry_unadj <- err_sp_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    }
    resultsH0 <- trialresults(covH0$mat_unadj_std, bdry_unadj)
    resultsH1 <- trialresults(covH1$mat_unadj_std, bdry_unadj)
  }
  
  # gather the results into a table
  stop_table_fut <- matrix(NA, nrow = 4, ncol = 6)
  stop_table_eff <- matrix(NA, nrow = 4, ncol = 6)
  rownames(stop_table_fut) <- c("H0_interim_futstop(%)", "H0_decision_accept(%)",
                            "H1_interim_futstop(%)", "H1_decision_accept(%)")
  rownames(stop_table_eff) <- c("H0_interim_effstop(%)", "H0_decision_reject(%)",
                                     "H1_interim_effstop(%)", "H1_decision_reject(%)")
  colnames(stop_table_fut) <- c("1st", "2nd", "3rd", "4th", "final", "sum")
  colnames(stop_table_eff) <- c("1st", "2nd", "3rd", "4th", "final", "sum")
  ESS <- matrix(NA, nrow = 1, ncol = 3)
  colnames(ESS) <- c("nmax", "ESS_H0", "ESS_H1")
  rownames(ESS) <- "sample size"
  ESS[1] <- nmax  
  
  # count interim stops and decision rejections for H0 and for H1
  for (H in 0:1){
    if (H == 0){stop_data <- resultsH0}
    if (H == 1){stop_data <- resultsH1}
    # temp matrix for counting stops
    stop_count_fut <- matrix(NA, nrow = 2, ncol = 6)
    stop_count_eff <- matrix(NA, nrow = 2, ncol = 6)
    for (istop in 1:5){
      stop_at_i <- subset(stop_data, whichstop == istop)
      stop_count_fut[1, istop] <- sum(stop_at_i$interimresult == 0)
      stop_count_fut[2, istop] <- sum(stop_at_i$decisionresult == 0)
      stop_count_eff[1, istop] <- sum(stop_at_i$interimresult == 1)
      stop_count_eff[2, istop] <- sum(stop_at_i$decisionresult == 1)
    } 
    stop_count_fut[1, 6] <- sum(stop_count_fut[1, 1:5])
    stop_count_fut[2, 6] <- sum(stop_count_fut[2, 1:5]) 
    stop_count_eff[1, 6] <- sum(stop_count_eff[1, 1:5])
    stop_count_eff[2, 6] <- sum(stop_count_eff[2, 1:5])   
    # compute percentage stopping
    stop_table_fut[(H*2+1):(H*2+2), ] <- stop_count_fut[1:2, ] / nsim * 100
    stop_table_eff[(H*2+1):(H*2+2), ] <- stop_count_eff[1:2, ] / nsim * 100

    tmp <- stop_table_fut[(H*2+1), 1:4] + stop_table_eff[(H*2+1), 1:4]
    tmp <- c(tmp, 100-sum(tmp))
    # This may not be correct! Needs to do with enrollment rate and number of pipeline patients
    ESS[H+2] <- tmp %*% num_W_obs_at_interim / 100
  }
  
  # print table
  cat("####################################\n")
  cat(sprintf("scnario: %s %s, nonbinding = %s, rho = %d.\n",
              as.vector(estimator[1]), scn, nonbinding, rho))
  print(round(ESS,0))
  if (print_more) {
    cat("Stopping for futility:\n")
    print(round(stop_table_fut, digits))
    cat("Stopping for efficacy:\n")
    print(round(stop_table_eff, digits))
    cat(paste0("\nalphas: ", paste(round(errs[1,],4), collapse = " "), "\n"))
    cat(paste0("betas:  ", paste(round(errs[2,],3), collapse = "  "), "\n"))
  }
  return(list(stop_table_fut = stop_table_fut,
              stop_table_eff = stop_table_eff,
              ESS = ESS,
              resultsH0 = resultsH0, resultsH1 = resultsH1))
}