rm(list = ls())

library(ltmle)

setwd("MISTIE/code_determine_nmax")


source("fcn_comppower.R")
source("fcn_nmax_unadj.R")

v <- list(nsim = 50000,
          nmax = NA,
          K = 5,
          enrollrate = 140,
          pY = 0.02955846, # calibration probability
          dir = "up", # for calibration
          e = 0.1475, # calibration probability
          p0 = 0.2222222, # for calibration, p0 = P(Y=1|A=0)
          p1 = 0.34375, # for calibration, p1 = P(Y=1|A=1)
          randomA = FALSE,
          exL = FALSE,
          exW = FALSE,
          A_to_L1 = 30 / 365,
          L1_to_Y = (180 - 30) / 365,
          alpha = 0.025,
          beta = 0.2,
          rho = 2,
          f_err = 
            function(x, .rho = 2, .alpha = 0.025){
              if (x < 0){
                f <- 0
              } else if (x < 1){
                f <- .alpha * x^.rho
              } else{
                f <- .alpha
              }
              return(f)
            },
          g_err =
            function(x, .rho = 2, .beta = 0.2){
              if (x < 0){
                g <- 0
              } else if (x < 1){
                g <- .beta * x^.rho
              } else{
                g <- .beta
              }
              return(g)
            }
)

nameH0s <- c("Result/nmax_H0_enrollrate140_4",
             "Result/nmax_H0exL_enrollrate140_3",
             "Result/nmax_H0exW_enrollrate140_5",
             "Result/nmax_H0exWL_enrollrate140_3",
             "Result/nmax_H0unadj_enrollrate140_1")
nameH1s <- c("Result/nmax_H1_enrollrate140_4",
             "Result/nmax_H1exL_enrollrate140_3",
             "Result/nmax_H1exW_enrollrate140_5",
             "Result/nmax_H1exWL_enrollrate140_3",
             "Result/nmax_H1unadj_enrollrate140_1")
nmaxs <- c(336, 340, 442, 436, 460)

Npara <- 500
nonbinding <- TRUE

for (iscn in 1:5){
  if (iscn %in% 1:4){ # for ltmle
    v$nmax <- nmaxs[iscn]
    resultnamesH0 <- paste0(nameH0s[iscn], "_", 1:Npara, ".rda") # result filenames from parallel jobs
    resultnamesH1 <- paste0(nameH1s[iscn], "_", 1:Npara, ".rda") # result filenames from parallel jobs
    
    nsim = v$nsim
    covH0 <- compcov(nameH0s[iscn], 1:Npara)
    covH1 <- compcov(nameH1s[iscn], 1:Npara)
    mvmean0 <- covH0$mean_ltmle_std
    mvmean1 <- covH1$mean_ltmle_std
    mvcov0 <- covH0$cov_ltmle_std
    mvcov1 <- covH1$cov_ltmle_std
    
    errs <- compute.errs(covH0$cov_ltmle, v)
    if (isTRUE(nonbinding)){
      bdry_ltmle <- err_spending_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    } else {
      bdry_ltmle <- err_spending_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    }
    results <- trialresults(covH1$mat_ltmle_std, bdry_ltmle) # to see early stops under H1
  } else { # for unadj
    v$nmax <- nmaxs[iscn]
    resultnamesH0 <- paste0(nameH0s[iscn], "_", 1:Npara, ".rda") # result filenames from parallel jobs
    resultnamesH1 <- paste0(nameH1s[iscn], "_", 1:Npara, ".rda") # result filenames from parallel jobs
    
    nsim = v$nsim
    covH0 <- compcov.unadjonly(nameH0s[iscn], 1:Npara)
    covH1 <- compcov.unadjonly(nameH1s[iscn], 1:Npara)
    mvmean0 <- covH0$mean_unadj_std
    mvmean1 <- covH1$mean_unadj_std
    mvcov0 <- covH0$cov_unadj_std
    mvcov1 <- covH1$cov_unadj_std
    
    errs <- compute.errs(covH0$cov_unadj, v)
    if (isTRUE(nonbinding)){
      bdry_unadj <- err_spending_bdry_nonbinding(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    } else {
      bdry_unadj <- err_spending_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
    }
    results <- trialresults(covH1$mat_unadj_std, bdry_unadj) # to see early stops under H1
  }

  print("######################")
  print(nameH1s[iscn])
  print(mean(results$whichstop))
}
# [1] "######################"
# [1] "Result/nmax_H1_enrollrate140_4"
# [1] 3.67478
# [1] "######################"
# [1] "Result/nmax_H1exL_enrollrate140_3"
# [1] 3.73982
# [1] "######################"
# [1] "Result/nmax_H1exW_enrollrate140_5"
# [1] 3.58824
# [1] "######################"
# [1] "Result/nmax_H1exWL_enrollrate140_3"
# [1] 3.7439
# [1] "######################"
# [1] "Result/nmax_H1unadj_enrollrate140_1"
# [1] 3.71406