rm(list = ls())

library(ltmle)
library(mvtnorm)


source("source/fcn_comppower.R")
source("source/fcn_reverse_prob.R")
source("source/fcn_evaltrial.R")




# collect results from parallel simulation --------------------------------

nsim <- 50000
K <- 5 # number of stages
est_mat_interim <- matrix(0, nrow = 0, ncol = K - 1)
est_mat_decision <- matrix(0, nrow = 0, ncol = K)

est_tmp_list <- list(interim = list(
  tx = list(ltmle = est_mat_interim, unadj = est_mat_interim),
  a0 = list(ltmle = est_mat_interim, unadj = est_mat_interim),
  a1 = list(ltmle = est_mat_interim, unadj = est_mat_interim)),
  decision = list(
    tx = list(ltmle = est_mat_decision, unadj = est_mat_decision),
    a0 = list(ltmle = est_mat_decision, unadj = est_mat_decision),
    a1 = list(ltmle = est_mat_decision, unadj = est_mat_decision)
  ))
# estimator for: tx: treatment effect; a0: A=0; a1: A=1

collected_est <- list(prognWL = list(H0 = est_tmp_list, H1 = est_tmp_list),
                      prognW = list(H0 = est_tmp_list, H1 = est_tmp_list),
                      prognL = list(H0 = est_tmp_list, H1 = est_tmp_list),
                      prognnon = list(H0 = est_tmp_list, H1 = est_tmp_list))

npara <- 10

interim_index <- c(1,3,5,7)
decision_index <- c(2,4,6,8,10)

for (case_id in 1:4) {
  
  if (case_id == 1) {
    exW <- FALSE
    exL <- FALSE
    nmax <- 300
  } else if (case_id == 2) {
    exW <- FALSE
    exL <- TRUE
    nmax <- 300
  } else if (case_id == 3) {
    exW <- TRUE
    exL <- FALSE
    nmax <- 480
  } else if (case_id == 4) {
    exW <- TRUE
    exL <- TRUE
    nmax <- 480
  }
  
  if (exL & exW) {
    scn <- "prognnon"
  } else if (exL) {
    scn <- "prognW"
  } else if (exW) {
    scn <- "prognL"
  } else {
    scn <- "prognWL"
  }
  
  taskIDs <- 1:npara
  
  filenames_H0 <- paste0("result/est_", scn, "_nmax", nmax, "_H0_", taskIDs, ".RDS")
  filenames_H1 <- paste0("result/est_", scn, "_nmax", nmax, "_H1_", taskIDs, ".RDS")
  
  for (filename in filenames_H0) {
    result <- readRDS(filename)
    collected_est[[scn]]$H0$interim$tx$ltmle <- rbind(collected_est[[scn]]$H0$interim$tx$ltmle,
                                                      result$ltmle_est[, interim_index])
    collected_est[[scn]]$H0$interim$a0$ltmle <- rbind(collected_est[[scn]]$H0$interim$a0$ltmle,
                                                      result$ltmle_est_0[, interim_index])
    collected_est[[scn]]$H0$interim$a1$ltmle <- rbind(collected_est[[scn]]$H0$interim$a1$ltmle,
                                                      result$ltmle_est_1[, interim_index])
    collected_est[[scn]]$H0$interim$tx$unadj <- rbind(collected_est[[scn]]$H0$interim$tx$unadj,
                                                      result$unadj_est[, interim_index])
    collected_est[[scn]]$H0$interim$a0$unadj <- rbind(collected_est[[scn]]$H0$interim$a0$unadj,
                                                      result$unadj_est_0[, interim_index])
    collected_est[[scn]]$H0$interim$a1$unadj <- rbind(collected_est[[scn]]$H0$interim$a1$unadj,
                                                      result$unadj_est_1[, interim_index])
    collected_est[[scn]]$H0$decision$tx$ltmle <- rbind(collected_est[[scn]]$H0$decision$tx$ltmle,
                                                       result$ltmle_est[, decision_index])
    collected_est[[scn]]$H0$decision$a0$ltmle <- rbind(collected_est[[scn]]$H0$decision$a0$ltmle,
                                                       result$ltmle_est_0[, decision_index])
    collected_est[[scn]]$H0$decision$a1$ltmle <- rbind(collected_est[[scn]]$H0$decision$a1$ltmle,
                                                       result$ltmle_est_1[, decision_index])
    collected_est[[scn]]$H0$decision$tx$unadj <- rbind(collected_est[[scn]]$H0$decision$tx$unadj,
                                                       result$unadj_est[, decision_index])
    collected_est[[scn]]$H0$decision$a0$unadj <- rbind(collected_est[[scn]]$H0$decision$a0$unadj,
                                                       result$unadj_est_0[, decision_index])
    collected_est[[scn]]$H0$decision$a1$unadj <- rbind(collected_est[[scn]]$H0$decision$a1$unadj,
                                                       result$unadj_est_1[, decision_index])
  }
  for (filename in filenames_H1) {
    result <- readRDS(filename)
    collected_est[[scn]]$H1$interim$tx$ltmle <- rbind(collected_est[[scn]]$H1$interim$tx$ltmle,
                                                      result$ltmle_est[, interim_index])
    collected_est[[scn]]$H1$interim$a0$ltmle <- rbind(collected_est[[scn]]$H1$interim$a0$ltmle,
                                                      result$ltmle_est_0[, interim_index])
    collected_est[[scn]]$H1$interim$a1$ltmle <- rbind(collected_est[[scn]]$H1$interim$a1$ltmle,
                                                      result$ltmle_est_1[, interim_index])
    collected_est[[scn]]$H1$interim$tx$unadj <- rbind(collected_est[[scn]]$H1$interim$tx$unadj,
                                                      result$unadj_est[, interim_index])
    collected_est[[scn]]$H1$interim$a0$unadj <- rbind(collected_est[[scn]]$H1$interim$a0$unadj,
                                                      result$unadj_est_0[, interim_index])
    collected_est[[scn]]$H1$interim$a1$unadj <- rbind(collected_est[[scn]]$H1$interim$a1$unadj,
                                                      result$unadj_est_1[, interim_index])
    collected_est[[scn]]$H1$decision$tx$ltmle <- rbind(collected_est[[scn]]$H1$decision$tx$ltmle,
                                                       result$ltmle_est[, decision_index])
    collected_est[[scn]]$H1$decision$a0$ltmle <- rbind(collected_est[[scn]]$H1$decision$a0$ltmle,
                                                       result$ltmle_est_0[, decision_index])
    collected_est[[scn]]$H1$decision$a1$ltmle <- rbind(collected_est[[scn]]$H1$decision$a1$ltmle,
                                                       result$ltmle_est_1[, decision_index])
    collected_est[[scn]]$H1$decision$tx$unadj <- rbind(collected_est[[scn]]$H1$decision$tx$unadj,
                                                       result$unadj_est[, decision_index])
    collected_est[[scn]]$H1$decision$a0$unadj <- rbind(collected_est[[scn]]$H1$decision$a0$unadj,
                                                       result$unadj_est_0[, decision_index])
    collected_est[[scn]]$H1$decision$a1$unadj <- rbind(collected_est[[scn]]$H1$decision$a1$unadj,
                                                       result$unadj_est_1[, decision_index])
  }
}

saveRDS(collected_est, file = "collected_est_nsim50000.RDS")





# compute relative efficiency ---------------------------------------------

for (scn in names(collected_est)) {
  for (hypothesis in c("H0", "H1")) {
    for (analysis_type in c("interim", "decision")) {
      for (estimand in c("tx", "a0", "a1")) {
        var_ltmle <- apply(collected_est[[scn]][[hypothesis]][[analysis_type]][[estimand]]$ltmle, 2, var)
        var_unadj <- apply(collected_est[[scn]][[hypothesis]][[analysis_type]][[estimand]]$unadj, 2, var)
        collected_est[[scn]][[hypothesis]][[analysis_type]][[estimand]]$RE <- var_unadj / var_ltmle
      }
    }
  }
}


##### For Table 3 in the paper: #####
# "RE from simulation" columns in the table

options(digits = 3)

cbind(collected_est$prognWL$H0$interim$tx$RE,
      collected_est$prognW$H0$interim$tx$RE,
      collected_est$prognL$H0$interim$tx$RE,
      collected_est$prognnon$H0$interim$tx$RE)

cbind(collected_est$prognWL$H0$decision$tx$RE,
      collected_est$prognW$H0$decision$tx$RE,
      collected_est$prognL$H0$decision$tx$RE,
      collected_est$prognnon$H0$decision$tx$RE)

cbind(collected_est$prognWL$H1$interim$tx$RE,
      collected_est$prognW$H1$interim$tx$RE,
      collected_est$prognL$H1$interim$tx$RE,
      collected_est$prognnon$H1$interim$tx$RE)

cbind(collected_est$prognWL$H1$decision$tx$RE,
      collected_est$prognW$H1$decision$tx$RE,
      collected_est$prognL$H1$decision$tx$RE,
      collected_est$prognnon$H1$decision$tx$RE)



##### For Table E.1 in the paper: #####
# "RE from simulation" columns in the table

options(digits = 3)

cbind(collected_est$prognWL$H0$interim$a0$RE,
      collected_est$prognW$H0$interim$a0$RE,
      collected_est$prognL$H0$interim$a0$RE,
      collected_est$prognnon$H0$interim$a0$RE)

cbind(collected_est$prognWL$H0$decision$a0$RE,
      collected_est$prognW$H0$decision$a0$RE,
      collected_est$prognL$H0$decision$a0$RE,
      collected_est$prognnon$H0$decision$a0$RE)

cbind(collected_est$prognWL$H0$interim$a1$RE,
      collected_est$prognW$H0$interim$a1$RE,
      collected_est$prognL$H0$interim$a1$RE,
      collected_est$prognnon$H0$interim$a1$RE)

cbind(collected_est$prognWL$H0$decision$a1$RE,
      collected_est$prognW$H0$decision$a1$RE,
      collected_est$prognL$H0$decision$a1$RE,
      collected_est$prognnon$H0$decision$a1$RE)



##### For Table E.2 in the paper: #####
# "RE from simulation" columns in the table

options(digits = 3)

cbind(collected_est$prognWL$H1$interim$a0$RE,
      collected_est$prognW$H1$interim$a0$RE,
      collected_est$prognL$H1$interim$a0$RE,
      collected_est$prognnon$H1$interim$a0$RE)

cbind(collected_est$prognWL$H1$decision$a0$RE,
      collected_est$prognW$H1$decision$a0$RE,
      collected_est$prognL$H1$decision$a0$RE,
      collected_est$prognnon$H1$decision$a0$RE)

cbind(collected_est$prognWL$H1$interim$a1$RE,
      collected_est$prognW$H1$interim$a1$RE,
      collected_est$prognL$H1$interim$a1$RE,
      collected_est$prognnon$H1$interim$a1$RE)

cbind(collected_est$prognWL$H1$decision$a1$RE,
      collected_est$prognW$H1$decision$a1$RE,
      collected_est$prognL$H1$decision$a1$RE,
      collected_est$prognnon$H1$decision$a1$RE)



# evaluate trials to get type I error, power, ESS, ED --------------------------------

# This is for Table 4 in the paper.

npara <- 10
RHO <- 2

source("source/v.R")
v$nsim <- nsim



# compute power for the temptation trials

scns <- c("unadj", "prognWL", "prognW", "prognL", "prognnon")
nmaxs <- c(480, 300, 300, 480,480)


# main table of nmax, type I error, power, expected sample size
output <- matrix(NA, nrow = length(scns), ncol = 5)
rownames(output) <- scns
colnames(output) <- c("nmax", "type I error", "power", "ESS_H0", "ESS_H1")
output[, "nmax"] <- nmaxs

# table of percentage early stop at each analysis
early_stop_eff <- list(NA, NA, NA, NA, NA)
names(early_stop_eff) <- paste0(scns, " nmax", nmaxs)
early_stop_fut <- early_stop_eff

# table of reversal probabilities
rv_prob <- list(NA, NA, NA, NA, NA)
names(rv_prob) <- paste0(scns, " nmax", nmaxs)

for (iscn in 1:length(scns)){
  scn <- scns[iscn]
  nmax <- nmaxs[iscn]
  v$nmax <- nmax
  
  cat("\n###########################################\n")
  cat(sprintf("Current scenario: %s. Current nmax: %d.\n\n", scn, nmax))
  cat("Compute type I error and power:\n")
  
  taskIDs <- 1:npara
  if (scn == "unadj") { # use simulated unadj estimator from prognnon setting (unadj under all four setting should be the same)
    filenames_H0 <- paste0("result/est_prognnon_nmax", nmax, "_H0_", taskIDs, ".RDS")
    filenames_H1 <- paste0("result/est_prognnon_nmax", nmax, "_H1_", taskIDs, ".RDS")
  } else {
    filenames_H0 <- paste0("result/est_", scn, "_nmax", nmax, "_H0_", taskIDs, ".RDS")
    filenames_H1 <- paste0("result/est_", scn, "_nmax", nmax, "_H1_", taskIDs, ".RDS")
  }
  
  
  power <- comppower(v, filenames_H0, filenames_H1, rho = RHO, nonbinding = FALSE)
  if (scn == "unadj"){
    output[iscn, c("type I error", "power")] <- power[3:4]
  } else {
    output[iscn, c("type I error", "power")] <- power[1:2]
  }
  
  if (scn == "unadj"){
    tmp <- percentage_stopping(v, filenames_H0, filenames_H1, nonbinding = FALSE, rho = 2,
                               estimator = "unadj", scn = scn)
  } else {
    tmp <- percentage_stopping(v, filenames_H0, filenames_H1, nonbinding = FALSE, rho = 2,
                               estimator = "ltmle", scn = scn)
  }
  
  output[iscn, c("ESS_H0", "ESS_H1")] <- tmp$ESS[1, 2:3]
  early_stop_eff[[iscn]] <- tmp$stop_table_eff
  early_stop_fut[[iscn]] <- tmp$stop_table_fut
  
  rv_prob[[iscn]] <- reverse_prob(v, filenames_H0, filenames_H1)
}

# dir.create("log", showWarnings = FALSE)
# save(output, early_stop, rv_prob, file = "log/evalresult.rda")
dir.create("table", showWarnings = FALSE)
write.csv(output, file = "table/table_typeIerror_power_ESS_nsim=50000.csv")

# Table 4 in the paper:
# > output
#          nmax type I error power ESS_H0 ESS_H1
# unadj     480       0.0250 0.811    318    382
# prognWL   300       0.0254 0.791    225    259
# prognW    300       0.0256 0.805    227    260
# prognL    480       0.0253 0.805    309    375
# prognnon  480       0.0248 0.811    321    384


# For Table D.2 in the appendix ---------------------------------------------------------------------

# The following is printed while running the above loop.

###########################################
# Current scenario: unadj. Current nmax: 480.
# 
# Compute type I error and power:
#   
#   ltmle: errors and boundaries at each stage:
#   alpha   0.0009   0.0029   0.0051   0.0069   0.0093
# beta    0.0073   0.0229   0.0406   0.0549   0.0743
# u       3.12   2.73   2.47   2.28   Inf
# c       1.24   1.50   1.72   1.91   2.04
# l       -1.16   -0.08   0.72   1.37   -Inf
# 
# unadj: errors and boundaries at each stage:
#   alpha   0.0010   0.0029   0.0051   0.0069   0.0092
# beta    0.0078   0.0232   0.0406   0.0548   0.0736
# u       3.10   2.72   2.47   2.28   Inf
# c       1.24   1.50   1.73   1.91   2.04
# l       -1.11   -0.06   0.74   1.39   -Inf
# ####################################
# scnario: unadj unadj, nonbinding = FALSE, rho = 2.
# nmax ESS_H0 ESS_H1
# sample size  480    318    382
# 
# ###########################################
# Current scenario: prognWL. Current nmax: 300.
# 
# Compute type I error and power:
#   
#   ltmle: errors and boundaries at each stage:
#   alpha   0.0009   0.0032   0.0051   0.0072   0.0086
# beta    0.0075   0.0254   0.0411   0.0573   0.0688
# u       3.11   2.71   2.47   2.27   Inf
# c       1.30   1.54   1.74   1.91   2.07
# l       -1.20   -0.08   0.70   1.36   -Inf
# 
# unadj: errors and boundaries at each stage:
#   alpha   0.0010   0.0030   0.0049   0.0069   0.0092
# beta    0.0078   0.0241   0.0396   0.0549   0.0736
# u       3.10   2.72   2.48   2.29   Inf
# c       1.18   1.36   1.53   1.67   2.11
# l       -1.37   -0.43   0.24   0.82   -Inf
# ####################################
# scnario: ltmle prognWL, nonbinding = FALSE, rho = 2.
# nmax ESS_H0 ESS_H1
# sample size  300    225    259
# 
# ###########################################
# Current scenario: prognW. Current nmax: 300.
# 
# Compute type I error and power:
#   
#   ltmle: errors and boundaries at each stage:
#   alpha   0.0008   0.0029   0.0050   0.0070   0.0093
# beta    0.0065   0.0234   0.0401   0.0558   0.0741
# u       3.16   2.73   2.48   2.28   Inf
# c       1.32   1.55   1.75   1.91   2.06
# l       -1.26   -0.14   0.66   1.33   -Inf
# 
# unadj: errors and boundaries at each stage:
#   alpha   0.0010   0.0030   0.0050   0.0072   0.0088
# beta    0.0077   0.0241   0.0403   0.0573   0.0706
# u       3.10   2.71   2.47   2.27   Inf
# c       1.17   1.36   1.53   1.68   2.11
# l       -1.38   -0.43   0.26   0.86   -Inf
# ####################################
# scnario: ltmle prognW, nonbinding = FALSE, rho = 2.
# nmax ESS_H0 ESS_H1
# sample size  300    227    260
# 
# ###########################################
# Current scenario: prognL. Current nmax: 480.
# 
# Compute type I error and power:
#   
#   ltmle: errors and boundaries at each stage:
#   alpha   0.0012   0.0032   0.0053   0.0075   0.0078
# beta    0.0093   0.0258   0.0426   0.0602   0.0621
# u       3.05   2.68   2.44   2.24   Inf
# c       1.24   1.50   1.72   1.91   2.06
# l       -0.99   0.04   0.80   1.47   -Inf
# 
# unadj: errors and boundaries at each stage:
#   alpha   0.0010   0.0029   0.0050   0.0070   0.0091
# beta    0.0078   0.0233   0.0400   0.0562   0.0727
# u       3.10   2.72   2.47   2.27   Inf
# c       1.24   1.50   1.72   1.90   2.05
# l       -1.11   -0.07   0.71   1.38   -Inf
# ####################################
# scnario: ltmle prognL, nonbinding = FALSE, rho = 2.
# nmax ESS_H0 ESS_H1
# sample size  480    309    375
# 
# ###########################################
# Current scenario: prognnon. Current nmax: 480.
# 
# Compute type I error and power:
#   
#   ltmle: errors and boundaries at each stage:
#   alpha   0.0009   0.0029   0.0051   0.0069   0.0093
# beta    0.0073   0.0229   0.0406   0.0549   0.0743
# u       3.12   2.73   2.47   2.28   Inf
# c       1.24   1.50   1.72   1.91   2.04
# l       -1.16   -0.08   0.72   1.37   -Inf
# 
# unadj: errors and boundaries at each stage:
#   alpha   0.0010   0.0029   0.0051   0.0069   0.0092
# beta    0.0078   0.0232   0.0406   0.0548   0.0736
# u       3.10   2.72   2.47   2.28   Inf
# c       1.24   1.50   1.73   1.91   2.04
# l       -1.11   -0.06   0.74   1.39   -Inf
# ####################################
# scnario: ltmle prognnon, nonbinding = FALSE, rho = 2.
# nmax ESS_H0 ESS_H1
# sample size  480    321    384

