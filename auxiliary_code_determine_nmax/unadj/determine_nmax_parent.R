rm(list = ls())

library(ltmle)

setwd("MISTIE/code_determine_nmax")

source("fcn_nmax_unadj.R")

nameH0 <- "Result/nmax_H0unadj_enrollrate140"
nameH1 <- "Result/nmax_H1unadj_enrollrate140"

nmax <- read.csv("nmax.csv")[, 1]
v <- list(nsim = 50000,
          nmax = nmax,
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


left_nmax <- 400
right_nmax <- 500

Npara <- 50
iter <- read.csv("iter.csv")[, 1]
finished <- 0

error_nmax <- 2
# want to find nmax within an error of 4
error_power <- 0.001
# or to find nmax such that power is close to 0.80 with an error of 0.01

while ( finished == 0 ){ # terminating flag
  
  resultnamesH0 <- paste0(nameH0, "_", iter, "_", 1:Npara, ".rda") # result filenames from parallel jobs
  resultnamesH1 <- paste0(nameH1, "_", iter, "_", 1:Npara, ".rda") # result filenames from parallel jobs
  
  ### wait until parallel jobs are done ###
  while ( sum(file.exists(resultnamesH0)) < Npara | sum(file.exists(resultnamesH1)) < Npara){
    Sys.sleep(10)
    print(paste(iter, "iteration, wait for parallel jobs to be done: H0",
                 sum(file.exists(resultnamesH0)), ";H1", sum(file.exists(resultnamesH1))))
  }
  
  ### read all the results from Npara parallel jobs ###
  caseH0 <- paste0(nameH0, "_", iter)
  caseH1 <- paste0(nameH1, "_", iter)
  
  power <- comppower.unadjonly(v, caseH0, caseH1, 1:Npara, nonbinding = FALSE)[2]
  
  
  ### print information for current stage:
  print("########################################################")
  print(paste0(iter, " iteration done."))
  print(paste("Current nmax value:", left_nmax, nmax, right_nmax))
  print(paste("Power:", power))
  
  ### binary search: update left and right thresholds, and new e ###
  
  if (abs(right_nmax - left_nmax) <= error_nmax | abs(power - 0.8) < error_power) { # found!
    write.table(paste0("enrollrate = 140, unadj, nmax is found to be ", nmax), file = "nmax_found.txt",
                row.names = FALSE, col.names = FALSE)
    print(paste0("!!!DONE!!! nmax is found to be ", nmax))
    finished <- 1
  } else if ( power > 0.8) { # need smaller sample size
    right_nmax <- nmax
    nmax <- round((left_nmax + right_nmax) / 2)
  } else if ( power < 0.8) { # need larger sample size
    left_nmax <- nmax
    nmax <- round((left_nmax + right_nmax) / 2)
  }
  if (nmax%%2 == 1) { # if nmax is odd, make it even
    nmax <- nmax + 1
  }
  
  
  ### print information for new value:
  print(paste("New nmax value:", left_nmax, nmax, right_nmax))
  print("########################################################")
  
  write.csv(nmax, "nmax.csv", row.names = FALSE) # to communicate with parallel jobs
  
  ### iteration number +1 ###
  iter <- iter + 1
  write.csv(iter, "iter.csv", row.names = FALSE) # to communicate with parallel jobs
  
  ### Create a file to indicate finish of one loop ###
  # to be used in bash shell script
  write.table(paste0("Current value:", left_nmax, nmax, right_nmax), file = "nmax_single_loop_done.txt",
              row.names = FALSE, col.names = FALSE)
}
