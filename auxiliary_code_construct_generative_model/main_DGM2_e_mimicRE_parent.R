rm(list = ls())

library(ltmle)

setwd("MISTIE/code")

source("source/fcn_DGM.R")

pY <- 0.02955846 # from main.DGM1_pY_mimicTE.R
ODatRE <- read.csv("log/ODatRE_nmax500.csv")[, 1] # from ODate_RE_5stg.csv

e <- read.csv("DGM2_e.csv")[,1]

left_e <- 0
right_e <- 0.3

Npara <- 500
iter <- read.csv("DGM2_iter.csv")[,1]
finished <- 0

error_e <- 0.001
# want to find e within an error of 0.01, such that augmented data mimics the relative efficiency
error_RE <- 0.01
# or to find e such that sum(diff(RE)) < 0.1

while ( finished == 0 ){ # terminating flag
  
  resultnames <- paste0("Result/DGM2_e_mimicRE_", iter, "_", 1:Npara, ".rda") # result filenames from parallel jobs
  
  ### wait until parallel jobs are done ###
  while ( sum(file.exists(resultnames)) < Npara){
    Sys.sleep(30)
    print(paste0(iter, " iteration, wait for parallel jobs to be done: ", sum(file.exists(resultnames))))
  }
  
  ### read all the results from 50 parallel jobs ###
  mat_ltmle <- matrix(NA, nrow = 1, ncol = 5)
  mat_unadj <- matrix(NA, nrow = 1, ncol = 5)  
  for(ifile in 1:Npara){
    load(resultnames[ifile])
    mat_ltmle <- rbind(mat_ltmle, test$ltmle_est)
    mat_unadj <- rbind(mat_unadj, test$unadj_est)
  }
  mat_ltmle <- mat_ltmle[-1, ]
  mat_unadj <- mat_unadj[-1, ]
  
  ### an operation of the result: compute relative efficiency ###
  var_ltmle <- apply(mat_ltmle, 2, var)
  var_unadj <- apply(mat_unadj, 2, var)
  
  RE <- var_unadj / var_ltmle # relative efficiency vector (stage 1 through stage 5)
  
  diffRE <- RE - ODatRE
  
  ### print information for current stage:
  print("########################################################")
  print(paste0(iter, " iteration done."))
  print(paste("Current e value:", left_e, e, right_e))
  print(paste("ODatRE:", paste(round(ODatRE,3), collapse = " ")))
  print(paste("Current RE:", paste(round(RE,3), collapse = " ")))
  
  ### binary search: update left and right thresholds, and new e ###
  
  if (abs(right_e - left_e) < error_e | abs(sum(diffRE)) < error_RE) { # found!
    write.table(paste0("nmax = 500, enrollrate = 100, e is found to be ", e), file = "parameter_found.txt",
                row.names = FALSE, col.names = FALSE)
    print(paste0("!!!DONE!!! e is found to be ", e))
    finished <- 1
  } else if (sum(diffRE) > 0) { # need more noise in Y (i.e. bigger e)
    left_e <- e
    e <- (left_e+right_e)/2
  } else if (sum(diffRE) < 0) { # need less noise in Y (i.e. smaller e)
    right_e <- e
    e <- (left_e+right_e)/2
  }
  
  ### print information for new value:
  print(paste("New e value:", left_e, e, right_e))
  print("########################################################")
  
  write.csv(e, "DGM2_e.csv", row.names = FALSE) # to communicate with parallel jobs
  
  ### iteration number +1 ###
  iter <- iter + 1
  write.csv(iter, "DGM2_iter.csv", row.names = FALSE) # to communicate with parallel jobs
  
  ### Create a file to indicate finish of one loop ###
  # to be used in bash shell script
  write.table(paste0("Current value:", left_e, e, right_e), file = "DGM2_single_loop_done.txt",
              row.names = FALSE, col.names = FALSE)
}
