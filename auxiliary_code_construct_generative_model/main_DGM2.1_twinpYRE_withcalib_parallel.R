rm(list = ls())

library(ltmle)

setwd("MISTIE/code_DGM")

source("fcn_DGM.R")
source("fcn_simtrials.R")

version_fut <- as.integer(Sys.getenv("SGE_TASK_ID"))

set.seed(version_fut)

dt <- read.csv("data100pts(updated).csv")

pY <- 0.02955846  # from main.DGM1_pY_mimicTE.R
direction <- "up" # from main.DGM1_pY_mimicTE.R

e <- 0.1422
iter <- 0.1422

###### WRONG!!!!! #####
## p0 and p1 should be computed in a more complicated way?
p0 <- 0.2222222 # = mean(subset(dt100, A == 0)$Y) = P(Y=1|A=0) 
p1 <- 0.34375   # = mean(subset(dt100, A == 1)$Y) = P(Y=1|A=1)

v <- list(nsim = 1000,
          nmax = 434,
          K = 5,
          enrollrate = 140,
          pY = pY, # calibration probability
          dir = direction, # for calibration
          e = e, # calibration probability
          p0 = p0, # for calibration, p0 = P(Y=1|A=0)
          p1 = p1, # for calibration, p1 = P(Y=1|A=1)
          randomA = FALSE,
          exL = FALSE,
          exW = FALSE,
          A_to_L1 = 30 / 365,
          L1_to_Y = (180 - 30) / 365,
          alpha = 0.025,
          beta = 0.2,
          rho = 2,
          f_err = 
               function(x, .rho = rho, .alpha = alpha){
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
               function(x, .rho = rho, .beta = beta){
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

test <- simtrials.nointerim(dt, v)

save(test, file = paste0("Result_RE/result_RE_", iter, "_", version_fut, ".rda"))