rm(list = ls())

library(ltmle)

setwd("MISTIE/code_determine_nmax/exW")

source("fcn_DGM.R")
source("fcn_simtrials.R")

version_fut <- as.integer(Sys.getenv("SGE_TASK_ID"))

set.seed(version_fut)

dt <- read.csv("data100pts(updated).csv")

iter <- read.csv("iter.csv")[, 1]
nmax <- read.csv("nmax.csv")[ ,1]

v <- list(nsim = 100,
          nmax = nmax,
          K = 5,
          enrollrate = 140,
          pY = 0.02955846, # calibration probability
          dir = "up", # for calibration
          e = 0.1422, # calibration probability
          p0 = 0.2222222, # for calibration, p0 = P(Y=1|A=0)
          p1 = 0.34375, # for calibration, p1 = P(Y=1|A=1)
          randomA = FALSE,
          exL = FALSE,
          exW = TRUE,
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

test <- simtrials(dt, v)

save(test, file = paste0("Result/nmax_H1exW_enrollrate", v$enrollrate, "_", iter, "_", version_fut, ".rda"))

v$randomA <- TRUE

test <- simtrials(dt, v)
save(test, file = paste0("Result/nmax_H0exW_enrollrate", v$enrollrate, "_", iter, "_", version_fut, ".rda"))