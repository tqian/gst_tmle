rm(list = ls())

library(ltmle)

setwd("MISTIE/code")

source("source/fcn_DGM.R")
source("source/fcn_simtrials.R")

taskID <- as.integer(Sys.getenv("SGE_TASK_ID"))

# read in parallel random seed
parallel_seeds <- data.matrix(read.csv("parallel_seed/parallel_seeds.csv"))
.Random.seed <- parallel_seeds[taskID, ]

dt <- read.csv("data/data100pts(updated).csv")

source("source/v.R")
v$nsim <- 500
v$nmax <- 500

test <- simtrials.Odat.nocalibration.nointerim(dt, v)

save(test, file = paste0("Result/ODatRE_nmax", v$nmax, "_", taskID, ".rda"))