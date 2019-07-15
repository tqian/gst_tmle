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
v$nsim <- 100
v$nmax <- 500

e <- read.csv("DGM2_e.csv")[, 1]
iter <- read.csv("DGM2_iter.csv")[, 1]

## p0 and p1 should be computed in a more complicated way?
p0 <- 0.2222222 # = mean(subset(dt100, A == 0)$Y) = P(Y=1|A=0) 
p1 <- 0.34375   # = mean(subset(dt100, A == 1)$Y) = P(Y=1|A=1)

v$p0 <- p0
v$p1 <- p1
v$e <- e

test <- simtrials.nointerim(dt, v)

save(test, file = paste0("Result/DGM2_e_mimicRE_", iter, "_", taskID, ".rda"))