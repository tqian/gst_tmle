rm(list = ls())

library(ltmle)

# Parallel setup ----------------------------------------------------------

on_HarvardRC <- TRUE

if (on_HarvardRC) {
  version_fut <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  setwd("~/gst_tmle")
} else {
  # on local computer
  setwd("~/Dropbox/Research/git_GST_TMLE")
}

source("source/fcn_simtrials.R")
source("source/v.R") # config of trial

dt100 <- read.csv("data/data100pts(updated).csv")

npara <- 10

parallel_seeds <- data.matrix(read.csv("parallel_seed/parallel_seeds.csv"))
taskID <- ((version_fut - 1) %% npara) + 1
case_id <- (version_fut - 1) %/% npara + 1


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


## Start simulation

# read in parallel random seed


v$nsim <- 5000
v$nmax <- nmax
v$enrollrate <- 140
v$exL <- exL
v$exW <- exW

print(paste("version_fut:", version_fut))
print(paste("taskID:", taskID))
print(paste("case_id:", case_id))
print(paste("nmax =", v$nmax))
print(paste("enrollrate =", v$enrollrate))
print(paste("exL", v$exL, "exW", v$exW))

if (exL & exW) {
  scn <- "prognnon"
} else if (exL) {
  scn <- "prognW"
} else if (exW) {
  scn <- "prognL"
} else {
  scn <- "prognWL"
}
v$scn <- scn

dir.create("result", showWarnings = FALSE)
filenameH0 <- paste0("result/est_", scn, "_nmax", v$nmax, "_H0_", taskID, ".RDS")
filenameH1 <- paste0("result/est_", scn, "_nmax", v$nmax, "_H1_", taskID, ".RDS")

# simulate under H0
.Random.seed <- parallel_seeds[taskID, ]
v$randomA <- TRUE # reset A randomly; i.e., H0
test <- simtrials(dt100, v, identical_L1Y = FALSE)  
test$v <- v
test$taskID <- taskID
test$random_seed <- .Random.seed
saveRDS(test, file = filenameH0)

# simulate under H1
.Random.seed <- parallel_seeds[taskID, ]
v$randomA <- FALSE # do not reset A; i.e., H1
test <- simtrials(dt100, v, identical_L1Y = FALSE)  
test$v <- v
test$taskID <- taskID
test$random_seed <- .Random.seed
saveRDS(test, file = filenameH1)