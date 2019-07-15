rm(list = ls())

library(ltmle)

source("fcn_DGM.R")

dt100 <- read.csv("data100pts(updated).csv")

dt <- dt100

## Generate twin data set

dtaug <- DGM.step0.twins(dt)

### Now, dt is the original data set, dtaug is the augmented data set

## compute trt effect for dt and dtaug

# without calibration, dtaug has larger treatment effect than dt
compute.unadj(dt) # 0.1215278
compute.ltmle(dt) # 0.1087564
compute.unadj(dtaug) # 0.11
compute.ltmle(dtaug) # 0.11


set.seed(123)

### Use root finding algorithm to find pY, to mimic treatment effect
find.pY(dtaug, nsim = 50000, true_trteff = 0.1215278, direction = "up")
# 0.02955846

# # look at combined
# nsim <- 10000
# dtaug_unadj <- rep(NA, nsim)
# 
# for (i in 1:nsim){
#   dtaug.cal <- calibrate.TE(dtaug, pY = 0.02955846, direction = "up")
#   
#   dtaug_unadj[i] <- compute.unadj(dtaug.cal)
# }
# 
# mean(dtaug_unadj)
# # 0.121467