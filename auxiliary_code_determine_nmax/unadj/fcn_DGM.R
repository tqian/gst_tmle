# wAge

compare.corr <- function(dt){
  fit_muhatA <- glm(Y ~ A, data = dt, family = "binomial")
  muhatA <- predict(fit_muhatA, type = "response")
  
  fit_muhatAW <- glm(Y ~ A + W1 + W2, data = dt, family = "binomial")
  muhatAW <- predict(fit_muhatAW, type = "response")
  
  fit_muhatAL <- glm(Y ~ A + L1, data = dt, family = "binomial")
  muhatAL <- predict(fit_muhatAL, type = "response")
  
  fit_muhatAWL <- glm(Y ~ A + W1 + W2 + L1, data = dt, family = "binomial")
  muhatAWL <- predict(fit_muhatAWL, type = "response")
  
  Y <- dt$Y
  
  cor_W <- sum((Y-muhatA)^2) / sum((Y-muhatAW)^2)
  cor_L <- sum((Y-muhatA)^2) / sum((Y-muhatAL)^2)
  cor_L_given_W <- sum((Y-muhatAW)^2) / sum((Y-muhatAWL)^2)
  
  return(list(cor_W = cor_W,
              cor_L = cor_L,
              cor_L_given_W = cor_L_given_W))
}



compute.unadj <- function(dt){
  # compute unadj estimator for a data frame
  unadj <- mean(subset(dt, A == 1)$Y) - mean(subset(dt, A == 0)$Y)
  return(unadj)
}



compute.ltmle <- function(dt){  
  # compute ltmle estimator for a data frame
  dttmle <- data.frame(S = dt$S, W1 = dt$W1, W2 = dt$W2, W3 = dt$W3, A = dt$A, L1 = dt$L1, Y = dt$Y)
  suppressMessages(ltmle.fit1 <- ltmle(dttmle, Anodes = "A", Lnodes = "L1", Ynodes = "Y", abar = 1, estimate.time=FALSE))
  suppressMessages(ltmle.fit0 <- ltmle(dttmle, Anodes = "A", Lnodes = "L1", Ynodes = "Y", abar = 0, estimate.time=FALSE))  
  ltmle <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]  
  return(ltmle)
}



relabel <- function(dat, pL1 = NULL, pY = NULL){
  # Generic function. Used in calibrate().
  # Reset L1 with probability pL
  # Reset Y with probability pY
  
  idxt <- which(dat$A == 1)
  idxc <- which(dat$A == 0)
  dat$L1[idxt] <- ifelse(rbinom(length(idxt),1,pL) == 1, 0, dat$L1[idxt])
  dat$L1[idxc] <- ifelse(rbinom(length(idxc),1,pL) == 1, 1, dat$L1[idxc])
  dat$Y[idxt]  <- ifelse(rbinom(length(idxt),1,pY) == 1, 0, dat$Y[idxt])
  dat$Y[idxc]  <- ifelse(rbinom(length(idxc),1,pY) == 1, 1, dat$Y[idxc])
  
  return(dat)
}



calibrate <- function(dt, # relabel probabilities
                      p1L1 = 0.07, p1Y = 0.07, p2L1 = 0.08, p2Y = 0.08, eps1 = 0.01,eps2 = 0.2){
  
  dt1 <- subset(dt, S == 1)
  dt2 <- subset(dt, S == 2)
  
  ## relabel L1 and Y to modify treatment effect
  idtwin1 <- which(dt1$twin == 1)
  dt1[idtwin1, ] <- relabel(dt1[idtwin1, ], p1L1, p1Y)  
  
  idtwin2 <- which(dt2$twin == 1)
  dt2[idtwin2, ] <- relabel(dt2[idtwin2, ], p2L1, p2Y)
  
  ## relabel Y to modify relative efficiency (doesn't change the treatment effect)
  pr1 <- ifelse(dt1$A == 1, mYt1, mYc1)
  rY1 <- rbinom(length(pr1), 1, pr1)
  dt1$Y <- ifelse(rbinom(nrow(dt1), 1, eps1) == 1, rY1, dt1$Y)
  
  pr2 <- ifelse(dt2$A == 1, mYt2, mYc2)
  rY2 <- rbinom(length(pr2), 1, pr2)
  dt2$Y <- ifelse(rbinom(nrow(dt2), 1, eps2) == 1, rY2, dt2$Y)
  
  return(rbind(dt1, dt2))
}

calibrate.TE <- function(dt, pY, direction = c("up", "down")){
  # pY = 0.02228878 for updated 208 patients
  
  dir <- direction
  
  ## relabel Y to modify treatment effect
  idtwin <- which(dt$twin == 1)
  dt[idtwin, ] <- relabel.pY(dt[idtwin, ], pY, direction = dir)  
  
  return(dt)
}

relabel.pY <- function(dt, pY, direction = c("up", "down")){
  # Generic function. Used in calibrate().
  # Reset Y with probability pY
  
  idA1 <- which(dt$A == 1)
  idA0 <- which(dt$A == 0)
  
  if (direction == "up"){
    # make more treated Y to be 1
    dt$Y[idA1]  <- ifelse(rbinom(length(idA1),1,pY) == 1, 1, dt$Y[idA1])
    dt$Y[idA0]  <- ifelse(rbinom(length(idA0),1,pY) == 1, 0, dt$Y[idA0])
  } else if (direction == "down"){
    # make more treated Y to be 0
    dt$Y[idA1]  <- ifelse(rbinom(length(idA1),1,pY) == 1, 0, dt$Y[idA1])
    dt$Y[idA0]  <- ifelse(rbinom(length(idA0),1,pY) == 1, 1, dt$Y[idA0])
  }
  
  return(dt)
}


### Use root finding algorithm to find pY, to mimic treatment effect
find.pY <- function(dtaug, nsim, true_trteff, direction = c("up", "down")){
  # direction: "up" means to increase the trteff of dtaug to match trt_trteff,
  #            "down" means to decrease the trteff of dtaug to match trt_trteff
  
  idtwin <- which(dtaug$twin == 1)  
  dir <- direction
  
  f <- function(pY){
    unadj <- rep(NA, nsim)
    for (i in 1:nsim){
      dtaug_tmp <- dtaug
      dtaug_tmp[idtwin, ] <- relabel.pY(dtaug_tmp[idtwin, ], pY, direction = dir)
      unadj[i] <- compute.unadj(dtaug_tmp)
    }
    return(mean(unadj) - true_trteff)
  }  
  result <- uniroot(f, c(0,1))  
  return(result)
}


DGM.step0.twins <- function(dt){
  fitL1 <- glm(L1 ~ S + A + W1 + W2 + W3,           family = binomial, data = dt)
  fitL2 <- glm(L2 ~ S + A + W1 + W2 + W3 + L1,      family = binomial, data = dt)
  fitY  <- glm(Y  ~ S + A + W1 + W2 + W3 + L1 + L2, family = binomial, data = dt)
  
  ## Generate twin data set
  
  dttwin <- dt
  dttwin$A <- ifelse(dt$A == 1, 0, 1) # swap treatment assigned
  
  prL1 <- predict(fitL1, newdata = dttwin, type="response")
  prL2 <- predict(fitL2, newdata = dttwin, type="response")
  prY  <- predict(fitY,  newdata = dttwin, type="response")
  
  dttwin$L1 <- as.numeric(prL1 > 0.5)
  dttwin$L2 <- as.numeric(prL2 > 0.5)
  dttwin$Y  <- as.numeric(prY  > 0.5)
  
  dtaug <- rbind(dt, dttwin)
  dtaug$twin <- c(rep(0, nrow(dt)), rep(1, nrow(dttwin))) # 1 from twin, 0 from original
  
  return(dtaug)
}


DGM.step1.pY.mimicTE <- function(dt, pY, direction = c("up", "down")){
  # pY = 0.02228878 for updated 208 patients
  
  dir <- direction
  
  ## relabel Y to modify treatment effect
  idtwin <- which(dt$twin == 1)
  dt[idtwin, ] <- relabel.pY(dt[idtwin, ], pY, direction = dir)  
  
  return(dt)
}


DGM.step2.e.mimicRE <- function(dt, e, p0, p1){
  
  id_resetY <- rbinom(nrow(dt), dt$twin == 1, e) # pick out newly added twins, wp e  
  p_resetY <- ifelse(dt$A == 0, p0, p1) # resetting probability for each person (if they are picked in id_resetY)
  
  resetY <- rbinom(sum(id_resetY), 1, p_resetY[which(id_resetY == 1)])
  dt[which(id_resetY == 1), "Y"] <- resetY
  
  return(dt)
}
