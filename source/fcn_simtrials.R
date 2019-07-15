library(ltmle)
source("source/fcn_DGM.R")




simtrials <- function(dt, v, identical_L1Y = FALSE, only_unadj = FALSE){
    # identical_L1Y: make L1 and Y identical, for debugging purpose
    
    # parameters set-up
    nsim <- v$nsim
    nmax <- v$nmax
    K <- v$K
    enrollrate <- v$enrollrate
    randomA <- v$randomA
    exL <- v$exL
    exW <- v$exW
    A_to_L1 <- v$A_to_L1
    L1_to_Y <- v$L1_to_Y
    interim_times  <- (1:K)/K * nmax/enrollrate + A_to_L1 + L1_to_Y # so that we have (k/K)*nmax patients with final outcome at k-th inteirm analysis
    decision_times <- pmin(interim_times + A_to_L1 + L1_to_Y, interim_times[K])
    # at K-th interim analysis, we have everyone's final outcome observed.
    
    # output initializaiton
    output <- list(ltmle_est = matrix(NA, nrow = nsim, ncol = 2*K),
                   ltmle_est_1 = matrix(NA, nrow = nsim, ncol = 2*K),
                   ltmle_est_0 = matrix(NA, nrow = nsim, ncol = 2*K),                 
                   unadj_est = matrix(NA, nrow = nsim, ncol = 2*K),
                   unadj_est_1 = matrix(NA, nrow = nsim, ncol = 2*K),
                   unadj_est_0 = matrix(NA, nrow = nsim, ncol = 2*K),
                   randomA = randomA,
                   exL = exL,
                   exW = exW)
    
    # Step 0: generate twins
    dtaug <- DGM.step0.twins(dt)  
    
    
    for( itrial in 1:nsim ){
        
        if( itrial %% 10 == 0){
            cat(paste0("\nTrial number ", itrial))
        }
        
        # Sample from the twin data with replacement, to have desired sample size
        ids <- sample(nrow(dtaug), nmax, replace = TRUE)
        dtc <- dtaug[ids, ]
        
        # Step 1: Calibrate Y with pY, to mimic Treatment Effect
        dtc <- DGM.step1.pY.mimicTE(dtc, pY = v$pY, direction = v$dir)
        
        # Step 2: Calibrate Y with e, to mimic Relative Efficiency
        dtc <- DGM.step2.e.mimicRE(dtc, e = v$e, p0 = v$p0, p1 = v$p1)
        
        if (identical_L1Y) {
            dtc$L1 <- dtc$Y
        }
        
        if (randomA) { # Reset A as Bernoulli(0.5), i.e. assume no effect
            dtc$A <- rbinom(nrow(dtc), 1, 0.5)      
        }
        if (exL) { # Set L exogeneous (complete permute)
            dtc$L1 <- sample(dtc$L1, nrow(dtc), replace = FALSE)
            dtc$L2 <- sample(dtc$L2, nrow(dtc), replace = FALSE)
        }
        if (exW) { # Set W exogeneous (complete permute)
            dtc$S <- sample(dtc$S, nrow(dtc), replace = FALSE)
            dtc$W1 <- sample(dtc$W1, nrow(dtc), replace = FALSE)
            dtc$W2 <- sample(dtc$W2, nrow(dtc), replace = FALSE)
            dtc$W3 <- sample(dtc$W3, nrow(dtc), replace = FALSE)
        }
        
        
        
        # Assign enrollment time for each patient
        
        dtc$enrolltime <- sample(1:nmax, nmax) / enrollrate
        dtc$L1time <- dtc$enrolltime + A_to_L1
        dtc$Ytime  <- dtc$L1time + L1_to_Y
        
        for( istage in 1:K ){ # Do each interim analysis and decision analysis, compute the estimators
            
            ## Interim analysis
            
            intrtime <- interim_times[istage]
            
            C0 <- as.numeric(dtc$enrolltime <= intrtime)
            C1 <- as.numeric(dtc$L1time     <= intrtime)
            C2 <- as.numeric(dtc$Ytime      <= intrtime)
            
            dttmle <- data.frame(C0 = C0, S = dtc$S, W3 = dtc$W3, A = dtc$A,
                                 C1 = C1, L1 = dtc$L1,
                                 C2 = C2, Y = dtc$Y)
            
            if (!only_unadj) {
                suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
                                                     Anodes = "A",
                                                     Cnodes = c("C0","C1","C2"),
                                                     Lnodes = c("S","W3","L1"),
                                                     Ynodes = "Y",
                                                     abar = 1,
                                                     estimate.time=FALSE))
                suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
                                                     Anodes = "A",
                                                     Cnodes = c("C0","C1","C2"),
                                                     Lnodes = c("S","W3","L1"),
                                                     Ynodes = "Y",
                                                     abar = 0,
                                                     estimate.time=FALSE))
                
                ltmle_fit1 <- ltmle.fit1$estimates[1]
                ltmle_fit0 <- ltmle.fit0$estimates[1]
                output$ltmle_est[itrial, 2*istage - 1] <- ltmle_fit1 - ltmle_fit0
                output$ltmle_est_1[itrial, 2*istage - 1] <- ltmle_fit1
                output$ltmle_est_0[itrial, 2*istage - 1] <- ltmle_fit0
            }
            
            unadj_fit1 <- mean(subset(dttmle, (C2 == 1) & (A == 1))$Y)
            unadj_fit0 <- mean(subset(dttmle, (C2 == 1) & (A == 0))$Y)
            output$unadj_est[itrial, 2*istage - 1] <- unadj_fit1 - unadj_fit0
            output$unadj_est_1[itrial, 2*istage - 1] <- unadj_fit1
            output$unadj_est_0[itrial, 2*istage - 1] <- unadj_fit0
            
            ## Decision analysis
            
            idx_decision <- which(C0==1) # decision analysis should be based on these people
            
            dt_decision <- dtc[idx_decision, ]
            
            dttmle <- data.frame(S = dt_decision$S, W3 = dt_decision$W3,
                                 A = dt_decision$A, Y = dt_decision$Y)
            
            if (!only_unadj) {
                suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
                                                     Anodes = "A",
                                                     Ynodes = "Y",
                                                     abar = 1,
                                                     estimate.time=FALSE))
                suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
                                                     Anodes = "A",
                                                     Ynodes = "Y",
                                                     abar = 0,
                                                     estimate.time=FALSE))
                
                ltmle_fit1 <- ltmle.fit1$estimates[1]
                ltmle_fit0 <- ltmle.fit0$estimates[1]
                output$ltmle_est[itrial, 2*istage] <- ltmle_fit1 - ltmle_fit0
                output$ltmle_est_1[itrial, 2*istage] <- ltmle_fit1
                output$ltmle_est_0[itrial, 2*istage] <- ltmle_fit0
            }
            
            unadj_fit1 <- mean(subset(dttmle, A == 1)$Y)
            unadj_fit0 <- mean(subset(dttmle, A == 0)$Y)
            output$unadj_est[itrial, 2*istage] <- unadj_fit1 - unadj_fit0
            output$unadj_est_1[itrial, 2*istage] <- unadj_fit1
            output$unadj_est_0[itrial, 2*istage] <- unadj_fit0
            
            
        } # end for loop: istage
    } # end for loop: itrial
    
    return(output)
}





simtrials_diffonly <- function(dt, v, identical_L1Y = FALSE){
  # identical_L1Y: make L1 and Y identical, for debugging purpose
  
  # parameters set-up
  nsim <- v$nsim
  nmax <- v$nmax
  K <- v$K
  enrollrate <- v$enrollrate
  randomA <- v$randomA
  exL <- v$exL
  exW <- v$exW
  A_to_L1 <- v$A_to_L1
  L1_to_Y <- v$L1_to_Y
  interim_times  <- (1:K)/K * nmax/enrollrate + A_to_L1 + L1_to_Y # so that we have (k/K)*nmax patients with final outcome at k-th inteirm analysis
  decision_times <- pmin(interim_times + A_to_L1 + L1_to_Y, interim_times[K])
  # at K-th interim analysis, we have everyone's final outcome observed.
  
  # output initializaiton
  output <- list(ltmle_est = matrix(NA, nrow = nsim, ncol = 2*K),  
                 unadj_est = matrix(NA, nrow = nsim, ncol = 2*K),
                 randomA = randomA,
                 exL = exL,
                 exW = exW)
  
  # Step 0: generate twins
  dtaug <- DGM.step0.twins(dt)  
  
  
  for( itrial in 1:nsim ){
    
    if( itrial %% 10 == 0){
      cat(paste0("\nTrial number ", itrial))
    }
    
    # Sample from the twin data with replacement, to have desired sample size
    ids <- sample(nrow(dtaug), nmax, replace = TRUE)
    dtc <- dtaug[ids, ]
    
    # Step 1: Calibrate Y with pY, to mimic Treatment Effect
    dtc <- DGM.step1.pY.mimicTE(dtc, pY = v$pY, direction = v$dir)
    
    # Step 2: Calibrate Y with e, to mimic Relative Efficiency
    dtc <- DGM.step2.e.mimicRE(dtc, e = v$e, p0 = v$p0, p1 = v$p1)
    
    if (identical_L1Y) {
      dtc$L1 <- dtc$Y
    }

    if (randomA == TRUE) { # Reset A as Bernoulli(0.5), i.e. assume no effect
      dtc$A <- rbinom(nrow(dtc), 1, 0.5)      
    }
    if (exL == TRUE) { # Set L exogeneous (complete permute)
      dtc$L1 <- sample(dtc$L1, nrow(dtc), replace = FALSE)
      dtc$L2 <- sample(dtc$L2, nrow(dtc), replace = FALSE)
    }
    if (exW == TRUE) { # Set W exogeneous (complete permute)
      dtc$S <- sample(dtc$S, nrow(dtc), replace = FALSE)
      dtc$W1 <- sample(dtc$W1, nrow(dtc), replace = FALSE)
      dtc$W2 <- sample(dtc$W2, nrow(dtc), replace = FALSE)
      dtc$W3 <- sample(dtc$W3, nrow(dtc), replace = FALSE)
    }
    

    
    # Assign enrollment time for each patient
    
    dtc$enrolltime <- sample(1:nmax, nmax) / enrollrate
    dtc$L1time <- dtc$enrolltime + A_to_L1
    dtc$Ytime  <- dtc$L1time + L1_to_Y
    
    for( istage in 1:K ){ # Do each interim analysis and decision analysis, compute the estimators
      
      ## Interim analysis
      
      intrtime <- interim_times[istage]
      
      C0 <- as.numeric(dtc$enrolltime <= intrtime)
      C1 <- as.numeric(dtc$L1time     <= intrtime)
      C2 <- as.numeric(dtc$Ytime      <= intrtime)
      
      dttmle <- data.frame(C0 = C0, S = dtc$S, W3 = dtc$W3, A = dtc$A,
                           C1 = C1, L1 = dtc$L1,
                           C2 = C2, Y = dtc$Y)
      
suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
                          Anodes = "A",
                          Cnodes = c("C0","C1","C2"),
                          Lnodes = c("S","W3","L1"),
                          Ynodes = "Y",
                          abar = 1,
                          estimate.time=FALSE))
suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
                          Anodes = "A",
                          Cnodes = c("C0","C1","C2"),
                          Lnodes = c("S","W3","L1"),
                          Ynodes = "Y",
                          abar = 0,
                          estimate.time=FALSE))
      
      output$ltmle_est[itrial, 2*istage - 1] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      output$unadj_est[itrial, 2*istage - 1] <- 
        mean(subset(dttmle, (C2 == 1) & (A == 1))$Y) - mean(subset(dttmle, (C2 == 1) & (A == 0))$Y)
      
      ## Decision analysis

      idx_decision <- which(C0==1) # decision analysis should be based on these people
      
      dt_decision <- dtc[idx_decision, ]
      
      dttmle <- data.frame(S = dt_decision$S, W3 = dt_decision$W3,
                           A = dt_decision$A, Y = dt_decision$Y)
      
suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
                                     Anodes = "A",
                                     Ynodes = "Y",
                                     abar = 1,
                                     estimate.time=FALSE))
suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
                                     Anodes = "A",
                                     Ynodes = "Y",
                                     abar = 0,
                                     estimate.time=FALSE))
      
      output$ltmle_est[itrial, 2*istage] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      output$unadj_est[itrial, 2*istage] <- 
        mean(subset(dttmle, A == 1)$Y) - mean(subset(dttmle, A == 0)$Y)
      
    } # end for loop: istage
  } # end for loop: itrial
  
  return(output)
}







simtrials.nointerim <- function(dt, v){
  # no interim analysis, for computing relative efficiency of decision analyses
  
  # parameters set-up
  nsim <- v$nsim
  nmax <- v$nmax
  K <- v$K
  enrollrate <- v$enrollrate
  randomA <- v$randomA
  exL <- v$exL
  exW <- v$exW
  A_to_L1 <- v$A_to_L1
  L1_to_Y <- v$L1_to_Y
  interim_times  <- (1:K)/K * nmax/enrollrate + A_to_L1 + L1_to_Y # so that we have (k/K)*nmax patients with final outcome at k-th inteirm analysis
  decision_times <- pmin(interim_times + A_to_L1 + L1_to_Y, interim_times[K])
  # at K-th interim analysis, we have everyone's final outcome observed.
  
  # output initializaiton
  output <- list(ltmle_est = matrix(NA, nrow = nsim, ncol = K),
                 unadj_est = matrix(NA, nrow = nsim, ncol = K),
                 randomA = randomA,
                 exL = exL,
                 exW = exW)
  
  # Step 0: generate twins
  dtaug <- DGM.step0.twins(dt)  
  
  
  for( itrial in 1:nsim ){
    
    if( itrial %% 10 == 0){
      cat(paste0("\nTrial number ", itrial))
    }
    
    # Sample from the twin data with replacement, to have desired sample size
    ids <- sample(nrow(dtaug), nmax, replace = TRUE)
    dtc <- dtaug[ids, ]
    
    # Step 1: Calibrate Y with pY, to mimic Treatment Effect
    dtc <- DGM.step1.pY.mimicTE(dtc, pY = v$pY, direction = v$dir)
    
    # Step 2: Calibrate Y with e0 and e1, to mimic Relative Efficiency
    dtc <- DGM.step2.e.mimicRE(dtc, e = v$e, p0 = v$p0, p1 = v$p1)
    
    #     if (randomA == TRUE) { # Reset A as Bernoulli(0.5), i.e. assume no effect
    #       dtc$A <- rbinom(nrow(dtc), 1, 0.5)      
    #     }
    #     if (exL == TRUE) { # Set L exogeneous (complete permute)
    #       dtc$L1 <- sample(dtc$L1, nrow(dtc), replace = FALSE)
    #     }
    #     if (exW == TRUE) { # Set W exogeneous (complete permute)
    #       dtc$W1 <- sample(dtc$W1, nrow(dtc), replace = FALSE)
    #       dtc$W2 <- sample(dtc$W2, nrow(dtc), replace = FALSE)
    #     }
    
    # Assign enrollment time for each patient
    
    dtc$enrolltime <- sample(1:nmax, nmax) / enrollrate
    dtc$L1time <- dtc$enrolltime + A_to_L1
    dtc$Ytime  <- dtc$L1time + L1_to_Y
    
    for( istage in 1:K ){ # Do each interim analysis and decision analysis, compute the estimators
      
            ## Interim analysis
            
            intrtime <- interim_times[istage]
            
            C0 <- as.numeric(dtc$enrolltime <= intrtime)
            C1 <- as.numeric(dtc$L1time     <= intrtime)
            C2 <- as.numeric(dtc$Ytime      <= intrtime)
      #       
      #       dttmle <- data.frame(C0 = C0, S = dtc$S, W1 = dtc$W1, W2 = dtc$W2, A = dtc$A,
      #                            C1 = C1, L1 = dtc$L1,
      #                            C2 = C2, Y = dtc$Y)
      #       
      # suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
      #                           Anodes = "A",
      #                           Cnodes = c("C0","C1","C2"),
      #                           Lnodes = c("S","W1","W2","L1"),
      #                           Ynodes = "Y",
      #                           abar = 1,
      #                           estimate.time=FALSE))
      # suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
      #                           Anodes = "A",
      #                           Cnodes = c("C0","C1","C2"),
      #                           Lnodes = c("S","W1","W2","L1"),
      #                           Ynodes = "Y",
      #                           abar = 0,
      #                           estimate.time=FALSE))
      #       
      #       output$ltmle_est[itrial, 2*istage - 1] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      #       output$unadj_est[itrial, 2*istage - 1] <- 
      #         mean(subset(dttmle, (C2 == 1) & (A == 1))$Y) - mean(subset(dttmle, (C2 == 1) & (A == 0))$Y)
      
      ## Decision analysis
      
      idx_decision <- which(C0==1) # decision analysis should be based on these people
      
      dt_decision <- dtc[idx_decision, ]
      
      dttmle <- data.frame(S = dt_decision$S, W3 = dt_decision$W3,
                           A = dt_decision$A, Y = dt_decision$Y)
      
      suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
                                           Anodes = "A",
                                           Ynodes = "Y",
                                           abar = 1,
                                           estimate.time=FALSE))
      suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
                                           Anodes = "A",
                                           Ynodes = "Y",
                                           abar = 0,
                                           estimate.time=FALSE))
      
      output$ltmle_est[itrial, istage] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      output$unadj_est[itrial, istage] <- 
        mean(subset(dttmle, A == 1)$Y) - mean(subset(dttmle, A == 0)$Y)
      
    } # end for loop: istage
  } # end for loop: itrial
  
  return(output)
}


simtrials.twinpY.nointerim <- function(dt, v){
  
  # parameters set-up
  nsim <- v$nsim
  nmax <- v$nmax
  K <- v$K
  enrollrate <- v$enrollrate
  randomA <- v$randomA
  exL <- v$exL
  exW <- v$exW
  A_to_L1 <- v$A_to_L1
  L1_to_Y <- v$L1_to_Y
  interim_times  <- (1:K)/K * nmax/enrollrate + A_to_L1 + L1_to_Y # so that we have (k/K)*nmax patients with final outcome at k-th inteirm analysis
  decision_times <- pmin(interim_times + A_to_L1 + L1_to_Y, interim_times[K])
  # at K-th interim analysis, we have everyone's final outcome observed.
  
  # output initializaiton
  output <- list(ltmle_est = matrix(NA, nrow = nsim, ncol = K),
                 unadj_est = matrix(NA, nrow = nsim, ncol = K),
                 randomA = randomA,
                 exL = exL,
                 exW = exW)
  
  # Step 0: generate twins
  dtaug <- DGM.step0.twins(dt)  
  
  
  for( itrial in 1:nsim ){
    
    if( itrial %% 10 == 0){
      cat(paste0("\nTrial number ", itrial))
    }
    
    # Sample from the twin data with replacement, to have desired sample size
    ids <- sample(nrow(dtaug), nmax, replace = TRUE)
    dtc <- dtaug[ids, ]
    
    # Step 1: Calibrate Y with pY, to mimic Treatment Effect
    dtc <- DGM.step1.pY.mimicTE(dtc, pY = v$pY, direction = v$dir)
    
#     # Step 2: Calibrate Y with e0 and e1, to mimic Relative Efficiency
#     dtc <- DGM.step2.e.mimicRE(dtc, e = v$e, p0 = v$p0, p1 = v$p1)
    
    #     if (randomA == TRUE) { # Reset A as Bernoulli(0.5), i.e. assume no effect
    #       dtc$A <- rbinom(nrow(dtc), 1, 0.5)      
    #     }
    #     if (exL == TRUE) { # Set L exogeneous (complete permute)
    #       dtc$L1 <- sample(dtc$L1, nrow(dtc), replace = FALSE)
    #     }
    #     if (exW == TRUE) { # Set W exogeneous (complete permute)
    #       dtc$W1 <- sample(dtc$W1, nrow(dtc), replace = FALSE)
    #       dtc$W2 <- sample(dtc$W2, nrow(dtc), replace = FALSE)
    #     }
    
    # Assign enrollment time for each patient
    
    dtc$enrolltime <- sample(1:nmax, nmax) / enrollrate
    dtc$L1time <- dtc$enrolltime + A_to_L1
    dtc$Ytime  <- dtc$L1time + L1_to_Y
    
    for( istage in 1:K ){ # Do each interim analysis and decision analysis, compute the estimators
      
            ## Interim analysis
            
            intrtime <- interim_times[istage]
            
            C0 <- as.numeric(dtc$enrolltime <= intrtime)
            C1 <- as.numeric(dtc$L1time     <= intrtime)
            C2 <- as.numeric(dtc$Ytime      <= intrtime)
      #       
      #       dttmle <- data.frame(C0 = C0, S = dtc$S, W1 = dtc$W1, W2 = dtc$W2, A = dtc$A,
      #                            C1 = C1, L1 = dtc$L1,
      #                            C2 = C2, Y = dtc$Y)
      #       
      # suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
      #                           Anodes = "A",
      #                           Cnodes = c("C0","C1","C2"),
      #                           Lnodes = c("S","W1","W2","L1"),
      #                           Ynodes = "Y",
      #                           abar = 1,
      #                           estimate.time=FALSE))
      # suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
      #                           Anodes = "A",
      #                           Cnodes = c("C0","C1","C2"),
      #                           Lnodes = c("S","W1","W2","L1"),
      #                           Ynodes = "Y",
      #                           abar = 0,
      #                           estimate.time=FALSE))
      #       
      #       output$ltmle_est[itrial, 2*istage - 1] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      #       output$unadj_est[itrial, 2*istage - 1] <- 
      #         mean(subset(dttmle, (C2 == 1) & (A == 1))$Y) - mean(subset(dttmle, (C2 == 1) & (A == 0))$Y)
      
      ## Decision analysis
      
      idx_decision <- which(C0==1) # decision analysis should be based on these people
      
      dt_decision <- dtc[idx_decision, ]
      
      dttmle <- data.frame(S = dt_decision$S, W3 = dt_decision$W3,
                           A = dt_decision$A, Y = dt_decision$Y)
      
      suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
                                           Anodes = "A",
                                           Ynodes = "Y",
                                           abar = 1,
                                           estimate.time=FALSE,
                                           variance.method = "ic"))
      suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
                                           Anodes = "A",
                                           Ynodes = "Y",
                                           abar = 0,
                                           estimate.time=FALSE,
                                           variance.method = "ic"))
      
      output$ltmle_est[itrial, istage] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      output$unadj_est[itrial, istage] <- 
        mean(subset(dttmle, A == 1)$Y) - mean(subset(dttmle, A == 0)$Y)      
    } # end for loop: istage
  } # end for loop: itrial
  
  return(output)
}



simtrials.Odat.nocalibration.nointerim <- function(dt, v){
  
  # parameters set-up
  nsim <- v$nsim
  nmax <- v$nmax
  K <- v$K
  enrollrate <- v$enrollrate
  randomA <- v$randomA
  exL <- v$exL
  exW <- v$exW
  A_to_L1 <- v$A_to_L1
  L1_to_Y <- v$L1_to_Y
  interim_times  <- (1:K)/K * nmax/enrollrate + A_to_L1 + L1_to_Y # so that we have (k/K)*nmax patients with final outcome at k-th inteirm analysis
  decision_times <- pmin(interim_times + A_to_L1 + L1_to_Y, interim_times[K])
  # at K-th interim analysis, we have everyone's final outcome observed.
  
  # output initializaiton
  output <- list(ltmle_est = matrix(NA, nrow = nsim, ncol = K),
                 unadj_est = matrix(NA, nrow = nsim, ncol = K),
                 randomA = randomA,
                 exL = exL,
                 exW = exW)
  
  dtaug <- dt # No twins! original data!
  
  
  for( itrial in 1:nsim ){
    
    if( itrial %% 10 == 0){
      cat(paste0("\nTrial number ", itrial))
    }
    
    # Sample from the twin data with replacement, to have desired sample size
    ids <- sample(nrow(dtaug), nmax, replace = TRUE)
    dtc <- dtaug[ids, ]
    
#     # Step 1: Calibrate Y with pY, to mimic Treatment Effect
#     dtc <- DGM.step1.pY.mimicTE(dtc, pY = v$pY, direction = v$dir)
#     
#     # Step 2: Calibrate Y with e0 and e1, to mimic Relative Efficiency
#     dtc <- DGM.step2.e.mimicRE(dtc, e = v$e, p0 = v$p0, p1 = v$p1)
    
    #     if (randomA == TRUE) { # Reset A as Bernoulli(0.5), i.e. assume no effect
    #       dtc$A <- rbinom(nrow(dtc), 1, 0.5)      
    #     }
    #     if (exL == TRUE) { # Set L exogeneous (complete permute)
    #       dtc$L1 <- sample(dtc$L1, nrow(dtc), replace = FALSE)
    #     }
    #     if (exW == TRUE) { # Set W exogeneous (complete permute)
    #       dtc$W1 <- sample(dtc$W1, nrow(dtc), replace = FALSE)
    #       dtc$W2 <- sample(dtc$W2, nrow(dtc), replace = FALSE)
    #     }
    
    # Assign enrollment time for each patient
    
    dtc$enrolltime <- sample(1:nmax, nmax) / enrollrate
    dtc$L1time <- dtc$enrolltime + A_to_L1
    dtc$Ytime  <- dtc$L1time + L1_to_Y
    
    for( istage in 1:K ){ # Do each interim analysis and decision analysis, compute the estimators
      
            ## Interim analysis
            
            intrtime <- interim_times[istage]
            
            C0 <- as.numeric(dtc$enrolltime <= intrtime)
            C1 <- as.numeric(dtc$L1time     <= intrtime)
            C2 <- as.numeric(dtc$Ytime      <= intrtime)
      #       
      #       dttmle <- data.frame(C0 = C0, S = dtc$S, W1 = dtc$W1, W2 = dtc$W2, A = dtc$A,
      #                            C1 = C1, L1 = dtc$L1,
      #                            C2 = C2, Y = dtc$Y)
      #       
      # suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
      #                           Anodes = "A",
      #                           Cnodes = c("C0","C1","C2"),
      #                           Lnodes = c("S","W1","W2","L1"),
      #                           Ynodes = "Y",
      #                           abar = 1,
      #                           estimate.time=FALSE))
      # suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
      #                           Anodes = "A",
      #                           Cnodes = c("C0","C1","C2"),
      #                           Lnodes = c("S","W1","W2","L1"),
      #                           Ynodes = "Y",
      #                           abar = 0,
      #                           estimate.time=FALSE))
      #       
      #       output$ltmle_est[itrial, 2*istage - 1] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      #       output$unadj_est[itrial, 2*istage - 1] <- 
      #         mean(subset(dttmle, (C2 == 1) & (A == 1))$Y) - mean(subset(dttmle, (C2 == 1) & (A == 0))$Y)
      
      ## Decision analysis
      
      idx_decision <- which(C0==1) # decision analysis should be based on these people
      
      dt_decision <- dtc[idx_decision, ]
      
      dttmle <- data.frame(S = dt_decision$S, W3 = dt_decision$W3,
                           A = dt_decision$A, Y = dt_decision$Y)
      
      suppressMessages(ltmle.fit1 <- ltmle(dttmle, 
                                           Anodes = "A",
                                           Ynodes = "Y",
                                           abar = 1,
                                           estimate.time=FALSE))
      suppressMessages(ltmle.fit0 <- ltmle(dttmle, 
                                           Anodes = "A",
                                           Ynodes = "Y",
                                           abar = 0,
                                           estimate.time=FALSE))
      
      output$ltmle_est[itrial, istage] <- ltmle.fit1$estimates[1]-ltmle.fit0$estimates[1]      
      output$unadj_est[itrial, istage] <- 
        mean(subset(dttmle, A == 1)$Y) - mean(subset(dttmle, A == 0)$Y)      
    } # end for loop: istage
  } # end for loop: itrial
  
  return(output)
}


