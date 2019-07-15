## functions: trialresults, run1trial


run1trial <- function(zstats, bdry){
  ## zstats: a vector of length 2K, with order Z_1, \tilde{Z_1}, ..., Z_K, \tilde{Z_K}
  ## bdry: a list of 3 vectors: u_k, l_k, c_k (output of err_spending_bdry)
  
  ## This function gives the output of at which interim analysis the trial stops,
  ## and the reason for the trial stops.
  
  output <- rep(NA, 3)
  # output:
  # 1) at which interim analysis the trial stops
  # 2) stops for efficacy or futility: 1 for efficacy, 0 for futility
  # 3) decision analysis result: 1 for significant result (reject H0), 0 for insignificant result.
  
  K <- length(zstats)/2
  z <- zstats[seq(1, 2*K-1, 2)] # Z_k's
  zz <- zstats[seq(2, 2*K, 2)] # \tilde{Z_k}'s
  
  u <- bdry$u
  l <- bdry$l
  c <- bdry$c
  
  ## Compute the stopping time. According to equation (1) of Hampson Jennison's paper, page 7.
  k <- 1 # index for interim analysis
  go <- 1 # continue to next interim analysis?
  while(go){
    if ( k <= (K-1) ) { # at interim analysis k <= K-1
      if(z[k] <= l[k]){ # early stop for futility
        output[1] <- k
        output[2] <- 0
        output[3] <- ifelse(zz[k] >= c[k], 1, 0) # decision analysis
        go <- 0 # stop the trial
      } else if (z[k] >= u[k]){ # early stop for efficacy
        output[1] <- k
        output[2] <- 1
        output[3] <- ifelse(zz[k] >= c[k], 1, 0) # decision analysis
        go <- 0 # stop the trial
      } else { # continue to next interim analysis
        k <- k+1
      }
    } else { # k = K
      output[1] <- K
      output[2] <- 0 # this doesn't matter when at interim analysis K
      output[3] <- ifelse(zz[K] >= c[K], 1, 0)
      go <- 0
    }
  } # end while
  
  return(output)
}

trialresults <- function(mat_zstats, bdry){
  ## mat_zstats: a nsim by 2K matrix, each row is the zstats within a single trial.
  ## Use output of compcov$mat_ltmle_std or compcov$mat_unadj_std.
  ## bdry: a list of 3 vectors: u_k, l_k, c_k (output of err_spending_bdry).
  results <- data.frame(t(apply(mat_zstats, 1, run1trial, bdry = bdry)))
  colnames(results) <- c("whichstop", "interimresult", "decisionresult")
  return(results)
}
