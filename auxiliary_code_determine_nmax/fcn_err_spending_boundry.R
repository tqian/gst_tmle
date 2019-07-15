# Compute error spending version of the testing boundries.
# err_spending_bdry: According to Sec 4.1.1 Method 1 of Hampson & Jennison's paper, equations (12)-(15).
# The paper: Hampson & Jennison, "Group Sequential Tests for Delayed Responses." JRSSB, 2013

library(mvtnorm)


err_spending_bdry.new <- function(mvmean0, mvcov0, mvmean1, mvcov1, errs, v = NULL){
  ## This is NONBINDING! ##
  # mvmean0, mvcov0: the mean and covariance matrix of the estimators simulated under H0
  # mvmean1, mvcov1: the mean and covariance matrix of the estimators simulated under H1
  # Note that, the vector of estimators is of length 2K, with order
  #      Z_1,\tilde{Z_1},...,Z_K,\tilde{Z_K}.
  # errs: a matrix of two rows, as output of compute.errs(). 1st row: alphas; 2nd row: betas.
  #       alphas, betas: the type 1,2 errors to spend at each interim analysis
  # v: the list of simulation scenario (only need the K in v); if not specified, let K = length(alphas)
  # According to Sec 4.1.1 Method 1 of Hampson & Jennison's paper, equations (12)-(15).
  
  ##### Note: The double uniroot might not return a solution, which may
  #####       indicate too wide threshold.
  alphas <- errs[1, ]
  betas <- errs[2, ]
  if( is.null(v) ){
    K = length(alphas)
  } else {
    K <- v$K
  }  
  if ( length(mvmean0) != (2*K) ){ stop("Length of mvmean0 should be 2K.") }
  if ( length(mvmean1) != (2*K) ){ stop("Length of mvmean1 should be 2K.") }
  if ( identical(dim(mvcov0), c(2*K, 2*K)) ){ stop("Dimension of mvcov0 should be 2K*2K.") }
  if ( identical(dim(mvcov1), c(2*K, 2*K)) ){ stop("Dimension of mvcov1 should be 2K*2K.") }
  if ( length(alphas) != K) { stop("Length of alphas should be K.") }
  if ( length(betas) != K) { stop("Length of betas should be K.") }
  
  # The thresholds: u[k] (efficacy), c[k] (decision analysis)
  u <- c <- as.numeric(rep(NA, K))
  
  for(k in 1:K){
    # do binary search to get all the u and l's
    # According to Sec 4.1.1 Method 1 of Hampson & Jennison's paper, equations (12)-(15).
    
    if (k == 1){
      
      u[1] <- uniroot( function(uu){ # outer
        
        cc <- uniroot( function(ccc){ # inner
          pmvnorm(lower = c(uu, -Inf), upper = c(Inf, ccc),
                  mean = mvmean1[1:2], sigma = mvcov1[1:2, 1:2]) - betas[1] },
        interval = c(2.5, 2.6), extendInt = "yes")$root
        
        return(pmvnorm(lower = c(uu, cc), upper = c(Inf, Inf),
                       mean = mvmean0[1:2], sigma = mvcov0[1:2, 1:2]) - alphas[1]) },
        interval = c(2.5, 2.6), extendInt = "yes")$root
      
      c[1] <- uniroot( function(ccc){
        pmvnorm(lower = c(u[1], -Inf), upper = c(Inf, ccc),
                mean = mvmean1[1:2], sigma = mvcov1[1:2, 1:2]) - betas[1] },
        interval = c(-10, 10), extendInt = "yes")$root
    } else if (k < K){ ##### Not finished!
      idx <- c(seq(from = 1, by = 2, length.out = k), 2*k) # idx for computing u[k] and c[k]
      
      u[k] <- uniroot( function(uu){ # outer
        
        cc <- uniroot( function(ccc){ # inner
          pmvnorm(lower = c(rep(-Inf,k-1), uu, -Inf), upper = c(u[1:(k-1)], Inf, ccc),
                  mean = mvmean1[idx], sigma = mvcov1[idx, idx]) - betas[k] },
          interval = c(-10, 10), extendInt = "yes")$root
        
        return(pmvnorm(lower = c(rep(-Inf,k-1), uu, cc), upper = c(u[1:(k-1)], Inf, Inf),
                       mean = mvmean0[idx], sigma = mvcov0[idx, idx]) - alphas[k]) },
        interval = c(0, 2), extendInt = "yes")$root
      
      c[k] <- uniroot( function(ccc){ # inner
        pmvnorm(lower = c(rep(-Inf,k-1), u[k], -Inf), upper = c(u[1:(k-1)], Inf, ccc),
                mean = mvmean1[idx], sigma = mvcov1[idx, idx]) - betas[k] },
        interval = c(-10, 10), extendInt = "yes")$root
    } else if (k == K) {
      idxc <- c(seq(from = 1, by = 2, length.out = (K-1)), 2*K)
      
      u[K] <- Inf
      l[K] <- -Inf
      c[K] <- uniroot( # equation (15)
        function(x){
          pmvnorm(lower = c(l[1:(K-1)], x), upper = c(u[1:(K-1)], Inf),
                  mean = mvmean0[idxc], sigma = mvcov0[idxc, idxc]) - alphas[K] },
        interval = c(-10, 10), extendInt = "yes")$root
    }
  }
  
  return(list(u = u, l = l, c = c))
}

compute.errs.new <- function(covmat, v){
  # cov: an original covariance matrix (of dim 2K by 2K)
  # Note that, the vector of estimators is of length 2K, with order
  #      Z_1,\tilde{Z_1},...,Z_K,\tilde{Z_K}.
  # Output: a two row matrix, 1st row is alphas, 2nd row is betas.
  #         Errors that will be spent at each analysis.
  ##### The errors will be based on \tilde{I_k}.
  # This will be used in the err_spending_bdry() function above.
  
  alpha <- v$alpha
  beta <- v$beta
  rho <- v$rho
  f_err <- v$f_err
  g_err <- v$g_err
  K <- v$K
  
  Imax <- 1 / covmat[nrow(covmat), ncol(covmat)]
  I <- rep(NA, K)
  for(k in 1:(K-1)){
    I[k] <- 1 / covmat[2*k, 2*k] # The errors will be based on \tilde{I_k}.
  }
  I[K] <- Imax
  Iratio <- I/Imax # argument to pass into f_err and g_err
  
  fk <- sapply(Iratio, f_err)
  fk_ <- sapply(c(0, Iratio[1:(K-1)]), f_err)
  alphas <- fk - fk_ # as equation (12)
  
  gk <- sapply(Iratio, g_err)
  gk_ <- sapply(c(0, Iratio[1:(K-1)]), g_err)
  betas <- gk - gk_ # as equation (13)
  
  return(rbind(alphas, betas))
}

# For debugging, try different .rho values
# .rho <- 10
# f_err <- function (x, .alpha = 0.025) {
#   if (x < 0) {
#     f <- 0
#   }
#   else if (x < 1) {
#     f <- .alpha * x^.rho
#   }
#   else {
#     f <- .alpha
#   }
#   return(f)
# }
# g_err <- function (x, .beta = 0.2) 
# {
#   if (x < 0) {
#     g <- 0
#   }
#   else if (x < 1) {
#     g <- .beta * x^.rho
#   }
#   else {
#     g <- .beta
#   }
#   return(g)
# }



err_spending_bdry <- function(mvmean0, mvcov0, mvmean1, mvcov1, errs, v = NULL){
  # mvmean0, mvcov0: the mean and covariance matrix of the estimators simulated under H0
  # mvmean1, mvcov1: the mean and covariance matrix of the estimators simulated under H1
  # Note that, the vector of estimators is of length 2K, with order
  #      Z_1,\tilde{Z_1},...,Z_K,\tilde{Z_K}.
  # errs: a matrix of two rows, as output of compute.errs(). 1st row: alphas; 2nd row: betas.
  #       alphas, betas: the type 1,2 errors to spend at each interim analysis
  # v: the list of simulation scenario (only need the K in v); if not specified, let K = length(alphas)
  # According to Sec 4.1.1 Method 1 of Hampson & Jennison's paper, equations (12)-(15).
  
  alphas <- errs[1, ]
  betas <- errs[2, ]
  if( is.null(v) ){
    K = length(alphas)
  } else {
    K <- v$K
  }  
  if ( length(mvmean0) != (2*K) ){ stop("Length of mvmean0 should be 2K.") }
  if ( length(mvmean1) != (2*K) ){ stop("Length of mvmean1 should be 2K.") }
  if ( identical(dim(mvcov0), c(2*K, 2*K)) ){ stop("Dimension of mvcov0 should be 2K*2K.") }
  if ( identical(dim(mvcov1), c(2*K, 2*K)) ){ stop("Dimension of mvcov1 should be 2K*2K.") }
  if ( length(alphas) != K) { stop("Length of alphas should be K.") }
  if ( length(betas) != K) { stop("Length of betas should be K.") }
  
  # The thresholds: u[k] (efficacy), l[k] (futility), c[k] (decision analysis)
  u <- l <- c <- as.numeric(rep(NA, K))
  
  for(k in 1:K){
    # do binary search to get all the u and l's
    # According to Sec 4.1.1 Method 1 of Hampson & Jennison's paper, equations (12)-(15).
    
    if (k == 1){
      u[1] <- qnorm(alphas[1], mvmean0[1], sqrt(mvcov0[1,1]), lower.tail = FALSE) # equation before (12)
      l[1] <- qnorm(betas[1], mvmean1[1], sqrt(mvcov1[1,1])) # equation before (12)
      c[1] <- uniroot( # equation (14)
        function(x){
          pmvnorm(lower = c(u[1], -Inf), upper = c(Inf, x),
                  mean = mvmean0[1:2], sigma = mvcov0[1:2, 1:2]) -
            pmvnorm(lower = c(-Inf, x), upper = c(l[1], Inf),
                    mean = mvmean0[1:2], sigma = mvcov0[1:2, 1:2]) },
        interval = c(-10, 10), extendInt = "yes")$root
    } else if (k < K){
      idx <- seq(from = 1, by = 2, length.out = k) # idx for computing u[k] and l[k]
      idxc <- c(idx, 2*k) # idx for computing c[k]
      
      u[k] <- uniroot( # equation (12)
        function(x){
          pmvnorm(lower = c(l[1:(k-1)], x), upper = c(u[1:(k-1)], Inf),
                  mean = mvmean0[idx], sigma = mvcov0[idx, idx]) - alphas[k] },
        interval = c(-10, 10), extendInt = "yes")$root 
      l[k] <- uniroot( # equation (13)
        function(x){
          pmvnorm(lower = c(l[1:(k-1)], -Inf), upper = c(u[1:(k-1)], x),
                  mean = mvmean1[idx], sigma = mvcov1[idx, idx]) - betas[k] },
        interval = c(-10, 10), extendInt = "yes")$root
      c[k] <- uniroot( # equation (14)
        function(x){
          pmvnorm(lower = c(l[1:(k-1)], u[k], -Inf), upper = c(u[1:(k-1)], Inf, x),
                  mean = mvmean0[idxc], sigma = mvcov0[idxc, idxc]) - 
            pmvnorm(lower = c(l[1:(k-1)], -Inf, x), upper = c(u[1:(k-1)], l[k], Inf),
                    mean = mvmean0[idxc], sigma = mvcov0[idxc, idxc]) },
        interval = c(-10, 10), extendInt = "yes")$root
    } else if (k == K) {
      idxc <- c(seq(from = 1, by = 2, length.out = (K-1)), 2*K)
      
      u[K] <- Inf
      l[K] <- -Inf
      c[K] <- uniroot( # equation (15)
        function(x){
          pmvnorm(lower = c(l[1:(K-1)], x), upper = c(u[1:(K-1)], Inf),
                  mean = mvmean0[idxc], sigma = mvcov0[idxc, idxc]) - alphas[K] },
        interval = c(-10, 10), extendInt = "yes")$root
    }
  }
  
  return(list(u = u, l = l, c = c))
}


err_spending_bdry_nonbinding <- function(mvmean0, mvcov0, mvmean1, mvcov1, errs, v = NULL){
  # mvmean0, mvcov0: the mean and covariance matrix of the estimators simulated under H0
  # mvmean1, mvcov1: the mean and covariance matrix of the estimators simulated under H1
  # errs: a matrix of two rows, as output of compute.errs(). 1st row: alphas; 2nd row: betas.
  #       alphas, betas: the type 1,2 errors to spend at each interim analysis
  # v: the list of simulation scenario (only need the K in v); if not specified, let K = length(alphas)
  # According to Sec 4.1.1 Method 1 of Hampson & Jennison's paper, equations (12)-(15).
  
  ##### THIS IS INAPPROPRIATE!
  ## Because the c_k's will be -10, except the last c_K.
  ## This is equivalent to that, at interim analyses we make decisions immediately, don't wait for the pipeline.
  
  alphas <- errs[1, ]
  betas <- errs[2, ]
  if( is.null(v) ){
    K = length(alphas)
  } else {
    K <- v$K
  }  
  if ( length(mvmean0) != (2*K) ){ stop("Length of mvmean0 should be 2K.") }
  if ( length(mvmean1) != (2*K) ){ stop("Length of mvmean1 should be 2K.") }
  if ( identical(dim(mvcov0), c(2*K, 2*K)) ){ stop("Dimension of mvcov0 should be 2K*2K.") }
  if ( identical(dim(mvcov1), c(2*K, 2*K)) ){ stop("Dimension of mvcov1 should be 2K*2K.") }
  if ( length(alphas) != K) { stop("Length of alphas should be K.") }
  if ( length(betas) != K) { stop("Length of betas should be K.") }
  
  # The thresholds: u[k] (efficacy), l[k] (futility), c[k] (decision analysis)
  u <- l <- c <- as.numeric(rep(NA, K))
  
  for(k in 1:K){
    # do binary search to get all the u and l's
    # According to Sec 4.1.1 Method 1 of Hampson & Jennison's paper, equations (12)-(15).
    
    if (k == 1){
      u[1] <- qnorm(alphas[1], mvmean0[1], sqrt(mvcov0[1,1]), lower.tail = FALSE) # equation before (12)
      l[1] <- -Inf
      c[1] <- uniroot( # equation (14)
        function(x){
          pmvnorm(lower = c(u[1], -Inf), upper = c(Inf, x),
                  mean = mvmean0[1:2], sigma = mvcov0[1:2, 1:2]) -
            pmvnorm(lower = c(-Inf, x), upper = c(l[1], Inf),
                    mean = mvmean0[1:2], sigma = mvcov0[1:2, 1:2]) },
        interval = c(-10, 10), extendInt = "yes")$root
    } else if (k < K){
      idx <- seq(from = 1, by = 2, length.out = k) # idx for computing u[k] and l[k]
      idxc <- c(idx, 2*k) # idx for computing c[k]
      
      u[k] <- uniroot( # equation (12)
        function(x){
          pmvnorm(lower = c(l[1:(k-1)], x), upper = c(u[1:(k-1)], Inf),
                  mean = mvmean0[idx], sigma = mvcov0[idx, idx]) - alphas[k] },
        interval = c(-10, 10), extendInt = "yes")$root 
      l[k] <- -Inf
      c[k] <- uniroot( # equation (14)
        function(x){
          pmvnorm(lower = c(l[1:(k-1)], u[k], -Inf), upper = c(u[1:(k-1)], Inf, x),
                  mean = mvmean0[idxc], sigma = mvcov0[idxc, idxc]) - 
            pmvnorm(lower = c(l[1:(k-1)], -Inf, x), upper = c(u[1:(k-1)], l[k], Inf),
                    mean = mvmean0[idxc], sigma = mvcov0[idxc, idxc]) },
        interval = c(-10, 10), extendInt = "yes")$root
    } else if (k == K) {
      idxc <- c(seq(from = 1, by = 2, length.out = (K-1)), 2*K)
      
      u[K] <- Inf
      l[K] <- -Inf
      c[K] <- uniroot( # equation (15)
        function(x){
          pmvnorm(lower = c(l[1:(K-1)], x), upper = c(u[1:(k-1)], Inf),
                  mean = mvmean0[idxc], sigma = mvcov0[idxc, idxc]) - alphas[K] },
        interval = c(-10, 10), extendInt = "yes")$root
    }
  }
  
  return(list(u = u, l = l, c = c))
}

compute.errs <- function(covmat, v){
  # cov: an original covariance matrix (of dim 2K by 2K)
  # Note that, the vector of estimators is of length 2K, with order
  #      Z_1,\tilde{Z_1},...,Z_K,\tilde{Z_K}.
  # Output: a two row matrix, 1st row is alphas, 2nd row is betas.
  #         Errors that will be spent at each analysis.
  # This will be used in the err_spending_bdry() function above.
  
  alphas <- sapply((1:5)/5, v$f_err, .rho = 4) - c(0, sapply((1:4)/5, v$f_err, .rho = 4))
  betas <- sapply((1:5)/5, v$g_err, .rho = 4) - c(0, sapply((1:4)/5, v$g_err, .rho = 4))
  errs <- rbind(alphas, betas)
  
  return(errs)
}


# The following is a should-be version of the classic error spending approach
# based on information levels. However, it produces strange results.
# Specificaly, alpha_5 < alpha_4, and beta_5 < beta_4.
# I thought the reasonable thing should be to spend most errors at the final stage?
#
# compute.errs <- function(covmat, v){
#   # cov: an original covariance matrix (of dim 2K by 2K)
#   # Note that, the vector of estimators is of length 2K, with order
#   #      Z_1,\tilde{Z_1},...,Z_K,\tilde{Z_K}.
#   # Output: a two row matrix, 1st row is alphas, 2nd row is betas.
#   #         Errors that will be spent at each analysis.
#   # This will be used in the err_spending_bdry() function above.
#   
#   alpha <- v$alpha
#   beta <- v$beta
#   rho <- v$rho
#   f_err <- v$f_err
#   g_err <- v$g_err
#   K <- v$K
#   
#   Imax <- 1 / covmat[nrow(covmat), ncol(covmat)]
#   I <- rep(NA, K)
#   for(k in 1:(K-1)){
#     I[k] <- 1 / covmat[2*k-1, 2*k-1]
#   }
#   I[K] <- Imax
#   Iratio <- I/Imax # argument to pass into f_err and g_err
#   
#   fk <- sapply(Iratio, f_err)
#   fk_ <- sapply(c(0, Iratio[1:(K-1)]), f_err)
#   alphas <- fk - fk_ # as equation (12)
#   
#   gk <- sapply(Iratio, g_err)
#   gk_ <- sapply(c(0, Iratio[1:(K-1)]), g_err)
#   betas <- gk - gk_ # as equation (13)
#   
#   return(rbind(alphas, betas))
# }


# ## For debug:
# 
# v <- list(nsim = 1000,
#              n1max = 500,
#              n2max = 500,
#              K = 5,
#              enrollrate1 = 140,
#              enrollrate2 = 140,
#              randomA = FALSE, # randomA = FALSE: H1
#              exL = FALSE,
#              exW = FALSE,
#              A_to_L1 = 30 / 365,
#              L1_to_Y = (180 - 30) / 365,
#              alpha = 0.025,
#              beta = 0.2,
#              rho = 2,
#              f_err = 
#                function(x, .rho = rho, .alpha = alpha){
#                  if (x < 0){
#                    f <- 0
#                  } else if (x < 1){
#                    f <- .alpha * x^.rho
#                  } else{
#                    f <- .alpha
#                  }
#                  return(f)
#                },
#              g_err =
#                function(x, .rho = rho, .beta = beta){
#                  if (x < 0){
#                    g <- 0
#                  } else if (x < 1){
#                    g <- .beta * x^.rho
#                  } else{
#                    g <- .beta
#                  }
#                  return(g)
#                }
# )
# 
# covH1 <- compcov("est_H1", 1:50)
# covH0 <- compcov("est_H0", 1:50)
# 
# mvmean0 <- covH0$mean_ltmle_std
# mvmean1 <- covH1$mean_ltmle_std
# mvcov0 <- covH0$cov_ltmle_std
# mvcov1 <- covH1$cov_ltmle_std
# 
# errs <- compute.errs(covH0$cov_ltmle, v)
# 
# bdry_ltmle <- err_spending_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)
# 
# mvmean0 <- covH0$mean_unadj_std
# mvmean1 <- covH1$mean_unadj_std
# mvcov0 <- covH0$cov_unadj_std
# mvcov1 <- covH1$cov_unadj_std
# 
# errs <- compute.errs(covH0$cov_unadj, v)
# 
# bdry_unadj <- err_spending_bdry(mvmean0, mvcov0, mvmean1, mvcov1, errs)






# # Error spending functions, not used in this code. Keep for record.

# .alpha <- 0.025
# .beta <- 0.2
# .K = 5
# 
# .f_err <- function(x, rho = 2, alpha = .alpha){
#   if (x < 0){
#     f <- 0
#   } else if (x < 1){
#     f <- alpha * x^rho
#   } else{
#     f <- alpha
#   }
#   return(f)
# }
# 
# .g_err <- function(x, rho = 2, beta = .beta){
#   if (x < 0){
#     g <- 0
#   } else if (x < 1){
#     g <- beta * x^rho
#   } else{
#     g <- beta
#   }
#   return(g)
# }