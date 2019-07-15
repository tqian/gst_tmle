# Gather all the parallel computing results from simtrials, then
# compute the covariance matrix and correlation matrix, for ltmle estimators and unadj estimators, respectively.

# Output: a list of: 
# 1) 4 matrices of all the estimators: ltmle & unadj, original & standardized (Z-statistic);
# 2) 4 covariance matrices of the above 4 matrices of estimators;
# 3) 4 mean vectors of the above 4 matices of estimators

compcov <- function(parallel_filenames){
  
  # initialize overall matrix to collect all parallel simulated esitmators
  mat_ltmle <- matrix(0, nrow = 0, ncol = 10)
  mat_unadj <- matrix(0, nrow = 0, ncol = 10)
  
  # rbind each resulting matrix to the overall matrix
  for (filename in parallel_filenames){
    result <- readRDS(filename)
    mat_ltmle <- rbind(mat_ltmle, result$ltmle_est)
    mat_unadj <- rbind(mat_unadj, result$unadj_est)
  }
  
  # Compute mean and covariance of ltmle and unadj
  mean_ltmle <- apply(mat_ltmle, 2, mean)
  cov_ltmle <- cov(mat_ltmle)
  mean_unadj <- apply(mat_unadj, 2, mean)
  cov_unadj <- cov(mat_unadj)
  
  # Standardize ltmle and unadj (divide by standard error, make it Z-statisitc)
  mat_ltmle_std <- mat_ltmle %*% diag(diag(cov_ltmle) ^ (-1/2))
  mat_unadj_std <- mat_unadj %*% diag(diag(cov_unadj) ^ (-1/2))
  
  # Compute mean and covariance of the standardized ones
  mean_ltmle_std <- apply(mat_ltmle_std, 2, mean)
  cov_ltmle_std <- cov(mat_ltmle_std)
  mean_unadj_std <- apply(mat_unadj_std, 2, mean)
  cov_unadj_std <- cov(mat_unadj_std)
    
  # gather the cov/corr results
  result <- list(mat_ltmle = mat_ltmle, mat_unadj = mat_unadj,
                 mean_ltmle = mean_ltmle, mean_unadj = mean_unadj,
                 cov_ltmle = cov_ltmle, cov_unadj = cov_unadj,
                 mat_ltmle_std = mat_ltmle_std, mat_unadj_std = mat_unadj_std,
                 mean_ltmle_std = mean_ltmle_std, mean_unadj_std = mean_unadj_std,
                 cov_ltmle_std = cov_ltmle_std, cov_unadj_std = cov_unadj_std)
  return(result)
}

# # Example:
# a <- compcov("est_H1", 1:50)
# lapply(a, summary)

# ## For debug
# generic_filename <- "est_H1"
# idx <- 1:3
# list_name = "test"
# directory = "Result/"


# ## For debug:
# 
# tst <- compcov("est_H1", 1:2)
# 
# # manually compute it to see if things agree
# load("Result/est_H1_1.rda")
# ltmle.mat <- test$ltmle_est
# unadj.mat <- test$unadj_est
# load("Result/est_H1_2.rda")
# ltmle.mat <- rbind(ltmle.mat, test$ltmle_est)
# unadj.mat <- rbind(unadj.mat, test$unadj_est)
# 
# tst.m <- list(apply(ltmle.mat, 2, mean), cov(ltmle.mat), cor(ltmle.mat),
#               apply(unadj.mat, 2, mean), cov(unadj.mat), cor(unadj.mat))
# 
# max(tst[[1]]-tst.m[[1]])
# max(tst[[2]]-tst.m[[2]])
# max(tst[[3]]-tst.m[[3]])
# max(tst[[4]]-tst.m[[4]])
# max(tst[[5]]-tst.m[[5]])
# max(tst[[6]]-tst.m[[6]])
# 
# # Good! All differences are small, so tst and tst.m are equal.