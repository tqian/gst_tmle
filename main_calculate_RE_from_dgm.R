#############################################################################
# 2017.08.15
# Tianchen Qian
#
# compute treatment effect heterogeneity and R^2_W, R^2_L for MISTIE data set
#
#############################################################################



rm(list = ls())

library(ltmle)

# Parallel setup ----------------------------------------------------------

source("source/fcn_simtrials.R")
source("source/v.R") # config of trial

dt100 <- read.csv("data/data100pts(updated).csv")


varnames <- c("scn", "Hypothesis", "rsq_w", "rsq_lgivenw", "rsq_l", "gamma",
              "rsq_w.a0", "rsq_lgivenw.a0", "rsq_l.a0",
              "rsq_w.a1", "rsq_lgivenw.a1", "rsq_l.a1")
output <- as.data.frame(matrix(NA, nrow = 8, ncol = length(varnames)))
names(output) <- varnames

irow <- 1

for (exW in c(FALSE, TRUE)) {
  for (exL in c(FALSE, TRUE)) {
    for (H0 in c(FALSE, TRUE)) {
      v$nsim <- 1
      v$nmax <- 1000000
      v$enrollrate <- 140
      v$exL <- exL
      v$exW <- exW
      
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
      
      
      ## Generate a single set from the generative model.
      ## The following code is copied from "fcn_simtrials".
      
      set.seed(20190725)
      
      dtaug <- DGM.step0.twins(dt100)
      ids <- sample(nrow(dtaug), v$nmax, replace = TRUE)
      dtc <- dtaug[ids, ]
      
      # Step 1: Calibrate Y with pY, to mimic Treatment Effect
      dtc <- DGM.step1.pY.mimicTE(dtc, pY = v$pY, direction = v$dir)
      
      # Step 2: Calibrate Y with e, to mimic Relative Efficiency
      dtc <- DGM.step2.e.mimicRE(dtc, e = v$e, p0 = v$p0, p1 = v$p1)
      
      if (H0) { # Reset A as Bernoulli(0.5), i.e. assume no effect
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
      
      
      
      dt_a1 <- subset(dtc, A == 1)
      dt_a0 <- subset(dtc, A == 0)
      dt_full <- dtc
      
      # W_vector <- "S + W1 + W2 + W3"
      W_vector <- "S + W3"
      
      fit_w.a1 <- glm(as.formula(paste("Y ~ ", W_vector)), data = dt_a1, family = binomial)
      fit_wl.a1 <- glm(as.formula(paste("Y ~ L1 + ", W_vector)), data = dt_a1, family = binomial)
      fit_l.a1 <- glm(Y ~ L1, data = dt_a1, family = binomial)
      
      fit_w.a0 <- glm(as.formula(paste("Y ~ ", W_vector)), data = dt_a0, family = binomial)
      fit_wl.a0 <- glm(as.formula(paste("Y ~ L1 + ", W_vector)), data = dt_a0, family = binomial)
      fit_l.a0 <- glm(Y ~ L1, data = dt_a0, family = binomial)
      
      var_fit_w.a1 <- var(predict(fit_w.a1, newdata = dt_a1, type = "response"))
      var_fit_lgivenw.a1 <- var(predict(fit_wl.a1, newdata = dt_a1, type = "response")
                                - predict(fit_w.a1, newdata = dt_a1, type = "response"))
      var_fit_l.a1 <- var(predict(fit_l.a1, newdata = dt_a1, type = "response"))
      
      var_fit_w.a0 <- var(predict(fit_w.a0, newdata = dt_a0, type = "response"))
      var_fit_lgivenw.a0 <- var(predict(fit_wl.a0, newdata = dt_a0, type = "response")
                                - predict(fit_w.a0, newdata = dt_a0, type = "response"))
      var_fit_l.a0 <- var(predict(fit_l.a0, newdata = dt_a0, type = "response"))
      var_y.a1 <- var(dt_a1$Y)
      var_y.a0 <- var(dt_a0$Y)
      
      rsq_w <- (var_fit_w.a1 + var_fit_w.a0) / (var_y.a1 + var_y.a0)
      rsq_lgivenw <- (var_fit_lgivenw.a1 + var_fit_lgivenw.a0) / (var_y.a1 + var_y.a0)
      rsq_l <- (var_fit_l.a1 + var_fit_l.a0) / (var_y.a1 + var_y.a0)
      gamma <- var(predict(fit_w.a1, newdata = dt_full, type = "response")
                   - predict(fit_w.a0, newdata = dt_full, type = "response")) /
        (var_y.a1 + var_y.a0)
      
      rsq_w.a0 <- var_fit_w.a0 / var_y.a0
      rsq_lgivenw.a0 <- var_fit_lgivenw.a0 / var_y.a0
      rsq_l.a0 <- var_fit_l.a0 / var_y.a0
      
      rsq_w.a1 <- var_fit_w.a1 / var_y.a1
      rsq_lgivenw.a1 <- var_fit_lgivenw.a1 / var_y.a1
      rsq_l.a1 <- var_fit_l.a1 / var_y.a1
      
      output$scn[irow] <- scn
      output$Hypothesis[irow] <- ifelse(H0, "H0", "H1")
      output[irow, c("rsq_w", "rsq_lgivenw", "rsq_l", "gamma")] <- c(rsq_w, rsq_lgivenw, rsq_l, gamma)
      output[irow, c("rsq_w.a0", "rsq_lgivenw.a0", "rsq_l.a0")] <- c(rsq_w.a0, rsq_lgivenw.a0, rsq_l.a0)
      output[irow, c("rsq_w.a1", "rsq_lgivenw.a1", "rsq_l.a1")] <- c(rsq_w.a1, rsq_lgivenw.a1, rsq_l.a1)
      
      irow <- irow + 1
    }
  }
}

options(digits = 8)
output
#        scn Hypothesis         rsq_w   rsq_lgivenw         rsq_l         gamma      rsq_w.a0 rsq_lgivenw.a0      rsq_l.a0      rsq_w.a1 rsq_lgivenw.a1      rsq_l.a1
# 1  prognWL         H1 3.6971712e-01 6.2576522e-02 2.9386239e-01 1.7515697e-02 4.5572165e-01  4.6995687e-02 3.2567592e-01 3.0256217e-01  7.4742509e-02 2.6902143e-01
# 2  prognWL         H0 3.4845194e-01 7.4992050e-02 3.0312715e-01 4.7891716e-07 3.4840124e-01  7.5146178e-02 3.0309111e-01 3.4850266e-01  7.4837861e-02 3.0316320e-01
# 3   prognW         H1 3.6971712e-01 5.5017571e-07 6.2391923e-07 1.7515697e-02 4.5572165e-01  4.7419035e-07 3.0348398e-09 3.0256217e-01  6.0950737e-07 1.1087245e-06
# 4   prognW         H0 3.4845194e-01 3.1848463e-06 4.4998966e-06 4.7891716e-07 3.4840124e-01  2.3632387e-06 6.0152617e-06 3.4850266e-01  4.0067804e-06 2.9839293e-06
# 5   prognL         H1 2.2014259e-07 2.9386264e-01 2.9386239e-01 1.8484712e-07 4.8991518e-07  3.2567592e-01 3.2567592e-01 9.4959954e-09  2.6902187e-01 2.6902143e-01
# 6   prognL         H0 1.5365859e-05 3.0312831e-01 3.0312715e-01 2.0026419e-05 1.9885576e-05  3.0309003e-01 3.0309111e-01 1.0844346e-05  3.0316661e-01 3.0316320e-01
# 7 prognnon         H1 5.9004334e-06 6.2379166e-07 6.2391923e-07 3.7053089e-06 1.1038899e-05  3.3316182e-09 3.0348398e-09 1.8881639e-06  1.1082656e-06 1.1087245e-06
# 8 prognnon         H0 1.5941369e-06 4.5012609e-06 4.4998966e-06 2.9190237e-06 1.4726559e-06  6.0217419e-06 6.0152617e-06 1.7156662e-06  2.9801757e-06 2.9839293e-06

options(digits = 3)
output
#        scn Hypothesis    rsq_w rsq_lgivenw    rsq_l    gamma rsq_w.a0 rsq_lgivenw.a0 rsq_l.a0 rsq_w.a1 rsq_lgivenw.a1 rsq_l.a1
# 1  prognWL         H1 3.70e-01    6.26e-02 2.94e-01 1.75e-02 4.56e-01       4.70e-02 3.26e-01 3.03e-01       7.47e-02 2.69e-01
# 2  prognWL         H0 3.48e-01    7.50e-02 3.03e-01 4.79e-07 3.48e-01       7.51e-02 3.03e-01 3.49e-01       7.48e-02 3.03e-01
# 3   prognW         H1 3.70e-01    5.50e-07 6.24e-07 1.75e-02 4.56e-01       4.74e-07 3.03e-09 3.03e-01       6.10e-07 1.11e-06
# 4   prognW         H0 3.48e-01    3.18e-06 4.50e-06 4.79e-07 3.48e-01       2.36e-06 6.02e-06 3.49e-01       4.01e-06 2.98e-06
# 5   prognL         H1 2.20e-07    2.94e-01 2.94e-01 1.85e-07 4.90e-07       3.26e-01 3.26e-01 9.50e-09       2.69e-01 2.69e-01
# 6   prognL         H0 1.54e-05    3.03e-01 3.03e-01 2.00e-05 1.99e-05       3.03e-01 3.03e-01 1.08e-05       3.03e-01 3.03e-01
# 7 prognnon         H1 5.90e-06    6.24e-07 6.24e-07 3.71e-06 1.10e-05       3.33e-09 3.03e-09 1.89e-06       1.11e-06 1.11e-06
# 8 prognnon         H0 1.59e-06    4.50e-06 4.50e-06 2.92e-06 1.47e-06       6.02e-06 6.02e-06 1.72e-06       2.98e-06 2.98e-06

saveRDS(output, "rsq_dgm.RDS")

Rsq <- output

##### compute asymptotic relative efficiency for estimating ATE and mean of Y in a single arm #####

are_ate <- function(rsq_w, rsq_lgivenw, gamma, py, pl) {
  (1 + (py/2) * gamma - rsq_w - (1 - py / pl) * rsq_lgivenw)^(-1)
}

are_single_arm <- function(rsq_w, rsq_lgivenw, py, pl) {
  (1 - (1 - py / 2) * rsq_w - (1 - py / pl) * rsq_lgivenw)^(-1)
}

# The following sample size numbers are from Table D.1 in the paper.
# For nmax = 480
interim_yobs <- c(96, 192, 288, 384)
interim_lobs <- c(57, 57, 57, 57)
interim_wobs <- c(12, 12, 12, 12)

interim_py_n480 <- interim_yobs / (interim_yobs + interim_lobs + interim_wobs)
interim_pl_n480 <- (interim_yobs + interim_lobs) / (interim_yobs + interim_lobs + interim_wobs)

decision_py_n480 <- rep(1, 5)
decision_pl_n480 <- rep(1, 5)

# For nmax = 300
interim_yobs <- c(60, 120, 180, 240)
interim_lobs <- c(57, 57, 57, 57)
interim_wobs <- c(12, 12, 12, 3)

interim_py_n300 <- interim_yobs / (interim_yobs + interim_lobs + interim_wobs)
interim_pl_n300 <- (interim_yobs + interim_lobs) / (interim_yobs + interim_lobs + interim_wobs)

decision_py_n300 <- rep(1, 5)
decision_pl_n300 <- rep(1, 5)



##### Calculate ARE from theory in Table 3 (estimating ATE) #####

## Under H0

interim_analysis_ARE <- matrix(NA, nrow = 4, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H0")
interim_analysis_ARE[, 1] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H0")
interim_analysis_ARE[, 2] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H0")
interim_analysis_ARE[, 3] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n480, interim_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H0")
interim_analysis_ARE[, 4] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n480, interim_pl_n480)

# > interim_analysis_ARE
#      [,1] [,2] [,3] [,4]
# [1,] 1.63 1.53 1.13    1
# [2,] 1.59 1.53 1.07    1
# [3,] 1.58 1.53 1.05    1
# [4,] 1.57 1.53 1.04    1

decision_analysis_ARE <- matrix(NA, nrow = 5, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H0")
decision_analysis_ARE[, 1] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H0")
decision_analysis_ARE[, 2] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H0")
decision_analysis_ARE[, 3] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n480, decision_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H0")
decision_analysis_ARE[, 4] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n480, decision_pl_n480)

# > decision_analysis_ARE
#      [,1] [,2] [,3] [,4]
# [1,] 1.53 1.53    1    1
# [2,] 1.53 1.53    1    1
# [3,] 1.53 1.53    1    1
# [4,] 1.53 1.53    1    1
# [5,] 1.53 1.53    1    1


## Under H1

interim_analysis_ARE <- matrix(NA, nrow = 4, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H1")
interim_analysis_ARE[, 1] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H1")
interim_analysis_ARE[, 2] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H1")
interim_analysis_ARE[, 3] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n480, interim_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H1")
interim_analysis_ARE[, 4] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], interim_py_n480, interim_pl_n480)

# > interim_analysis_ARE
#      [,1] [,2] [,3] [,4]
# [1,] 1.66 1.58 1.12    1
# [2,] 1.62 1.57 1.07    1
# [3,] 1.61 1.57 1.05    1
# [4,] 1.60 1.57 1.04    1


decision_analysis_ARE <- matrix(NA, nrow = 5, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H1")
decision_analysis_ARE[, 1] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H1")
decision_analysis_ARE[, 2] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H1")
decision_analysis_ARE[, 3] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n480, decision_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H1")
decision_analysis_ARE[, 4] <- are_ate(Rsq$rsq_w[rowid_in_Rsq], Rsq$rsq_lgivenw[rowid_in_Rsq], Rsq$gamma[rowid_in_Rsq], decision_py_n480, decision_pl_n480)

# > decision_analysis_ARE
#      [,1] [,2] [,3] [,4]
# [1,] 1.56 1.56    1    1
# [2,] 1.56 1.56    1    1
# [3,] 1.56 1.56    1    1
# [4,] 1.56 1.56    1    1
# [5,] 1.56 1.56    1    1





##### Calculate ARE from theory in Table G.2 (estimating mean Y in one arm under H0) #####

## under H0, estimating E(Y|A=0)

interim_analysis.a0 <- matrix(NA, nrow = 4, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H0")
interim_analysis.a0[, 1] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H0")
interim_analysis.a0[, 2] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H0")
interim_analysis.a0[, 3] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n480, interim_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H0")
interim_analysis.a0[, 4] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n480, interim_pl_n480)

# > interim_analysis.a0
#      [,1] [,2] [,3] [,4]
# [1,] 1.44 1.36 1.13    1
# [2,] 1.36 1.31 1.07    1
# [3,] 1.32 1.29 1.05    1
# [4,] 1.29 1.26 1.04    1

decision_analysis.a0 <- matrix(NA, nrow = 5, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H0")
decision_analysis.a0[, 1] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H0")
decision_analysis.a0[, 2] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H0")
decision_analysis.a0[, 3] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n480, decision_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H0")
decision_analysis.a0[, 4] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n480, decision_pl_n480)

# > decision_analysis.a0
#      [,1] [,2] [,3] [,4]
# [1,] 1.21 1.21    1    1
# [2,] 1.21 1.21    1    1
# [3,] 1.21 1.21    1    1
# [4,] 1.21 1.21    1    1
# [5,] 1.21 1.21    1    1


## under H0, estimating E(Y|A=1)

interim_analysis.a1 <- matrix(NA, nrow = 4, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H0")
interim_analysis.a1[, 1] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H0")
interim_analysis.a1[, 2] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H0")
interim_analysis.a1[, 3] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n480, interim_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H0")
interim_analysis.a1[, 4] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n480, interim_pl_n480)

# > interim_analysis.a1
#      [,1] [,2] [,3] [,4]
# [1,] 1.44 1.37 1.13    1
# [2,] 1.35 1.31 1.07    1
# [3,] 1.32 1.29 1.05    1
# [4,] 1.29 1.26 1.04    1

decision_analysis.a1 <- matrix(NA, nrow = 5, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H0")
decision_analysis.a1[, 1] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H0")
decision_analysis.a1[, 2] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H0")
decision_analysis.a1[, 3] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n480, decision_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H0")
decision_analysis.a1[, 4] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n480, decision_pl_n480)

# > decision_analysis.a1
#      [,1] [,2] [,3] [,4]
# [1,] 1.21 1.21    1    1
# [2,] 1.21 1.21    1    1
# [3,] 1.21 1.21    1    1
# [4,] 1.21 1.21    1    1
# [5,] 1.21 1.21    1    1




##### Calculate ARE from theory in Table G.3 (estimating mean Y in one arm under H1) #####

## under H1, estimating E(Y|A=0)

interim_analysis.a0 <- matrix(NA, nrow = 4, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H1")
interim_analysis.a0[, 1] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H1")
interim_analysis.a0[, 2] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H1")
interim_analysis.a0[, 3] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n480, interim_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H1")
interim_analysis.a0[, 4] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], interim_py_n480, interim_pl_n480)

# > interim_analysis.a0
#      [,1] [,2] [,3] [,4]
# [1,] 1.59 1.54 1.14    1
# [2,] 1.48 1.45 1.08    1
# [3,] 1.43 1.41 1.06    1
# [4,] 1.39 1.38 1.04    1

decision_analysis.a0 <- matrix(NA, nrow = 5, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H1")
decision_analysis.a0[, 1] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H1")
decision_analysis.a0[, 2] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H1")
decision_analysis.a0[, 3] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n480, decision_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H1")
decision_analysis.a0[, 4] <- are_single_arm(Rsq$rsq_w.a0[rowid_in_Rsq], Rsq$rsq_lgivenw.a0[rowid_in_Rsq], decision_py_n480, decision_pl_n480)

# > decision_analysis.a0
#      [,1] [,2] [,3] [,4]
# [1,]  1.3  1.3    1    1
# [2,]  1.3  1.3    1    1
# [3,]  1.3  1.3    1    1
# [4,]  1.3  1.3    1    1
# [5,]  1.3  1.3    1    1


## under H1, estimating E(Y|A=1)

interim_analysis.a1 <- matrix(NA, nrow = 4, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H1")
interim_analysis.a1[, 1] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H1")
interim_analysis.a1[, 2] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n300, interim_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H1")
interim_analysis.a1[, 3] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n480, interim_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H1")
interim_analysis.a1[, 4] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], interim_py_n480, interim_pl_n480)

# > interim_analysis.a1
#      [,1] [,2] [,3] [,4]
# [1,] 1.37 1.30 1.11    1
# [2,] 1.30 1.26 1.07    1
# [3,] 1.27 1.24 1.05    1
# [4,] 1.24 1.22 1.04    1

decision_analysis.a1 <- matrix(NA, nrow = 5, ncol = 4)
rowid_in_Rsq <- which(Rsq$scn == "prognWL" & Rsq$Hypothesis == "H1")
decision_analysis.a1[, 1] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognW" & Rsq$Hypothesis == "H1")
decision_analysis.a1[, 2] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n300, decision_pl_n300)
rowid_in_Rsq <- which(Rsq$scn == "prognL" & Rsq$Hypothesis == "H1")
decision_analysis.a1[, 3] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n480, decision_pl_n480)
rowid_in_Rsq <- which(Rsq$scn == "prognnon" & Rsq$Hypothesis == "H1")
decision_analysis.a1[, 4] <- are_single_arm(Rsq$rsq_w.a1[rowid_in_Rsq], Rsq$rsq_lgivenw.a1[rowid_in_Rsq], decision_py_n480, decision_pl_n480)

# > decision_analysis.a1
#      [,1] [,2] [,3] [,4]
# [1,] 1.18 1.18    1    1
# [2,] 1.18 1.18    1    1
# [3,] 1.18 1.18    1    1
# [4,] 1.18 1.18    1    1
# [5,] 1.18 1.18    1    1