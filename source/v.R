v <- list(nsim = 500,
          nmax = 500,
          K = 5,
          enrollrate = 140,
          pY = 0.02955846,
          dir = "up",
          e = 0.1640625,
          p0 = 0.2222222, # = mean(subset(dt100, A == 0)$Y) = P(Y=1|A=0) 
          p1 = 0.34375, # # = mean(subset(dt100, A == 1)$Y) = P(Y=1|A=1)
          randomA = FALSE,
          exL = FALSE,
          exW = FALSE,
          A_to_L1 = 30 / 365,
          L1_to_Y = (180 - 30) / 365,
          alpha = 0.025,
          beta = 0.2,
          rho = 2,
          f_err = 
            function(x, .rho = 2, .alpha = 0.025){
              if (x < 0){
                f <- 0
              } else if (x < 1){
                f <- .alpha * x^.rho
              } else{
                f <- .alpha
              }
              return(f)
            },
          g_err =
            function(x, .rho = 2, .beta = 0.2){
              if (x < 0){
                g <- 0
              } else if (x < 1){
                g <- .beta * x^.rho
              } else{
                g <- .beta
              }
              return(g)
            }
)
