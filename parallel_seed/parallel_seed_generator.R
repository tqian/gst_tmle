rm(list = ls())
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

Npara <- 1000

parallel_seeds <- matrix(NA, nrow = Npara, ncol = 7)

s <- .Random.seed
for (i in 1:Npara) {
  parallel_seeds[i,] <- s
  s <- nextRNGStream(s)
  # send s to worker i as .Random.seed
}

write.csv(parallel_seeds, "parallel_seeds.csv", row.names = F)