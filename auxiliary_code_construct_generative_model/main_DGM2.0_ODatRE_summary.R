rm(list = ls())

library(ltmle)

setwd("MISTIE/code")

Npara <- 100

filename <- "ODatRE_nmax500"

resultnames <- paste0("Result/", filename, "_", 1:Npara, ".rda") # result filenames from parallel jobs

mat_ltmle <- matrix(NA, nrow = 1, ncol = 5)
mat_unadj <- matrix(NA, nrow = 1, ncol = 5)  
for(ifile in 1:Npara){
  load(resultnames[ifile])
  mat_ltmle <- rbind(mat_ltmle, test$ltmle_est)
  mat_unadj <- rbind(mat_unadj, test$unadj_est)
}
mat_ltmle <- mat_ltmle[-1, ]
mat_unadj <- mat_unadj[-1, ]

var_ltmle <- apply(mat_ltmle, 2, var)
var_unadj <- apply(mat_unadj, 2, var)

RE <- var_unadj / var_ltmle # relative efficiency vector (stage 1 through stage 5)

write.csv(RE, file = paste0("log/", filename, ".csv"), row.names = FALSE)