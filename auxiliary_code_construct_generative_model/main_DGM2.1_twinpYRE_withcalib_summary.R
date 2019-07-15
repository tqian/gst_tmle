rm(list = ls())

library(ltmle)

Nparal <- 50


resultnames <- paste0("Result_RE/", "result_RE_0.1422", "_", 1:Nparal, ".rda") # result filenames from parallel jobs
  
mat_ltmle <- matrix(NA, nrow = 1, ncol = 5)
mat_unadj <- matrix(NA, nrow = 1, ncol = 5)  
for(ifile in 1:Nparal){
  load(resultnames[ifile])
  mat_ltmle <- rbind(mat_ltmle, test$ltmle_est)
  mat_unadj <- rbind(mat_unadj, test$unadj_est)
}
mat_ltmle <- mat_ltmle[-1, ]
mat_unadj <- mat_unadj[-1, ]

apply(mat_ltmle, 2, mean)
# 0.1211771 0.1214903 0.1214574 0.1214358 0.1214755
apply(mat_unadj, 2, mean)
# 0.1215537 0.1216033 0.1215617 0.1215140 0.1215288

var_ltmle <- apply(mat_ltmle, 2, var)
var_unadj <- apply(mat_unadj, 2, var)

RE <- var_unadj / var_ltmle # relative efficiency vector (stage 1 through stage 5)
# 1.578433 1.581305 1.574506 1.568643 1.554664
RE
