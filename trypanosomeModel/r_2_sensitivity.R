require("sensitivity")
library(parallel)
library(gridExtra)
library(grid)
library(plyr)
library(scales)
library(xlsx)

#************r scripts for model and parameters************
source("r_1_model.R")
source("r_1_params.R")
#********************************************************
randSnd <- read.csv("LHSParms.csv") # read in parameter values from LHS
randSnd <- randSnd[,-1]
#*****************************************************************

#******************simulate*******************************
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

start.time <- Sys.time()

trypspop <- parRapply(cl,randSnd,function(x){
  source("r_1_model.R")
  source("r_1_params.R")
  trans <- rep(seq(-19,10,1)*1000,1)
  
  sim.res <- tsetse_run_mod(tsetse_params(md.p=0.00000478
                                          ,mu.a.b=0.03
                                          ,alph=0.5
                                          ,mo.val=0#0.152
                                          ,b=x[8]
                                          ,pvit=x[1]
                                          ,pvig=x[2]
                                          ,ti=x[3]
                                          ,hw=x[4]
                                          ,phi=x[5]
                                          ,hi=x[6]
                                          ,hr=x[7]))
  dat <- output_func(sim.res)
  dat$run <- x[9]
  perc.vector <- ((dat$ivlt+dat$evlt1 + dat$evlt2 + dat$evlt3)/rowSums(dat[,1:6]))*100
  perc.host <- ((dat$ihlt)/(rowSums(dat[,7:10])))*100
  num.tsetse <- rowSums(dat[,1:6])
  output <- cbind.data.frame(run=dat$run,perc.vector,perc.host,num.tsetse,trans)
  
  return(output[10,])
})

stopCluster(cl)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

dat <- do.call(rbind.data.frame,trypspop)
write.csv(dat,"trypanosomeSensitivityResults.csv")
