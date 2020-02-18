library(lhs)
library(spartan)
require("lattice")
require("MASS")
require("sensitivity")
library(parallel)
library(gridExtra)
library(grid)
library(plyr)
library(scales)

source("r_1_model.R")
source("r_1_params.R")

# Tb parameter values
tbrp <- read.csv("tbSens.csv")
 
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)


start.time <- Sys.time()

trypspop <- parRapply(cl,tbrp,function(x){
  source("r_1_model.R")
  source("r_1_params.R")
  trans <- seq(-19,10,1)*1000

   sim.res <- tsetse_run_mod(tsetse_params(md.p=0.00000478
                                           ,mu.a.b=0.03
                                           ,alph=0.5
                                           ,mo.val=0.152
                                           ,b=x[10]
                                           ,pvit=x[3]
                                           ,pvig=x[4]
                                           ,ti=x[5]
                                           ,hw=x[6]
                                           ,phi=x[7]
                                           ,hi=x[8]
                                           ,hr=x[9]))
   dat <- output_func(sim.res)
   dat$run <- x[11]
   perc.vector <- ((dat$ivlt+dat$evlt1 + dat$evlt2 + dat$evlt3)/rowSums(dat[,1:6]))*100
   perc.host <- ((dat$ihlt)/(rowSums(dat[,7:10])))*100
   num.tsetse <- rowSums(dat[,1:6])
   output <- cbind.data.frame(run=dat$run,perc.vector,perc.host,num.tsetse,trans)
   
   return(output[10:25,])
 })
 
 stopCluster(cl)
 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

dat <- do.call(rbind.data.frame,trypspop)
write.csv(dat,"sensTbOutput.csv")

# Tco parameter values

tcp <- read.csv("tcSens.csv")

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)


start.time <- Sys.time()

trypspop <- parRapply(cl,tcp,function(x){
  source("r_1_model.R")
  source("r_1_params.R")
  trans <- seq(-19,10,1)*1000
  
  sim.res <- tsetse_run_mod(tsetse_params(md.p=0.00000478
                                          ,mu.a.b=0.03
                                          ,alph=0.5
                                          ,mo.val=0.152
                                          ,b=x[10]
                                          ,pvit=x[3]
                                          ,pvig=x[4]
                                          ,ti=x[5]
                                          ,hw=x[6]
                                          ,phi=x[7]
                                          ,hi=x[8]
                                          ,hr=x[9]))
  dat <- output_func(sim.res)
  dat$run <- x[11]
  perc.vector <- ((dat$ivlt+dat$evlt1 + dat$evlt2 + dat$evlt3)/rowSums(dat[,1:6]))*100
  perc.host <- ((dat$ihlt)/(rowSums(dat[,7:10])))*100
  num.tsetse <- rowSums(dat[,1:6])
  output <- cbind.data.frame(run=dat$run,perc.vector,perc.host,num.tsetse,trans)
  
  return(output[10:25,])
})

stopCluster(cl)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

dat <- do.call(rbind.data.frame,trypspop)
write.csv(dat,"sensTcOutput.csv")
