library(lhs)
library(spartan)
require(lattice)
require(MASS)
require(sensitivity)
library(parallel)
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(scales)
library(RColorBrewer)

 
source("r_1_model.R")
source("r_1_params.R")

# # Set ranges for parameter values 
 larvi <- c(1/12,1/8)/2
 emerg <- c(1/50,1/30)/2
 pden <-  c(0.00001,0.00009)/2
 pmort <- c(0.0025,0.01)/2
 amort <- c(0.01,0.03)/2
 difcoeff <- c(0.2,1)/2 # 1 / day
 mo.out <- c(0.03,0.5)/2
#
 params <- c("larvi"
             ,"emerg"
             ,"pden"
             ,"pmort"
             ,"amort"
             ,"difcoeff"
             ,"mo.out"  )
 
# # select random sets of parameter values within parameter value ranges
 randSnd <- data.frame(lhc_generate_lhc_sample(getwd()
                                    ,PARAMETERS=params
                                    ,NUMSAMPLES=1000
                                    ,PMIN=c(larvi[1]
                                            ,emerg[1]
                                            ,pden[1]
                                            ,pmort[1]
                                            ,amort[1]
                                            ,difcoeff[1]
                                            ,mo.out[1])
                                    ,PMAX=c(larvi[2]
                                            ,emerg[2]
                                            ,pden[2]
                                            ,pmort[2]
                                            ,amort[2]
                                            ,difcoeff[2]
                                            ,mo.out[2])
                                    ,ALGORITHM="normal"))
 
 
 # ********run model for each parameter set******
 no_cores <- detectCores() - 1
 cl <- makeCluster(no_cores)
 
 
 flypop <- parRapply(cl,randSnd,function(x){
   source("r_1_model.R")
   source("r_1_params.R")
   
   sim.res <- tsetse_run_mod(tsetse_params(l=x[1]
                                           ,p=x[2]
                                           ,md.p=x[3]
                                           ,mup=x[4]
                                           ,mu.a.b=x[5]
                                           ,alph=x[6]
                                           ,mo.val=x[7]))
   
   
   
   output <- output_func(sim.res)
   return(output)
 })
 
 stopCluster(cl)
 
 
 flypopout <- sapply(1:length(flypop),function(x){
   temp <- flypop[[x]]
   temp <- temp[,-1]
   tots <- rowSums(temp)
   perc <- round(tots[21]/tots[1]*100,2)
   return(perc)
 })
 
 randSnd$flypop <- flypopout
 
 write.csv(randSnd,"tsetseSensitivityResults.csv") 

 rm(list=ls())
