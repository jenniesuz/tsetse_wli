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

#******************Model constants************************
area.length = 30       # length and width of total area
area.width = 4
length.inside <- 20      # length inside protected area 
length.outside <- 10     # length outside protected area

starting.adults <- 5000  # adult tsetse at day 1 / km2
starting.pupae <- 5000 # pupae at day 1 / km2

days <- 3*365*2           # days to run simulation - allow enough time for population to reach carrying capacity and tryps equilibrium
daysseq <- seq(2,days,1) # run 0.5 days
#**********************************************************
#*****************Model parameters*************************
tsetse_params <- function(   l = (1/10)/2           ## proportion of adults producing larvae per day
                             , p = (1/30)/2         ## proportion of pupae that emerge as adults  
                             , md.p = (0.00006)/2   ## 0.00006 pupal density-dependent mortality
                             , mup = (0.006)/2      ## daily probability of pupal mortality
                             , mu.a.b =(0.03)/2     ## 1/40 0.018 ## adult mortality
                             , alph = (0.5)/2       ## adult diffusion coefficient
                             , mo.val=(0.06)/2      ## adult mortality rate outside protected area
)       
   return(as.list(environment()))
#*************************************************************
#*****************summarise_output***************

output_func <- function(mod_output){
   
   #*****************subset vector data*************
   v<-mod_output[[1]]     # each vector
   p<-mod_output[[2]]
   vl <- v[,,days]  # take the last time point (end of model run)
   pl <- p[,,days]
   vlt <- vl[1,]        # take one row of the 2D world - one transect from inside to outside protected area
   plt <- pl[1,]
   trans<-seq(-19,10,1)*1000
   return(cbind.data.frame(trans,vlt,plt))
}
#**************************************************************************


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
                                    ,NUMSAMPLES=1
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
   perc <- round(temp$vlt/max(temp$vlt)*100,2)
   return(c(perc[21]))
 })
 
 randSnd$flypop <- flypopout
 
 write.csv(randSnd,"sensitAllParms200214.csv") 

 rm(list=ls())
