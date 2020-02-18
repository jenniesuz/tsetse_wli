library(lhs)
library(spartan)
library(lattice)
library(MASS)
library(sensitivity)
library(parallel)
library(gridExtra)
library(grid)
library(plyr)
library(scales)


#************Set ranges for parameter values**************
# trypanosome parameters
ptrans.ht <- c(0,0.08)
ptrans.hg <- c(0,0.005)
vinc <- c(1/30,1/15)
wane <- c(1/100,1)
ptrans.vh <- c(0.2,0.8)
hinc <- c(1/15,1/5)
hrec <- c(1/120,1/10)
feed <- c(1/3,1/2)

params <- c("ptrans.ht"
            ,"ptrans.hg"
            ,"vinc"
            ,"wane"
            ,"ptrans.vh"
            ,"hinc"
            ,"hrec"
            ,"feed"
)
#**********generate combinations of parameter values***********************

 randSnd <- data.frame(lhc_generate_lhc_sample(getwd()
                                     ,PARAMETERS=params
                                     ,NUMSAMPLES=1000
                                     ,PMIN=c(ptrans.ht[1]
                                             ,ptrans.hg[1]
                                             ,vinc[1]
                                             ,wane[1]
                                             ,ptrans.vh[1]
                                             ,hinc[1]
                                             ,hrec[1]
                                             ,feed[1]
                                             )
                                     ,PMAX=c(ptrans.ht[2]
                                             ,ptrans.hg[2]
                                             ,vinc[2]
                                             ,wane[2]
                                             ,ptrans.vh[2]
                                             ,hinc[2]
                                             ,hrec[2]
                                             ,feed[2]
                                             )
                                     ,ALGORITHM="normal"))
 
 
 randSnd$run <- seq(1,length(randSnd[,1]),1)
write.csv(randSnd,"LHSParms.csv")


#***************************************************************
#*****************************************************************
 