library(minpack.lm)
library(ggplot2)
library(gridExtra)
library(beepr)
library(plyr)

source("r_1_model.R")

#******************Model constants************************
area.length = 60         # length and width of total area each cell of side 500 m
area.width = 4
length.inside <- 40      # length inside protected area 
length.outside <- 20     # length outside protected area

starting.adults <- 2500  # adult tsetse at day 1 / km2
starting.pupae <- 2500   # pupae at day 1 / km2

days <- 3*365*4          # 0.25 day time steps days to run simulation - allow enough time for population to reach carrying capacity
daysseq <- seq(2,days,1) # run 0.5 days
#**********************************************************
#*****************Model parameters*************************
tsetse_params <- function(   l = (1/10)/4          ## proportion of adults producing larvae per day
                             , p = (1/30)/4        ## proportion of pupae that emerge as adults  
                             , md.p = (0.00006)/4  ## pupal density-dependent mortality
                             , mup = (0.006)/4     ## daily probability of pupal mortality
                             , mu.a.b =(0.03)/4    ## 1/40 0.018 ## adult mortality
                             , alph = (0.5)/4      ## adult diffusion coefficient
                             , mo.val=(0.06)/4     ## adult mortality rate outside protected area
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
  vlt <- vl[1,]        # take one row 
  plt <- pl[1,]
  trans<-seq(-19.5,10,0.5)*1000
  return(cbind.data.frame(trans,vlt,plt))
}
#**************************************************************************

#****************function for fitting*********************************
fit <- function(par,dat=data){
    # parameters for the estimation routine
    a <- par[1] 
    b <- par[2]

    sim.res <- tsetse_run_mod(tsetse_params(mo.val=a
                                            ,md.p=b
                                            ,l=(1/12)/4# l=1/12
                                            ,alph=(0.25*2)/4
                                            ))
    out <- output_func(sim.res)
    # filter data that contains distances where data are available
    out.dat <- out[out$trans %in% dat$dist,]
    # evaluate residual
   diff <-  sapply(unique(dat$dist),function(x){
      subdat <- dat[dat$dist %in% x,]
      p <- out.dat[out.dat$trans %in% x,]
      p <- rep(p$vlt,length(subdat$dist))
      return(log10(p+1) - log10(subdat$tot*100+1))
    })
    
    ssqres <- unlist(diff)#log10(out.dat$vlt+1) - log10(dat$tot*100+1)
    return(ssqres)
}

#********************************************************************
#************************read in data******************************
data <- read.csv("tsetseCounts.csv")
data <- newdat[!newdat$dist < -5000,]

#********starting parameters***********
params <- c(0.1/4,0.0002/4)
#********Fit*******************

fm <- nls.lm(par=params,fn=fit,lower=c(0.05/4,.000001/4),upper=c(0.9/4,0.001/4),control = nls.lm.control(maxfev = integer(), 
#**********************************************************                                                                                         maxiter = 50, nprint=100))

summary(fm)
coef(fm)  

beep()


