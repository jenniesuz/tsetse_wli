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
                             , p = (1/30)/2      ## proportion of pupae that emerge as adults  
                             , md.p = (0.00006)/2  ## 0.00006 pupal density-dependent mortality
                             , mup = (0.006)/2     ## daily probability of pupal mortality
                             , mu.a.b =(0.03)/2  ## 1/40 0.018 ## adult mortality
                             , alph = (0.5)/2   ## adult diffusion coefficient
                             , mo.val=(0.06)/2      ## adult mortality rate outside protected area
)       
return(as.list(environment()))
#*************************************************************
#*****************summarise_output***************

output_func <- function(mod_output){
  
  #*****************subset vector data*************
  v<-mod_output[[1]]     # each vector
  vl <- v[,,days]  # take the last time point (end of model run)
  trans<-seq(-19,10,1)*1000
  vl <- matrix(vl,ncol=4,byrow=T)
  return(cbind.data.frame(trans,vl))
}
