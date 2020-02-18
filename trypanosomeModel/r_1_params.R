

#******************Model constants************************
area.length = 30         # length and width (km) of total area
area.width = 4
length.inside <- 20      # length inside protected area
length.outside <- 10     # length outside protected area

starting.adults <- 5000  # adult tsetse at day 1 / km2
starting.pupae <- 5000  # pupae at day 1 / km2
starting.hosts.inside <- 30
starting.hosts.outside <- 30   # hosts at day 1 / km2 #****if host density constant across area*****

days <- 3*365         # days to run simulation - allow enough time for population to reach carrying capacity and tryps equilibrium
daysseq <- seq(2,days,1)
#**********************************************************
#*****************Model parameters*************************
tsetse_params <- function(   l = (1/10)          ## proportion of adults producing larvae per day
                             , p = (1/30)        ## proportion of pupae that emerge as adults  
                             , md.p = (0.00006)  ## 0.00006 pupal density-dependent mortality
                             , mup = (0.006)     ## daily probability of pupal mortality
                             , mu.a.b = (0.03)   ## 1/40 0.018 ## adult mortality
                             , alph = (0.9)      ## adult diffusion coefficient
                             , b = (1/4)         ## proportion biting
                             , pvit = 0.065      ## probability of vector infection
                             , pvig = 0.065
                             , ti = (1/25)        ## extrinsic incubation 
                             , hm = (1/(10*365))  ## host birth and death 
                             , hw = (1/50)        ## waining immunity 
                             , phi = 0.62         ## probability of vector to host transmission 
                             , hi = (1/12)        ## incubation 
                             , hr = (1/50)        ## host recovery 
                             , mo.val=(0.2) )       
return(as.list(environment()))
#*************************************************************
