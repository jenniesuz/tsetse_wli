
#***********Model of tsetse population and trypanosome transmission dynamics*********
#****************Function to create reflective boundaries********************
#********for tsetse movement************
shift.matrix <- function( mat=adults ) {     # matrix of numbers of tsetse
  
  matN <- rbind( mat[1,], mat[-nrow(mat),] ) # shift matrix down 1 cell - append the first row onto the top - reflective boundary
  matS <- rbind( mat[-1,], mat[nrow(mat),] ) # shift matrix up 1 cell
  matE <- cbind( mat[,-1], mat[,ncol(mat)] ) # shift matrix left 1 cell
  matW <- cbind( mat[,1], mat[,-ncol(mat)] ) # shift matrix right 1 cell
  
  return(array( data=c(matN,matS,matE,matW), # stack the matrices into an array
                
                dim=c(length(matN[,1]), length(matN[1,]),4),
                
                dimnames=list(NULL,NULL,c("N","S","E","W")))) # name each in direction shifted
}
#***************************************************************************
#************************Model***********************************************
tsetse.func <- function(s.tsetse                    # matrix of susceptible adult tsetse 
                        ,s.tsetse.shift             # shifted matrix of susceptible adult tsetse         
                        ,e.tsetse1                   # exposed
                        ,e.tsetse1.shift
                        ,e.tsetse2                   # exposed
                        ,e.tsetse2.shift
                        ,e.tsetse3                   # exposed
                        ,e.tsetse3.shift
                        ,i.tsetse                   # infectious
                        ,i.tsetse.shift
                        ,g.tsetse                   # non-teneral - for brucei
                        ,g.tsetse.shift  
                        ,pupae                      # matrix of pupae                   
                        ,larv=l                     ###***see tsetse params in r_2_params.R*****
                        ,pup=p
                        ,mud.p=md.p
                        ,mu.p=mup
                        ,mu.b=mu.a.b               
                        ,alpha=alph
                        ,bite=b
                        ,prob.vtinf=pvit
                        ,prob.vginf=pvig
                        ,t.inc=ti
                        ,h.mort=hm
                        ,h.wan=hw
                        ,prob.hinf = phi
                        ,h.inc = hi
                        ,h.rec= hr
                        ,s.hosts                   # matrix of susceptible hosts
                        ,e.hosts                   # exposed hosts
                        ,i.hosts                   # infected
                        ,r.hosts
                        ,mu.a.out = mu.out
) {  
  
  hosts <- (s.hosts + e.hosts + i.hosts + r.hosts)      # total hosts in each grid cell
  adults <- (s.tsetse + e.tsetse1 + e.tsetse2 + e.tsetse3 + i.tsetse + g.tsetse) # total adult vectors in each grid cell
  alpha <- alpha/4                                      # tsetse movement per day divided by number of directions may move in (n,s,e,w)
  
  # susceptible tsetse:
  diff.sus <- s.tsetse + 
    alpha * s.tsetse.shift[,,"N"] +
    alpha * s.tsetse.shift[,,"S"] +
    alpha * s.tsetse.shift[,,"E"] +
    alpha * s.tsetse.shift[,,"W"] -
    alpha*4*s.tsetse
  
  # non-teneral resistant tsetse:
  diff.g <- g.tsetse + 
    alpha * g.tsetse.shift[,,"N"] +
    alpha * g.tsetse.shift[,,"S"] +
    alpha * g.tsetse.shift[,,"E"] +
    alpha * g.tsetse.shift[,,"W"] -
    alpha*4*g.tsetse
  
  
    st <- (diff.sus - diff.sus*mu.b)*(1 - mu.a.out) - bite*(diff.sus - diff.sus*mu.b)*(1 - mu.a.out) + pupae*pup


  # exposed tsetse:
  diff.exp1 <- e.tsetse1 + 
    alpha * e.tsetse1.shift[,,"N"] +
    alpha * e.tsetse1.shift[,,"S"] +
    alpha * e.tsetse1.shift[,,"E"] +
    alpha * e.tsetse1.shift[,,"W"] -
    alpha*4*e.tsetse1
  
  et1 <- (diff.exp1 - diff.exp1*mu.b)*(1 - mu.a.out)*(1 - t.inc*3) + bite*prob.vtinf*(i.hosts/hosts)*(diff.sus - diff.sus*mu.b)*(1 - mu.a.out) + bite*prob.vginf*(i.hosts/hosts)*(diff.g - diff.g*mu.b)*(1 - mu.a.out)
  
  diff.exp2 <- e.tsetse2 + 
    alpha * e.tsetse2.shift[,,"N"] +
    alpha * e.tsetse2.shift[,,"S"] +
    alpha * e.tsetse2.shift[,,"E"] +
    alpha * e.tsetse2.shift[,,"W"] -
    alpha*4*e.tsetse2
  
  et2 <- (diff.exp2 - diff.exp2*mu.b)*(1 - mu.a.out)*(1 - t.inc*3) + 3*t.inc*(diff.exp1 - diff.exp1*mu.b)*(1 - mu.a.out)
  
  
  diff.exp3 <- e.tsetse3 + 
    alpha * e.tsetse3.shift[,,"N"] +
    alpha * e.tsetse3.shift[,,"S"] +
    alpha * e.tsetse3.shift[,,"E"] +
    alpha * e.tsetse3.shift[,,"W"] -
    alpha*4*e.tsetse3

    et3 <- (diff.exp3 - diff.exp3*mu.b)*(1 - mu.a.out)*(1 - t.inc*3) + 3*t.inc*(diff.exp2 - diff.exp2*mu.b)*(1 - mu.a.out)
  
    
  # infected tsetse:
  diff.inf <- i.tsetse + 
    alpha * i.tsetse.shift[,,"N"] +
    alpha * i.tsetse.shift[,,"S"] +
    alpha * i.tsetse.shift[,,"E"] +
    alpha * i.tsetse.shift[,,"W"] -
    alpha*4*i.tsetse
  
  it <- (diff.inf - diff.inf*mu.b)*(1 - mu.a.out) + 3*t.inc*(diff.exp3 - diff.exp3*mu.b)*(1 - mu.a.out)
  
  gt <- (diff.g - diff.g*mu.b)*(1 - mu.a.out) + bite*(diff.sus - diff.sus*mu.b)*(1 - mu.a.out) - bite*prob.vtinf*(i.hosts/hosts)*(diff.sus - diff.sus*mu.b)*(1 - mu.a.out) - bite*prob.vginf*(i.hosts/hosts)*(diff.g - diff.g*mu.b)*(1 - mu.a.out)

  
  p <- pupae*(1 - pup)*(1 - mu.p)
  p <- p - p*p*mud.p + adults*larv/2                            
  

  new_s.hosts <- h.mort*hosts + s.hosts*(1 - h.mort) - ((1 - ((1 - 1/hosts)^(i.tsetse*bite*prob.hinf))))*s.hosts*(1 - h.mort) + r.hosts*(1 - h.mort)*(h.wan) 
  new_e.hosts <- e.hosts*(1 - h.mort)*(1 - h.inc) + ((1 - ((1 - 1/hosts)^(i.tsetse*bite*prob.hinf))))*s.hosts*(1 - h.mort) 
  new_i.hosts <- i.hosts*(1 - h.mort)*(1 - h.rec) + h.inc*e.hosts*(1 - h.mort)
  new_r.hosts <- r.hosts*(1 - h.mort)*(1 - h.wan)  + h.rec*i.hosts*(1 - h.mort)

  return(list(s.tsetse=st
              ,e.tsetse1=et1
              ,e.tsetse2=et2
              ,e.tsetse3=et3
              ,i.tsetse=it
              ,g.tsetse=gt
              ,pupae=p
              ,s.hosts=new_s.hosts
              ,e.hosts=new_e.hosts
              ,i.hosts=new_i.hosts
              ,r.hosts=new_r.hosts
            ))
}
#*************************************************************************
#*****************Run Model********************************************
tsetse_run_mod <- function(parms=tsetse_params())with(c(parms), {
  
  mu.out <- matrix(nrow=area.width,ncol=area.length)            
  for( i in 1:length(mu.out[,1])){
    mu.out[i,] <- rep(c(rep(0,length.inside)                     # proportion cattle inside (0)
                             ,rep(mo.val,length.outside))) 
  }

  
  s.tsetse_output <- array(0,dim=c(area.width,area.length,days))       # array to store results - create new for each run
  e1.tsetse_output <- array(0,dim=c(area.width,area.length,days))
  e2.tsetse_output <- array(0,dim=c(area.width,area.length,days))
  e3.tsetse_output <- array(0,dim=c(area.width,area.length,days))
  i.tsetse_output <- array(0,dim=c(area.width,area.length,days))
  g.tsetse_output <- array(0,dim=c(area.width,area.length,days))
  pupae_output <- array(0,dim=c(area.width,area.length,days))
  
  s.tsetse_output[,,1] <- matrix(rep(starting.adults,area.length*area.width),
                                 nrow=area.width,ncol=area.length) 
  
  
  pupae_output[,,1] <- matrix(rep(starting.pupae,area.length*area.width),
                              nrow=area.width,ncol=area.length)
  
  s.hosts_output <- array(0,dim=c(area.width,area.length,days))
  e.hosts_output <- array(0,dim=c(area.width,area.length,days))
  i.hosts_output <- array(0,dim=c(area.width,area.length,days))
  r.hosts_output <- array(0,dim=c(area.width,area.length,days))
  
  sh <- matrix(nrow=area.width,ncol=area.length)
  for(i in 1:length(sh[,1])){
    sh[i,] <- rep(c(rep(starting.hosts.inside,length.inside)
                    ,rep(starting.hosts.outside,length.outside)))
  }
  s.hosts_output[,,1] <- sh
  
  
  i.hosts_output[,,1] <- matrix(rep(1,area.length*area.width),
                                nrow=area.width,ncol=area.length) 
  output <- list(s.tsetse_output
                 ,e1.tsetse_output
                 ,e2.tsetse_output
                 ,e3.tsetse_output
                 ,i.tsetse_output
                 ,g.tsetse_output
                 ,pupae_output
                 ,s.hosts_output
                 ,e.hosts_output
                 ,i.hosts_output
                 ,r.hosts_output)
  
  for(i in 2:days){
    s.tsetse_output <- output[[1]]
    e1.tsetse_output <- output[[2]]
    e2.tsetse_output <- output[[3]]
    e3.tsetse_output <- output[[4]]
    i.tsetse_output <- output[[5]]
    g.tsetse_output <- output[[6]]
    pupae_output <- output[[7]]
    s.hosts_output <- output[[8]]
    e.hosts_output <- output[[9]]
    i.hosts_output <- output[[10]]
    r.hosts_output <- output[[11]]

    s.tsetseoutput.shift <- shift.matrix(s.tsetse_output[,,i-1])
    e.tsetseoutput1.shift <- shift.matrix(e1.tsetse_output[,,i-1])
    e.tsetseoutput2.shift <- shift.matrix(e2.tsetse_output[,,i-1])
    e.tsetseoutput3.shift <- shift.matrix(e3.tsetse_output[,,i-1])
    i.tsetseoutput.shift <- shift.matrix(i.tsetse_output[,,i-1])
    g.tsetseoutput.shift <- shift.matrix(g.tsetse_output[,,i-1])
    
    temp.list <- tsetse.func(s.tsetse = s.tsetse_output[,,i-1]                # matrix of susceptible adult tsetse 
                             ,s.tsetse.shift = s.tsetseoutput.shift           # shifted matrix of susceptible adult tsetse         
                             ,e.tsetse1  =  e1.tsetse_output[,,i-1]             # exposed
                             ,e.tsetse1.shift =  e.tsetseoutput1.shift
                             ,e.tsetse2  =  e2.tsetse_output[,,i-1]             # exposed
                             ,e.tsetse2.shift =  e.tsetseoutput2.shift  
                             ,e.tsetse3  =  e3.tsetse_output[,,i-1]             # exposed
                             ,e.tsetse3.shift =  e.tsetseoutput3.shift  
                             ,i.tsetse   = i.tsetse_output[,,i-1]                # infectious
                             ,i.tsetse.shift = i.tsetseoutput.shift  
                             ,g.tsetse  = g.tsetse_output[,,i-1]                   # non-teneral
                             ,g.tsetse.shift  = g.tsetseoutput.shift  
                             ,pupae=pupae_output[,,i-1]                  
                             ,s.hosts  = s.hosts_output[,,i-1]                 # susceptible hosts
                             ,e.hosts   = e.hosts_output[,,i-1]                 # exposed hosts
                             ,i.hosts   = i.hosts_output[,,i-1]                 # infected
                             ,r.hosts  = r.hosts_output[,,i-1]                  # recovered
                             ,larv=l
                             ,pup=p
                             ,mud.p=md.p
                             ,mu.p=mup
                             ,mu.b=mu.a.b               # matrix of adult mortalities from mu.a.b and mu.a.u
                             ,alpha=alph
                             ,bite=b
                             ,prob.vtinf=pvit
                             ,prob.vginf=pvig
                             ,t.inc=ti
                             ,h.mort=hm
                             ,h.wan=hw
                             ,prob.hinf = phi
                             ,h.inc = hi
                             ,h.rec= hr
                             ,mu.a.out = mu.out
)
    
    output[[1]][,,i] <- temp.list[[1]]
    output[[2]][,,i] <- temp.list[[2]]
    output[[3]][,,i] <- temp.list[[3]]
    output[[4]][,,i] <- temp.list[[4]]
    output[[5]][,,i] <- temp.list[[5]]
    output[[6]][,,i] <- temp.list[[6]]
    output[[7]][,,i] <- temp.list[[7]]
    output[[8]][,,i] <- temp.list[[8]]
    output[[9]][,,i] <- temp.list[[9]]
    output[[10]][,,i] <- temp.list[[10]]
    output[[11]][,,i] <- temp.list[[11]]

  }
  
  #return( output_flies[[1]][,,days])
  return(output)
})
#*************************************************************
#*****************summarise_output***************

output_func <- function(mod_output){
  #**********subset host data***********
sh<-mod_output[[8]]
eh<-mod_output[[9]]     # each host class
rh<-mod_output[[11]]
ih<-mod_output[[10]]

ihl <- ih[,,days]  # take the last time point (end of model run)
ehl <- eh[,,days]
shl <- sh[,,days]
rhl <- rh[,,days]

ihlt <- ihl[1,]   # take one row of the 2D world - one transect from inside to outside protected area
rhlt <- rhl[1,]
ehlt <- ehl[1,]
shlt <- shl[1,]

#********************************************************

#*****************subset vector data*************
sv<-mod_output[[1]]     # each vector
ev1<-mod_output[[2]]
ev2<-mod_output[[3]]
ev3<-mod_output[[4]]
iv<-mod_output[[5]]
gv<-mod_output[[6]]

svl <- sv[,,days]  # take the last time point (end of model run)
evl1 <- ev1[,,days]
evl2 <- ev2[,,days]
evl3 <- ev3[,,days]
ivl <- iv[,,days]
gvl <- gv[,,days]

svlt <- svl[1,]        # take one row of the 2D world - one transect from inside to outside protected area
evlt1 <- evl1[1,]
evlt2 <- evl2[1,]
evlt3 <- evl3[1,]
ivlt <- ivl[1,]
gvlt <- gvl[1,]

return(cbind.data.frame(svlt,evlt1,evlt2,evlt3,ivlt,gvlt,ihlt,rhlt,ehlt,shlt))
}
