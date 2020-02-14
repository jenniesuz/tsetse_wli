#***********************Model of tsetse population dynamics*******************

#****************Function to create reflecting boundaries********************
#********for tsetse movement************
shift.matrix <- function( mat=adults ) {     # matrix of numbers of tsetse
  
  matN <- rbind( mat[1,], mat[-nrow(mat),] ) # shift matrix down 1 cell - append the first row onto the top 
  matS <- rbind( mat[-1,], mat[nrow(mat),] ) # shift matrix up 1 cell
  matE <- cbind( mat[,-1], mat[,ncol(mat)] ) # shift matrix left 1 cell
  matW <- cbind( mat[,1], mat[,-ncol(mat)] ) # shift matrix right 1 cell
  
  return(array( data=c(matN,matS,matE,matW), # stack the matrices into an array
                
                dim=c(length(matN[,1]), length(matN[1,]),4),
                
                dimnames=list(NULL,NULL,c("N","S","E","W")))) # name each in direction shifted
}
#***************************************************************************
#************************Model***********************************************
tsetse.func <- function(tsetse                    # matrix of adult tsetse 
                        ,tsetse.shift             # shifted matrix of adult tsetse         
                        ,pupae                    # matrix of pupae                  
                        ,larv=l                   ###***see tsetse params in r_1_params.R*****
                        ,pup=p
                        ,mud.p=md.p
                        ,mu.p=mup
                        ,mu.b=mu.a.b               
                        ,alpha=alph
                        ,mu.a.out=mo
) {  
  
  alpha <- alpha/4                                 # tsetse movement per day divided by number of directions may move in (n,s,e,w)

  diff.tsetse <- tsetse + 
    alpha * tsetse.shift[,,"N"] +
    alpha * tsetse.shift[,,"S"] +
    alpha * tsetse.shift[,,"E"] +
    alpha * tsetse.shift[,,"W"] -
    alpha*4*tsetse

    adults <- (diff.tsetse)*(1 - mu.b)*(1 - mu.a.out) + pupae*pup
    pupae <- (pupae)*(1 - pup)*(1 - mu.p)*(1 - pupae*(1 - pup)*(1 - mu.p)*mud.p) + adults*larv/2

  return(list(tsetse=adults
              ,pupae=pupae))
}
#*************************************************************************
#*****************Run Model********************************************
tsetse_run_mod <- function(parms=tsetse_params())with(c(parms), {
  
  mo <- matrix(nrow=area.width,ncol=area.length)            
  for( i in 1:length(mo[,1])){
    mo[i,] <- rep(c(rep(0,length.inside)                    
                             ,rep(mo.val,length.outside)))  
  }
  
  
  tsetse_output <- array(0,dim=c(area.width,area.length,days))       # array to store results - create new for each run
  pupae_output <- array(0,dim=c(area.width,area.length,days))
  
  tsetse_output[,,1] <- matrix(rep(starting.adults,area.length*area.width),
                                 nrow=area.width,ncol=area.length) 
  
  
  pupae_output[,,1] <- matrix(rep(starting.pupae,area.length*area.width),
                              nrow=area.width,ncol=area.length) 
  
  output <- list(tsetse_output
                 ,pupae_output)
  
  diff.coef <- matrix(rep(alph,area.length*area.width),nrow=area.width,ncol=area.length)
  
  for(i in 2:days){
    tsetse_output <- output[[1]]
    pupae_output <- output[[2]]

    tsetseoutput.shift <- shift.matrix(tsetse_output[,,i-1])
    
    temp.list <- tsetse.func(tsetse = tsetse_output[,,i-1]                # matrix of susceptible adult tsetse 
                             ,tsetse.shift = tsetseoutput.shift           # shifted matrix of susceptible adult tsetse         
                             ,pupae=pupae_output[,,i-1]                  
                             ,larv=l
                             ,pup=p
                             ,mud.p=md.p
                             ,mu.p=mup
                             ,mu.b=mu.a.b               # matrix of adult mortalities from mu.a.b and mu.a.u
                             ,alpha=alph
                             ,mu.a.out=mo
)
    
    output[[1]][,,i] <- temp.list[[1]]
    output[[2]][,,i] <- temp.list[[2]]

  }

  return(output)
})
#*************************************************************

