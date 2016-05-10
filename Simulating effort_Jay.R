
  ###########################################################################################
  # Simulation of fishing effort by using change of effort over time
  # Change of effort was described by a linear model with slope s and intercept i
  # s and i is random by using Monte Carlo method
  # Change of effort was integrated to effort over time
  #
  # Original from Carruthers et al. (2012) - function built using equations provided in appendix A1
  # modified by Yi-Jay Chang (included analytical integration)
  ###########################################################################################
  
  
  getFhist<-function(nsim,Esd,nyears,dFmin,dFmax,i){
  
  ne<-nsim*10                                           # Number of simulated effort datasets
  dVfinal<-rnorm(ne,0,1)                                # Sample the final gradient in effort
  s<-(dVfinal-i)/nyears                                 # Derive slope to get there from intercept
  s<-array(s,dim=c(ne,nyears))                          # Slope array
  i<-array(i,dim=c(ne,nyears))                          # Intercept array
  y<-array(rep(1:nyears,each=ne),dim=c(ne,nyears))      # Year array
  dV<-s*y+i                                             # Change in effort
  
  E<-array(NA,dim=c(ne,nyears))                         # Define total effort array
  
  # old code -------------------------------------------------
  #E[,1]<-dE[,1]
  #for(y in 2:nyears){
  #  E[,y]<-apply(dE[,1:y],1,sum)
  #}
  # old code -------------------------------------------------
  
  # new code by Jay-------------------------------------------
  E<-0.5*s*y^2+i*y
  # new code by Jay-------------------------------------------
  
  E<-E/array(apply(E,1,mean),dim=c(ne,nyears))          # Standardise Effort to average 1
  cond<-apply(E,1,min)>0
  pos<-(1:ne)[cond]
  pos<-pos[1:nsim]
    
  E<-E[pos,]                                            # Sample only those without negative effort
  Emu<--0.5*Esd^2
  Eerr<-array(exp(rnorm(nyears*nsim,rep(Emu,nyears),rep(Esd,nyears))),c(nsim,nyears))
  E*Eerr
  }

 nyears<-60
 nsim<-30

 Edat<-getFhist(nsim=nsim,Esd=0.1,nyears=nyears,dFmin=-0.1,dFmax=0.1,i=0.1)
 matplot(t(Edat),type='l',ylim=c(0,max(Edat)),ylab="Simulated effort",xlab="Year")
