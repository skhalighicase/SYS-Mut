## Last Modification 6/6/2020
### The core of the model

SYS_Model<-function(CNV_Meth_Pa,RNA_GII,RNA_TGG,PGs_Weights,Miss_Com_Indxx,weightsGI){
 
  library(rjags)
  library(R2jags)
  library(coda)
  
  if(!is.matrix(CNV_Meth_Pa)){
    CNV_Meth_Pa<-t(as.matrix(CNV_Meth_Pa))
  }
  datastruct <-list("GI"=as.vector(RNA_GII),
                    "PG"= t(CNV_Meth_Pa),
                    "n"=length(RNA_TGG),
                    "K"=nrow(CNV_Meth_Pa),
                    "y"=as.vector(RNA_TGG),
                    "weightsP"=PGs_Weights,
                    "R"=Miss_Com_Indxx,
                    "weightsGI"=weightsGI)
  niter=15#1500 
  nburnin =10#750
  nchains =3
  
  SYSmodel <- "model{
  
  ######################
  for( i in 1 : n ) {
  y[i] ~ dnorm( mu[i] , sigma )
  mu[i] <- GI[i] * (Mutation[R[i]])+ PG[i,1:K] %*% BetaP[1:K]
  # loglik[i] <- logdensity.norm(y[i], mu[i], sigma)
  }
  ##
  for (j in 1:K){
  BetaP[j] ~ dnorm(weightsP[j], lambdaP)
  }
  
  ### Priors
  lambdaP ~ dgamma(0.01,0.01)
  
  ### Priors
  sigma ~ dgamma(0.01,0.01)
  
  ## Depending on the number of mutation we can increase
  for (S in 1:2){
  Mutation[S] ~ dnorm(Mu_Muta[S] , SigMut)
  }
  ### Priors
  Mu_Muta[1]~dnorm(weightsGI,BigSigMut)
  Mu_Muta[2]~dnorm(weightsGI,BigSigMut) ## WeightsGI is an input
  
  BigSigMut ~ dgamma(0.01,0.01)
   
  SigMut~dgamma( sh , ra )
 # parameterized by mean (m) and standard deviation (sd)
  sh <- pow(m,2) / pow(sd,2)
  ra <-     m    / pow(sd,2)
  m ~ dunif(0,100)
  sd ~ dunif(0,100)
 
  Invsigma <- 1/sqrt(sigma)
  }"
  
  #Initialization
  inits1 <- list(BetaP=rnorm((nrow(CNV_Meth_Pa))*1,mean=0,sd=dgamma(1,2,0.5)), Mutation=rnorm(2*1,mean=0,sd=dgamma(1,2,0.5)))
  inits2 <- list(BetaP=rnorm((nrow(CNV_Meth_Pa))*1,mean=0,sd=dgamma(1,2,0.5)), Mutation=rnorm(2*1,mean=0,sd=dgamma(1,2,0.5)))
  inits3 <- list(BetaP=rnorm((nrow(CNV_Meth_Pa))*1,mean=0,sd=dgamma(1,2,0.5)), Mutation=rnorm(2*1,mean=0,sd=dgamma(1,2,0.5)))
  angell_inits <- list(inits1, inits2, inits3)
  
  ## modelend
  bayes.fit <- jags(data=datastruct,
                    parameters.to.save = c("Mutation","Invsigma"),
                    model.file = textConnection(SYSmodel),
                    inits=angell_inits,
                    n.chains = nchains,
                    n.iter = niter,
                    n.burnin = nburnin,
                    n.thin = 1,
                    DIC=1)
  
  # print(bayes.fit)
  # plot(bayes.fit)
  # browser()
  
  Mu_Phi2<-bayes.fit$BUGSoutput$mean$Mutation
  Var_Phi2<-bayes.fit$BUGSoutput$sd$Mutation
  Var_G2<-bayes.fit$BUGSoutput$mean$Invsigma
  PostMutation=data.frame(Mu_Phi2=Mu_Phi2,Var_G2=Var_G2,Var_Phi2=Var_Phi2)
  return(PostMutation)
}