#' Routine for simulating the dynamic R (different routines for different models)
#'
#' @param R_sim a matrix of dimension (unit.samp+1)*B with the first row giving the current R values.
#' @param Suff_R a matrix of dimesion B*nsuff where nsuff is the number of sufficient statistics for updating parameters
#' @param dyn_model The model for dynamic R, currently only RW (random walk) is implemented
#' @param param_dynamicR List with parameters related to dynamic model for R. For RW, mean and standard deviation in inverse  Gamma distribution for variance in dynamic R distribution. Default is (1,1).
#' @param mu_AMP R0 The amplification factor AMP is assumed to follow a lognormal distribution with expectation mu_AMP. Default is 2.8.
#' @param sd_AMP The standard deviation in the lognormal distribution for AMP. Default is 0.0 in which case no updating of AMP is performed.
#' @param prob_h Probability for hospitalisation. Default to 0.036.
#' @param alpha_hosp First parameter in beta-binomial distribution for hospitalization. If NULL, the model reduces to a binomial distribution with 'prob_h' being the probability for hospitalization. In general, the second (beta) parameter is calculated as 'alpha_hosp * (1 - prob_h) / prob_h'. Default is NULL.
#' @param alpha_test First parameter in beta-binomial distribution for testdata. If NULL, the model reduces to a binomial distribution with 'prob_h' being the probability for hospitalization. In general, the second (beta) parameter is calculated as 'alpha_hosp * (1 - prob_h) / prob_h'. Default is NULL.
#' @param pi_test Detection probability for test data
#' @param lambda Expected length from sickness to hospitalization. Default to 9.08.
#' @param lambda_size Size parameter in the Negative binomial distribution for length until hospitalization. Default is NULL in which case a Poisson distribution is used.
#' @param lambda_ch If 'lambda' and 'lambda_size' are vectors, this correspond to the dates for changes in these parameters
#' @return A list containing particles of latent variables ('latent'), particles of parameters ('param') and estimate of marginal log-likelihood.
smc_sim_R_RW = function(R_sim,Suff_R,unit.samp,param_R)
{
  B = ncol(R_sim)
  sig.R = 1/sqrt(rgamma(B,shape=param_R[1]+0.5*Suff_R[,1],rate=param_R[2]+0.5*Suff_R[,2]))
  for (i in 1:unit.samp)
  {
    eps = rnorm(B, 0, sig.R)
    R_sim[i+1,]=R_sim[i,]*exp(eps)    #Note first new value in column 2! -> ZERO ?### ALF
    
    Suff_R[,2] = Suff_R[,2] + eps * eps                  #Sufficient statistic for sig.R
    Suff_R[,1] = Suff_R[,1] + 1
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=matrix(sig.R,ncol=1))
}

smc_sim_R_AR = function(R_sim,Suff_R,unit.samp,param_R)
{
  sampvec = function(p)
  {
    sample(1:length(p),1,prob=p)
  }
  #Suff1 is number of timepoints
  #Suff2 is sum of log(R)
  #Suff3 is sum of lagged log(R)
  #Suff4 is sum of log(R) squared
  #Suff5 is sum of log(R) and lagged log(R)
  #Suff6 is sum of lagged log(R) squared
  mu0 = param_R[1]
  kappa0 = param_R[2]
  alpha = param_R[3]
  beta = param_R[4]
  B = ncol(R_sim)
  nadiscr = 25
  adiscr = seq(0.01,0.99,length=nadiscr)
  p.a = matrix(NA,nrow=B,ncol=nadiscr)
  for(i in 1:nadiscr)
  {
    a = adiscr[i]
    num = kappa0+Suff_R[,1]*(1-a)^2
    mu.hat = (kappa0*mu0+(1-a)*(Suff_R[,2]-a*Suff_R[,3]))/num
    g = kappa0*mu0^2+Suff_R[,4]-2.0*a*Suff_R[,6]+a^2*Suff_R[,5]-num*mu.hat^2
    if(min(g)<0)
    {
      print("Something is wrong")
      show(a)
      show(mu.hat)
      show(g)
    }
    p.a[,i] = -0.5*log(num)-(alpha+0.5*Suff_R[,1])*log(param_R[2]+0.5*g)
  }
  p.a = t(apply(p.a,1,smc_trans_to_prob))
  foo = apply(p.a,1,sampvec)
  a.sim = adiscr[foo]
  num = kappa0+Suff_R[,1]*(1-a.sim)^2
  mu.hat = (kappa0*mu0+(1-a.sim)*(Suff_R[,2]-a.sim*Suff_R[,3]))/num
  g = kappa0*mu0^2+Suff_R[,4]-2.0*a.sim*Suff_R[,6]+a.sim^2*Suff_R[,5]-num*mu.hat^2
  if(min(g)<0)
  {
    print("Something is wrong")
    show(a.sim)
    show(mu.hat)
    show(g)
  }
 
  tau.sim = rgamma(B,alpha+0.5*Suff_R[,1],beta+0.5*g)
  sig.sim = 1/sqrt(tau.sim)
  mu.sim = rnorm(B,mu.hat,sig.sim/sqrt(num))
  for (i in 1:unit.samp)
  {
    eps = rnorm(B, 0, sig.sim)
    R_sim[i+1,]=exp(mu.sim+a.sim*(log(R_sim[i,])-mu.sim)+eps)    #Note first new value in column 2! -> ZERO ?### ALF
    
    Suff_R[,1] = Suff_R[,1] + 1
    Suff_R[,2] = Suff_R[,2] + log(R_sim[i+1,])
    Suff_R[,3] = Suff_R[,3] + log(R_sim[i,])
    Suff_R[,4] = Suff_R[,4] + log(R_sim[i+1,])^2
    Suff_R[,5] = Suff_R[,5] + log(R_sim[i,])^2
    Suff_R[,6] = Suff_R[,6] + log(R_sim[i,])*log(R_sim[i+1,])              
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=cbind(a.sim,sig.sim,mu.sim))
}

smc_sim_R_AR_mu_0 = function(R_sim,Suff_R,unit.samp,param_R)
{
  sampvec = function(p)
  {
    sample(1:length(p),1,prob=p)
  }
  #Suff1 is number of timepoints
  #Suff2 is sum of log(R)
  #Suff3 is sum of lagged log(R)
  #Suff4 is sum of log(R) squared
  #Suff5 is sum of log(R) and lagged log(R)
  #Suff6 is sum of lagged log(R) squared
  a0 = param_R[1]
  kappa0 = param_R[2]
  alpha = param_R[3]
  beta = param_R[4]
  B = ncol(R_sim)
  num = kappa0+Suff_R[,5]
  a.hat = (kappa0*a0+Suff_R[,6])/num
  g = kappa0*a0^2+Suff_R[,4]-num*a.hat^2
  show(range(a.hat))
  if(min(g)<0)
  {
    print("Something is wrong")
    show(a.hat)
    show(g)
  }
  
  tau.sim = rgamma(B,alpha+0.5*Suff_R[,1],beta+0.5*g)
  sig.sim = 1/sqrt(tau.sim)
  a.sim = rnorm(B,a.hat,sig.sim/sqrt(num))
  mu.sim = rep(0,B)
  for (i in 1:unit.samp)
  {
    eps = rnorm(B, 0, sig.sim)
    R_sim[i+1,]=exp(mu.sim+a.sim*(log(R_sim[i,])-mu.sim)+eps)    #Note first new value in column 2! -> ZERO ?### ALF
    
    Suff_R[,1] = Suff_R[,1] + 1
    Suff_R[,2] = Suff_R[,2] + log(R_sim[i+1,])
    Suff_R[,3] = Suff_R[,3] + log(R_sim[i,])
    Suff_R[,4] = Suff_R[,4] + log(R_sim[i+1,])^2
    Suff_R[,5] = Suff_R[,5] + log(R_sim[i,])^2
    Suff_R[,6] = Suff_R[,6] + log(R_sim[i,])*log(R_sim[i+1,])              
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=cbind(a.sim,sig.sim,mu.sim))
}

smc_sim_R_AR_fix_par = function(R_sim,Suff_R,unit.samp,param_R)
{
  mu.sim = rep(0,B)
  a.sim = rep(param_R[1],B)
  tau.sim = rep(param_R[3]/param_R[4],B)
  sig.sim = 1/sqrt(tau.sim)
  for (i in 1:unit.samp)
  {
    eps = rnorm(B, 0, sig.sim)
    R_sim[i+1,]=exp(mu.sim+a.sim*(log(R_sim[i,])-mu.sim)+eps)    #Note first new value in column 2! -> ZERO ?### ALF
    
    Suff_R[,1] = Suff_R[,1] + 1
    Suff_R[,2] = Suff_R[,2] + log(R_sim[i+1,])
    Suff_R[,3] = Suff_R[,3] + log(R_sim[i,])
    Suff_R[,4] = Suff_R[,4] + log(R_sim[i+1,])^2
    Suff_R[,5] = Suff_R[,5] + log(R_sim[i,])^2
    Suff_R[,6] = Suff_R[,6] + log(R_sim[i,])*log(R_sim[i+1,])              
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=cbind(a.sim,sig.sim,mu.sim))
}

smc_sim_R_AR_a8 = function(R_sim,Suff_R,unit.samp,param_R)
{
  sampvec = function(p)
  {
    sample(1:length(p),1,prob=p)
  }
  #Suff1 is number of timepoints
  #Suff2 is sum of log(R)
  #Suff3 is sum of lagged log(R)
  #Suff4 is sum of log(R) squared
  #Suff5 is sum of log(R) and lagged log(R)
  #Suff6 is sum of lagged log(R) squared
  mu0 = param_R[1]
  kappa0 = param_R[2]
  alpha = param_R[3]
  beta = param_R[4]
  B = ncol(R_sim)
  nadiscr = 25
  adiscr = seq(0.8,0.99,length=nadiscr)
  p.a = matrix(NA,nrow=B,ncol=nadiscr)
  for(i in 1:nadiscr)
  {
    a = adiscr[i]
    num = kappa0+Suff_R[,1]*(1-a)^2
    mu.hat = (kappa0*mu0+(1-a)*(Suff_R[,2]-a*Suff_R[,3]))/num
    g = kappa0*mu0^2+Suff_R[,4]-2.0*a*Suff_R[,6]+a^2*Suff_R[,5]-num*mu.hat^2
    if(min(g)<0)
    {
      print("Something is wrong")
      show(a)
      show(mu.hat)
      show(g)
    }
    p.a[,i] = -0.5*log(num)-(alpha+0.5*Suff_R[,1])*log(param_R[2]+0.5*g)
  }
  p.a = t(apply(p.a,1,smc_trans_to_prob))
  foo = apply(p.a,1,sampvec)
  a.sim = adiscr[foo]
  num = kappa0+Suff_R[,1]*(1-a.sim)^2
  mu.hat = (kappa0*mu0+(1-a.sim)*(Suff_R[,2]-a.sim*Suff_R[,3]))/num
  g = kappa0*mu0^2+Suff_R[,4]-2.0*a.sim*Suff_R[,6]+a.sim^2*Suff_R[,5]-num*mu.hat^2
  if(min(g)<0)
  {
    print("Something is wrong")
    show(a.sim)
    show(mu.hat)
    show(g)
  }
  
  tau.sim = rgamma(B,alpha+0.5*Suff_R[,1],beta+0.5*g)
  sig.sim = 1/sqrt(tau.sim)
  mu.sim = rnorm(B,mu.hat,sig.sim/sqrt(num))
  for (i in 1:unit.samp)
  {
    eps = rnorm(B, 0, sig.sim)
    R_sim[i+1,]=exp(mu.sim+a.sim*(log(R_sim[i,])-mu.sim)+eps)    #Note first new value in column 2! -> ZERO ?### ALF
    
    Suff_R[,1] = Suff_R[,1] + 1
    Suff_R[,2] = Suff_R[,2] + log(R_sim[i+1,])
    Suff_R[,3] = Suff_R[,3] + log(R_sim[i,])
    Suff_R[,4] = Suff_R[,4] + log(R_sim[i+1,])^2
    Suff_R[,5] = Suff_R[,5] + log(R_sim[i,])^2
    Suff_R[,6] = Suff_R[,6] + log(R_sim[i,])*log(R_sim[i+1,])              
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=cbind(a.sim,sig.sim,mu.sim))
}

smc_sim_R_AR_fix_a = function(R_sim,Suff_R,unit.samp,param_R)
{
  sampvec = function(p)
  {
    sample(1:length(p),1,prob=p)
  }
  #Suff1 is number of timepoints
  #Suff2 is sum of log(R)
  #Suff3 is sum of lagged log(R)
  #Suff4 is sum of log(R) squared
  #Suff5 is sum of log(R) and lagged log(R)
  #Suff6 is sum of lagged log(R) squared
  mu0 = param_R[1]
  kappa0 = param_R[2]
  alpha = param_R[3]
  beta = param_R[4]
  B = ncol(R_sim)
  a.sim = rep(0.95,B)
  num = kappa0+Suff_R[,1]*(1-a.sim)^2
  mu.hat = (kappa0*mu0+(1-a.sim)*(Suff_R[,2]-a.sim*Suff_R[,3]))/num
  g = kappa0*mu0^2+Suff_R[,4]-2.0*a.sim*Suff_R[,6]+a.sim^2*Suff_R[,5]-num*mu.hat^2
  if(min(g)<0)
  {
    print("Something is wrong")
    show(a.sim)
    show(mu.hat)
    show(g)
  }
  
  tau.sim = rgamma(B,alpha+0.5*Suff_R[,1],beta+0.5*g)
  sig.sim = 1/sqrt(tau.sim)
  mu.sim = rnorm(B,mu.hat,sig.sim/sqrt(num))
  for (i in 1:unit.samp)
  {
    eps = rnorm(B, 0, sig.sim)
    R_sim[i+1,]=exp(mu.sim+a.sim*(log(R_sim[i,])-mu.sim)+eps)    #Note first new value in column 2! -> ZERO ?### ALF
    
    Suff_R[,1] = Suff_R[,1] + 1
    Suff_R[,2] = Suff_R[,2] + log(R_sim[i+1,])
    Suff_R[,3] = Suff_R[,3] + log(R_sim[i,])
    Suff_R[,4] = Suff_R[,4] + log(R_sim[i+1,])^2
    Suff_R[,5] = Suff_R[,5] + log(R_sim[i,])^2
    Suff_R[,6] = Suff_R[,6] + log(R_sim[i,])*log(R_sim[i+1,])              
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=cbind(a.sim,sig.sim,mu.sim))
}

smc_sim_R_CP = function(R_sim,Suff_R,unit.samp,param_R)
{
  #Suff1 is number of changepoints
  #Suff2 is sum of non-changepoints
  #Suff3 is sum of eps^2
  B = ncol(R_sim)
  phi.c = rbeta(B,param_R[1]+Suff_R[,1],param_R[2]+Suff_R[,2])
  tau.R = rgamma(B,param_R[3]+0.5*Suff_R[,1],param_R[4]+0.5*Suff_R[,3])
  sig.R = 1/sqrt(tau.R)
  for (i in 1:unit.samp)
  {
     cp = rbinom(B,1,phi.c)
     eps = rnorm(B, 0, sig.R)
     R_sim[i+1,]=R_sim[i,]*exp(cp*sig.R*eps)    #Note first new value in column 2! -> ZERO ?### ALF
     #cp = as.numeric(abs(log(R_sim[i+1,])-log(R_sim[i,]))>0.000001)
     #eps = log(R_sim[i+1,])-log(R_sim[i,])
     Suff_R[,1] = Suff_R[,1] + cp
     Suff_R[,2] = Suff_R[,2] + 1-cp
     Suff_R[,3] = Suff_R[,3] + cp * eps * eps                  
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=cbind(phi.c,sig.R))
}

smc_sim_R_TCP = function(R_sim,Suff_R,unit.samp,param_R)
{
  #Suff1 is number of changepoints
  #Suff2 is sum of non-changepoints
  #Suff3 is sum of eps^2
  B = ncol(R_sim)
  phi.c = rbeta(B,param_R[1]+Suff_R[,1],param_R[2]+Suff_R[,2])
  tau.R = rgamma(B,param_R[3]+0.5*Suff_R[,1],param_R[4]+0.5*Suff_R[,3])
  sig.R = 1/sqrt(tau.R)
  for (i in 1:unit.samp)
  {
    cp = rbinom(B,1,phi.c)
    sc = rchisq(B,2)
    eps = rnorm(B, 0, sig.R)
    R_sim[i+1,]=R_sim[i,]*exp(cp*sig.R*eps/sc)    #Note first new value in column 2! -> ZERO ?### ALF
    #cp = as.numeric(abs(log(R_sim[i+1,])-log(R_sim[i,]))>0.000001)
    #eps = log(R_sim[i+1,])-log(R_sim[i,])
    Suff_R[,1] = Suff_R[,1] + cp
    Suff_R[,2] = Suff_R[,2] + 1-cp
    Suff_R[,3] = Suff_R[,3] + cp * eps * eps                  
  }
  list(R_sim=R_sim,Suff_R=Suff_R,par.R=cbind(phi.c,sig.R))
}


