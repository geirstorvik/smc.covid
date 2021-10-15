#' Particle MCMC routine for updating parameters
#'
#' Uses estimates of marginal likelihoods using the smc_loglik routine
#' Make Gaussian random walk moves on log-scale
#' Also assumes flat prior on log-scale
#' So far only using default input to smc_loglik, should be modified
#' 
#' @param start_date starting date for initial simulations
#' @param init_date date until particles are simulated before calibration towards data
#' @param hosp_data data frame containing hospidal data, with components date and N.
#' @param sig.R0 Standard deviation in dynamic random-walk model for. If NULL, static model is used. Default to NULL
#' @param B Number of particles. Default to 100.
#' @param resamp.frac When effective number of samples below resamp.fram, resampling is performed. Default to 1. in which case no resampling is performed.
#' @param seed Seed number. Default to 123.
#' @param parallel Logical variable, if true parallel computing is performed.
#' @param DEBUG If true, some output are written during running
#' @param pmcmc Logical value, if true smc_particle_sim_pmcmc is used for simulation of particles
#' @return data frame with all iterations including parameters and estimated marginal likelihood
#' @export
#' 
smc_PMCMC = function(
      mu_R0=3.3,sd_R0=0.1,mu_AMP=2.8,sd_AMP=0.0,
      dyn_model="RW",
      param_dynamicR = list(alpha=2.4,beta=0.28,phi.a=100,phi.b=900,mu0=1,kappa0=1),
      lambda = 9.08,lsize=3, lch=NULL,
      prob_h = c(0.036,0.022),change_prob=c(as.Date("2020-05-01")), 
      alpha_hosp=8,alpha_test=8,q_test=c(0,0.05,0.080,0.16,0.4,0.3,0.01),pi_test=c(-1.175,7.4e-05),
      start_date=as.Date("2020-02-17"),dynR_date=as.Date("2020-03-08"),
      end_date = as.Date("2020-05-05"),pred_date = NULL,
      hosp_data = hospitalisations_df,  test_data=test_data,
      pop = pop, loc.seed=seeding,mob_mat=mobility_matrices,resdir="./",
      sig.prop=c(0.2,1,1,1,0.01),samp_ind = 1:3,
      B = 20,Bmem=B,B2 = NULL,unit.samp = 3,resamp.frac = 1.1,smo.lag = 25,parallel = TRUE,pmcmc=TRUE,
      VERBOSE=FALSE,DEBUG = FALSE,output_latent = FALSE,seed = 123,Nmcmc=10,N_cores=4)
{
  logit = function(p){log(p)-log(1-p)}
  invlogit = function(x){1/(1+exp(-x))}
  #Priors on transformed scale
  f1 = function(x,a){log(x-a)}
  f2 = function(x,a,b){dnorm(log(x),a,b)}
  f3 = function(x,mu,sigma){dnorm(x,mu,sigma,log=TRUE)}
  if(VERBOSE)
    show(c(alpha_hosp,lsize,param_dynamicR$alpha,param_dynamicR$beta))
  #Estimate of likelihood for initial parameters
  param = rep(NA,5)
  param[1] = log(alpha_hosp-4)
  param[2] = log(param_dynamicR$alpha)
  param[3] = log(param_dynamicR$beta)
  param[4] = pi_test[1]
  param[5] = pi_test[2]
  logp.hat = smc_loglik_mem(
        mu_R0=mu_R0,sd_R0=sd_R0,mu_AMP=mu_AMP,sd_AMP=sd_AMP,
        dyn_model=dyn_model,
        param_dynamicR=param_dynamicR,
        lambda=lambda,lambda_size=lsize,lambda_ch=lch,
        prob_h=prob_h,change_prob=change_prob,alpha_hosp=alpha_hosp,
        q_test=q_test,alpha_test=alpha_test,pi_test=pi_test,
        start_date=start_date,dynR_date=dynR_date,end_date=end_date,
        hosp_data = hosp_data, test_data = test_data, pop = pop, 
        loc.seed = loc.seed, mob_mat = mob_mat, resdir = resdir, 
        B=B,unit.samp=unit.samp,smo.lag=smo.lag,resamp.frac=resamp.frac,parallel=parallel,
        VERBOSE=TRUE,DEBUG=FALSE,pmcmc=pmcmc,output_latent=FALSE,
        seed=seed,N_cores=N_cores)$logp.hat
  pri = f1(alpha_hosp,4)+f2(param_dynamicR$alpha,0.1,1)+f2(param_dynamicR$beta,0.1,1)+
        f3(pi_test[1],0,5)+f3(pi_test[2],0,0.1)
  if(VERBOSE)
    cat("alpha=",alpha_hosp,"alpha=",param_dynamicR$alpha,"beta=",param_dynamicR$beta,"pi=",pi_test,"loglik=",logp.hat,"\n")
  acc = 0
  res = c(0,0,alpha_hosp,param_dynamicR$alpha,param_dynamicR$beta,pi_test,logp.hat,acc)
  write(as.vector(round(res,7)),file="res_mcmc_param.txt",9)
  names(res) = c("m","j","alpha_hosp","alpha","beta","pi0","pi1","Log-lik","acc")
  #Choose parameters to sample:
  #Sample which parameter to be changes, if dynamic R, only amp is changes
  for(m in 1:Nmcmc)
   {
      #Sample proposal, only chaching one variable
      j = sample(samp_ind,1)
      param.prop = param
      param.prop[j] = param[j] + sig.prop[j]*rnorm(1)
      alpha_hosp = 4 + exp(param.prop[1])
      param_dynamicR$alpha = exp(param.prop[2])
      param_dynamicR$beta = exp(param.prop[3])
      pi_test[1] = param.prop[4]
      pi_test[2] = param.prop[5]
      pri.prop = f1(alpha_hosp,4)+f2(param_dynamicR$alpha,0.1,1)+f2(param_dynamicR$beta,0.1,1)+
                  f3(pi_test[1],0,5)+f3(pi_test[2],0,0.01)
      start = Sys.time()
      if(0)
       logp.hat.prop = logp.hat + 0.1*rnorm(1)
      else
      logp.hat.prop = smc_loglik_mem(
        mu_R0=mu_R0,sd_R0=sd_R0,mu_AMP=mu_AMP,sd_AMP=sd_AMP,
        dyn_model=dyn_model,
        param_dynamicR=param_dynamicR,
        lambda=lambda,lambda_size=lsize,lambda_ch=lch,
        prob_h=prob_h,change_prob=change_prob,alpha_hosp=alpha_hosp,
        q_test=q_test,alpha_test=alpha_test,pi_test=pi_test,
        start_date=start_date,dynR_date=dynR_date,end_date=end_date,
        hosp_data = hosp_data, test_data = test_data, pop = pop, 
        loc.seed = loc.seed, mob_mat = mob_mat, resdir =resdir, 
        B=B,unit.samp=unit.samp,smo.lag=smo.lag,resamp.frac=resamp.frac,parallel=parallel,
        VERBOSE=VERBOSE,DEBUG=DEBUG,pmcmc=pmcmc,output_latent=FALSE,
        seed=seed+m,N_cores=N_cores)$logp.hat
    if(is.infinite(logp.hat.prop))
        logp.hat.prop = -9999.99        #Somewhat ad hoc
      if(VERBOSE)
        cat("m=",m,"j=",j,"alpha=",alpha_hosp,"alpha=",param_dynamicR$alpha,"beta=",param_dynamicR$beta,"pi=",pi_test,"loglik=",logp.hat.prop,"\n")
      res.prop = c(m,j,alpha_hosp,param_dynamicR$alpha,param_dynamicR$beta,pi_test,logp.hat.prop)
      acc.cur = 0
      if(runif(1)<exp(logp.hat.prop+pri.prop-logp.hat-pri))
        {
	        acc.cur = 1
          if(VERBOSE)
            cat("Accepted\n")
          param = param.prop
          logp.hat = logp.hat.prop
      }
      acc = acc + acc.cur
      res.prop = c(res.prop,acc.cur)
      names(res.prop) = c("m","j","alpha_hosp","alpha","beta","pi0","pi1","Log-lik","acc")
      show(res.prop)
      alpha_hosp = 4 + exp(param[1])
      param_dynamicR$alpha = exp(param[2])
      param_dynamicR$beta = exp(param[3])
      pi_test[1] = param[4]
      pi_test[2] = param[5]
      end = Sys.time()
      res = c(m,j,alpha_hosp,param_dynamicR$alpha,param_dynamicR$beta,pi_test,logp.hat.prop,acc.cur)
      if(VERBOSE)
        cat("Time for iteration=",end-start,"\n")
      write(round(res.prop,7),file="res_mcmc_param.txt",9,append=TRUE)
   }
 res = data.frame(res)
 #if(VERBOSE)
   cat("Acceptance-rate=",acc/Nmcmc,"\n")
 res
}
