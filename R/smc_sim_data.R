#' Marginal loglikelihood estimation by sequential Monte Carlo
#'
#' R is assumed to be constant until 'dynR_date', thereafter following a random walk distribution on log-scale
#' The variance for the random walk is estimated by using the method of sufficient statistics (Storvik-2002).
#' Sequential Monte Carlo is performed through the Bootstrap filter, that is the latent variables (particles) are simulated from their prior.
#' Importance weights are calculated based on hospital data which are assumed to follow a Beta-binomial distribution for hospitalization and
#' a Negative binomial distribution for the length from symptons to hospitalization.
#' 
#' Latent variables are simulated by the prior model (the spread model) using the run_model routine
#' First initial particles are generated until an initial date using the smc_init_sim routine.
#' Within the dynamic period, the particles are simulated using the smc_particle_sim routine.
#' Thereafter simulations are calibrated towards the hospital data.
#' Three possible predictions are made available: Filtering through \eqn{p(x_t|y_{1:t})}, fixed-lag smoothing through \eqn{p(x_t|y_{1:t+u})} and
#' full smoothing through \eqn{p(x_t|y_{1:T})}. All incidence predictions are only given by full smoothing.
#' Simulation results related to R and AMP are written into some output files. All predictions are returned as a list.
#' 
#' Global variables used
#' 'hospD'
#'
#' There are currently many global parameter used in the underlying routine 'run_model'. No control on that
#'
#' @param prob_h Probability for hospitalisation. Default to 0.036.
#' @param alpha_hosp First parameter in beta-binomial distribution for hospitalization. If NULL, the model reduces to a binomial distribution with 'prob_h' being the probability for hospitalization. In general, the second (beta) parameter is calculated as 'alpha_hosp * (1 - prob_h) / prob_h'. Default is NULL.
#' @param alpha_test First parameter in beta-binomial distribution for testdata. If NULL, the model reduces to a binomial distribution with 'prob_h' being the probability for hospitalization. In general, the second (beta) parameter is calculated as 'alpha_hosp * (1 - prob_h) / prob_h'. Default is NULL.
#' @param pi_test Detection probability for test data
#' @param mob_mat list of mobility matrices for each (6 hours) time-interval
#' @param pop data table describing population structure
#' @param lambda Expected length from sickness to hospitalization. Default to 9.08.
#' @param lambda_size Size parameter in the Negative binomial distribution for length until hospitalization. Default is NULL in which case a Poisson distribution is used.
#' @param lambda_ch If 'lambda' and 'lambda_size' are vectors, this correspond to the dates for changes in these parameters
#' @param start_date starting date for initial simulations
#' @param end_date final date for doing simulation (note that hospital data used after this date)
#' @param parallel TRUE of parallelization during computation is to be used. Should always be used (if not in a debugging setting). Default is TRUE.
#' @param VERBOSE If TRUE, some output will be written during excecution. Default is TRUE
#' @param DEBUG If TRUE, extra information is written during excecution. Should only be used during debugging. Default is FALSE.
#' @param output.latent If FALSE, the latent structure is not saved, only information related to the dynamic R and AMP. Can save some memory during testing. Default is FALSE.
#' @param seed Seed number. Default to 123. Note that the 'run_model' function that is called frequently do not use this seed, so that it is currently not possible to obtain a full repetition of the simulations.
#' @param output_latent If FALSE, output of particles for latent variables is suppressed in order to save memory. Default is FALSE
#' @return A list containing particles of latent variables ('latent'), particles of parameters ('param') and estimate of marginal log-likelihood.
#' @export
#' @examples
#' res.smc = smc_loglik(params=c(R0,Reff,AMP),B=10,dynamic.R=TRUE)
smc_sim_data = function(R,
                        mu_AMP = 2.8,sd_AMP = 0.0,
                        lambda = 9.08,lambda_size=NULL, lambda_ch=NULL,
                        prob_h = c(0.036,0.022),change_prob=c(as.Date("2020-05-01")), alpha_hosp = NULL,
                        alpha_test=NULL, q_test=c(0,0.05,0.080,0.16,0.4,0.3,0.01),pi_test=c(-1.175,7.4e-05),
                        start_date = as.Date("2020-02-17"),
                        loc.seed=seeding,mob_mat = mobility_matrices,pop=pop,
                        unit.samp = 3,parallel = TRUE,
                        output.call="call.txt",
                        VERBOSE=TRUE,DEBUG = FALSE,output_latent = FALSE,seed = 123,N_cores=4)
{
  print("Calling smc_sim_data:")
  capture.output(as.list(match.call()),file=output.call)
  #Vector of hospitalisation probabilities
  prob_h<-prob_h[findInterval(seq.Date(from=start_date,to=end_date,by=1),change_prob)+1]
  

  ###############
  ##DELAY OF REPORTING ACCOUNTED AS A CHANGE IN THE HOSPITALISATION PROBABILITY DURING THE LAST 4 DAYS REPORTED
  
  #NOTE THAT THESE ARE IN REVERSE ORDER COMPARED TO EXPAND HOSP!!! (HERE -> DAY-4:LASTDAY)
  delays.sunday<-c(0.86,0.68,0.49,0.32)
  delays.general<-c(0.91,0.82,0.77,0.53)
  period=(end_date-start_date)+1
  if (lubridate::wday(end_date)==1){
    prob_h[(period-3):period]<-prob_h[(period-3):period]*delays.sunday
  } else {
    prob_h[(period-3):period]<-prob_h[(period-3):period]*delays.general
  }
  #################
  
  
  
  #Number of days
  ndays = length(R)
  days <- as.Date(start_date, origin = "1970-01-01") + lubridate::days(1:ndays)
  end_date = tail(days,1)
  nloc = nrow(pop)
  set.seed(seed)
  
  
  #Calculation of probabilities for length until hospitalisations, we need the hasard
  #P(Z=k|Z>k-1), either from a Poisson or a Negative binomial distribution
  #Define maximum number of days to hospitalization
  
  lambda_cur = lambda[1]
  if(length(lambda_ch)!=(length(lambda)-1))
  {
    show(lambda)
    show(lambda_ch)
    print("smc_loglik_mem: Missmatch between length of lambda and length of lambda_ch")
    return(NULL)
  }
  if(is.null(lambda_ch))
    lambda_ch = NA
  if(!is.null(lambda_size))
  {
    lambda_size_cur = lambda_size[1]
    if(length(lambda_size)!=length(lambda))
    {
      show(lambda)
      show(lambda_size)
      print("smc_loglik_mem: Missmatch between length of lambda and length of lambda_size")
      return(NULL)
    }
  }
  
  pbar = smc_calc_lambda(lambda_cur,lambda_size_cur)
  MaxdaysH = length(pbar)-1
  
  #Calculation of probabilities for length until testing, we need the hasard
  #P(Z=k|Z>k-1)
  #Define maximum number of days to hospitalization
  MaxdaysT = length(q_test)
  qbar = q_test
  for(i in 2:length(q_test))
    qbar[i] = qbar[i]/sum(qbar[i:length(qbar)])
  
  AMP_sim = mu_AMP * exp(rnorm(1, -0.5 * sd_AMP*sd_AMP, sd_AMP))
  #Then start simulations one-day ahead calculating importance weights
  
  if(0)
  {
  inshist = matrix(0,nrow=MaxdaysH,ncol=MaxdaysH)
  instest = matrix(0,nrow=MaxdaysT,ncol=MaxdaysT)
  
  smc.hospD = as.data.table(data.frame(date=days,N=NA))
  smc.testD = as.data.table(data.frame(date=days,Negative=NA,Positive=NA,PositiveIkkeUtlandet=NA))
    input.particle = pop
  input = list(particle=pop,R_sim=R,AMP_sim=AMP_sim,
               inshist=inshist,instest=instest)
  res = smc_particle_sim(input,input_date=start_date,unit.samp=ndays,day.numb=1,
                             loc.seed=loc.seed,mob_mat=mob_mat,hospD=smc.hospD,testD=smc.testD,pbar=pbar,qbar=qbar,prob_h=prob_h,
                             alpha_hosp=alpha_hosp,alpha_test=alpha_test,pi_test=pi_test,MaxdaysH=MaxdaysH,MaxdaysT=MaxdaysT)
  }

  particle <- smc_run_model(mob_mat, 
                     start_date = start_date, 
                     end_date = end_date,
                     N_sim = 1,
                     R_0 = R[1],
                     reff_change_dates=days[-1],
                     Reffs = R[-1],
                     parallel = TRUE,
                     seeding=loc.seed,
                     se1e2iiar_pop=pop,
                     AMP_factor = AMP_sim)[[1]]
  insH0 = particle[,.(ins=sum(c_symp_incidence+c_asymp_incidence)),by=day]$ins
  insT0 = insH0
  N = length(insH0)
  insH = matrix(nrow=N,ncol=MaxdaysH)
  insT = matrix(nrow=N,ncol=MaxdaysT)
  for(i in 1:N)
  {
    #for(j in 1:min(MaxdaysH,N-i+1))
    for(j in 1:MaxdaysH)
    {
      insH[i,j] = rbinom(1,insH0[i],pbar[j])
      insH0[i] = insH0[i]-insH[i,j]
    }
    for(j in 1:MaxdaysT)
    {
      insT[i,j] = rbinom(1,insT0[i],qbar[j])
      insT0[i] = insT0[i]-insT[i,j]
    }
  }
  insH = cbind(insH0,insH)
  insT = cbind(insT0,insT)
  
  N = ndays
  for(i in 1:ndays)
  {
    NN = min(MaxdaysH,i)
    ind = 1+c(1:NN)*(N-1)  #Index to pick out right values of insH
    #num[i] = sum(insH[M+unit.samp+1-i+ind],na.rm=TRUE)
    num[i] = sum(insH[M+i+ind],na.rm=TRUE)
    if(is.null(alpha_hosp))
          hospD$N[i] = rbinom(1,num[i],prob_h[M+i])
        else
          hospD$N[i] = rbbinom(1,num[i],alpha=alpha_hosp,beta=beta_hosp[M+i])
       numT[i] = sum(insT[M+i+ind],na.rm=TRUE)
        eta.t = pi_test[1]+pi_test[2]*testD$N[M+i]
        pi.t = exp(eta.t)/(1+exp(eta.t))
        #print("sim_particle_sim")
        #show(c(M+i,numT[i]))
        #show(testD[M+i,])
        if(testD$Positive[M+i]>numT[i])
        {
          #l =  -9999                #Somewhat ad hoc
          #logw.new = logw.new
          l = -9999
          #print("Warning: number of tests larger than number of incidences")
          #show(c(M+i,numT[i],pi.t))
          #show(as.matrix(testD[day.numb+1-1,]))
        }
        else
        {
          if(is.null(alpha_hosp))
            l = dbinom(testD$Positive[M+i],numT[i],pi.t,log=TRUE)
          else
            l = dbbinom(testD$Positive[M+i],numT[i],
                        alpha=alpha_test,beta=alpha_test*(1-pi.t)/pi.t,log=TRUE)
          if(is.na(l))
            show(c(l))
        }
        logw.new = logw.new+l
        if(is.na(logw.new))
          print("smc_particle_sim: Something is wrong")
        ll[i] = l
      }
  return(list(particle=particle,insH=insH,insT=insT))
}
