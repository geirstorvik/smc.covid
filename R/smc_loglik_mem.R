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
#' There are currently many global parameter used in the underlying routine 'run_model'. No control on that
#'
#' @param mu_R0 R0 in the first static period is assumed to follow a lognormal distribution with expectation mu_R0. Default is 3.3.
#' @param sd_R0 The standard deviation in the lognormal distribution for R0. Default is 0.1.
#' @param dyn_model The model for dynamic R, currently only RW (random walk) is implemented
#' @param param_dynamicR List with parameters related to dynamic model for R. alpha,beta are shape and rate parameters for gamma distribution of precision, a,b are parameters in beta distribution for prob of changepoint, mu0,kappa0 are mean and proportional precision for mean
#' @param mu_AMP R0 The amplification factor AMP is assumed to follow a lognormal distribution with expectation mu_AMP. Default is 2.8.
#' @param sd_AMP The standard deviation in the lognormal distribution for AMP. Default is 0.0 in which case no updating of AMP is performed.
#' @param prob_h Probability for hospitalisation. Default to 0.036.
#' @param alpha_hosp First parameter in beta-binomial distribution for hospitalization. If NULL, the model reduces to a binomial distribution with 'prob_h' being the probability for hospitalization. In general, the second (beta) parameter is calculated as 'alpha_hosp * (1 - prob_h) / prob_h'. Default is NULL.
#' @param alpha_test First parameter in beta-binomial distribution for testdata. If NULL, the model reduces to a binomial distribution with 'prob_h' being the probability for hospitalization. In general, the second (beta) parameter is calculated as 'alpha_hosp * (1 - prob_h) / prob_h'. Default is NULL.
#' @param pi_test Detection probability for test data
#' @param lambda Expected length from sickness to hospitalization. Default to 9.08.
#' @param lambda_size Size parameter in the Negative binomial distribution for length until hospitalization. Default is NULL in which case a Poisson distribution is used.
#' @param lambda_ch If 'lambda' and 'lambda_size' are vectors, this correspond to the dates for changes in these parameters
#' @param start_date starting date for initial simulations
#' @param init_date date until particles are simulated before calibration towards data
#' @param end_date final date for doing simulation (note that hospital data used after this date)
#' @param hosp_data data frame containing hospital data, with components date and N.
#' @param test_data data frame containing test data, with components ...If NULL, only hospital data is used. Can contain NA
#' @param mob_mat list of mobility matrices for each (6 hours) time-interval
#' @param pop data table describing population structure
#' @param output.R Filename for writing out dynamic R simulations
#' @param output.R_smo Filename for writing out dynamic R simulations smoothing backwards.
#' @param output.R_lag Filename for writing out dynamic R simulations smoothing 'smo.lag' units backwards.
#' @param output.C Filename for writing out c_symp+c_asymp
#' @param output.C_lag Lagged version of previous
#' @param output.C_smo Smoothed version of previous
#' @param output.phat Writing estimate of marginal likelihood to file
#' @param B Number of particles. Default to 20.
#' @param B2 The storage of all the latent variable could be problematic due to memory. A smaller number of particles to be stored can be specified through B2 in which case only B2 particles are sampled at the end and returned. Default is NULL is which case no downscaling is performed.
#' @param unit.samp The number of days for which simulation is performed within each iteration (in the dynamic period). When call to the run_model routine is made within each iteration, some of the latent variables (mobility information) is not transfered, making the simulations somewhat approiximate. A large unit.samp value is preferred from this perspective. However, importance weights are calculated based on the number of days that are simulated, giving more extreme weights for a larger number of unit.samp. A trade-off is therefore needed. Values of 3 and 5 seems to work ok. Default is 3.
#' @param resamp.frac Resampling will increase variability at the current time (but might decrease later). When the effective sample size is smaller than B*resamp.frac, resampling is performed. For the current implementation a value larger than 1 is recommended, in which case resampling will always be performed. Default is 1.1.
#' @param smo.lag Fixed lag smoothing is based on inference from \eqn{p(x_t|y_{t+h})} where $h$ corresponds to 'smo.lag'. A large value increase the data information but due to resampling may leed to very few unique samples. Recommended (and default) value is 25.
#' @param parallel TRUE of parallelization during computation is to be used. Should always be used (if not in a debugging setting). Default is TRUE.
#' @param VERBOSE If TRUE, some output will be written during excecution. Default is TRUE
#' @param DEBUG If TRUE, extra information is written during excecution. Should only be used during debugging. Default is FALSE.
#' @param output_latent If FALSE, the latent structure is not saved, only information related to the dynamic R and AMP. Can save some memory during testing. Default is FALSE.
#' @param seed Seed number. Default to 123. Note that the 'run_model' function that is called frequently do not use this seed, so that it is currently not possible to obtain a full repetition of the simulations.
#' @param output_latent If FALSE, output of particles for latent variables is suppressed in order to save memory. Default is FALSE
#' @param fix.R If non-null, contains an ndays x B matrix with fixed R's (so they are not simulated)
#' @param sim.data If TRUE, data are simulated, and replaces the input data
#' @return A list containing particles of latent variables ('latent'), particles of parameters ('param') and estimate of marginal log-likelihood.
#' @export
#' @examples
#' res.smc = smc_loglik(params=c(R0,Reff,AMP),B=10,dynamic.R=TRUE)
smc_loglik_mem = function(mu_R0 = 3.3,sd_R0 = 0.1,mu_AMP = 2.8,sd_AMP = 0.0,
                          dyn_model="RW",
                          param_dynamicR = list(alpha=2.4,beta=0.28,phi.a=100,phi.b=900,a0=0.8,mu0=0,kappa0=1),
                          lambda = 9.08,lambda_size=NULL, lambda_ch=NULL,
                          prob_h = c(0.036,0.022),change_prob=c(as.Date("2020-05-01")), alpha_hosp = NULL,
                          alpha_test=NULL, q_test=c(0,0.05,0.080,0.16,0.4,0.3,0.01),pi_test=c(-1.175,7.4e-05),
                          start_date = as.Date("2020-02-17"),dynR_date = as.Date("2020-03-08"),end_date = as.Date("2020-05-05"),
                          hosp_data = hospitalisations_df,test_data=NULL,pop=pop,
                          loc.seed=seeding,mob_mat = mobility_matrices,
                          resdir = "./",
                          output.call="call.txt",
                          output.R = "smc_dynamic_R.txt",output.R_smo = "smc_dynamic_R_smo.txt",output.R_lag = "smc_dynamic_R_lag.txt",
                          output.W = "smc_weights.txt",output.AMP = "smc_AMP.txt",output.SdR = "smc_SdR.txt",
                          output.C = "smc_dynamic_C.txt",output.C_smo = "smc_dynamic_C_smo.txt",output.C_lag = "smc_dynamic_C_lag.txt",
                          output.phat = "smc_phat.txt",
                          B = 2000,Bmem=B,B2 = NULL,unit.samp = 3,resamp.frac = 1.1,smo.lag = 25,parallel = TRUE,pmcmc=FALSE,
                          fix.R=NULL,sim.data=FALSE,
                          VERBOSE=TRUE,DEBUG = FALSE,output_latent = FALSE,seed = 123,N_cores=4)
{
 print("Calling smc_loglik_mem:")
 output.call = paste(resdir,output.call,sep="")
 capture.output(as.list(match.call()),file=output.call)
 
 output.R = paste(resdir,output.R,sep="")
 file.create(output.R)
 output.R_smo = paste(resdir,output.R_smo,sep="")
 file.create(output.R_smo)
 output.R_lag = paste(resdir,output.R_lag,sep="")
 file.create(output.R_lag)
 output.W = paste(resdir,output.W,sep="")
 file.create(output.W)
 output.AMP = paste(resdir,output.AMP,sep="")
 file.create(output.AMP)
 output.SdR = paste(resdir,output.SdR,sep="")
 file.create(output.SdR)
 output.C = paste(resdir,output.C,sep="")
 file.create(output.C)
 output.C_smo = paste(resdir,output.C_smo,sep="")
 file.create(output.C_smo)
 output.C_lag = paste(resdir,output.C_lag,sep="")
 file.create(output.C_lag)
 output.phat = paste(resdir,output.phat,sep="")
 
 #Vector of hospitalisation probabilities
 prob_h<-prob_h[findInterval(seq.Date(from=start_date,to=end_date,by=1),change_prob)+1]
 
 if(!pmcmc)
 {
  smc_init_sim_sel = smc_init_sim
  smc_particle_sim_sel = smc_particle_sim
 }
 else
 {
  #smc_init_sim_sel = smc_init_sim_pmcmc
  #smc_particle_sim_sel = smc_particle_sim_pmcmc
  smc_init_sim_sel = smc_init_sim
  smc_particle_sim_sel = smc_particle_sim
  
 }
 
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
 
 if (unit.samp > smo.lag)
  cat("Warning: smo.lag should be larger than unit.samp\n")
 
 #Transform data to suiteable form
 smc.hospD = hosp_data[(date >= start_date) & (date <= end_date)]
 if(smc.hospD[1,]$date > start_date)
 {
  #Fill inn 
  l = smc.hospD[1,]$date-start_date
  dd = start_date+lubridate::days(1:l)-1
  smc.hospD=rbind(data.table(data.frame(date=dd,N=NA)),smc.hospD)
 }
 n = ncol(smc.hospD)
 if(smc.hospD[n,]$date < start_date)
 {
  #Fill inn 
  l = end_date-smc.hospD[1,]$date
  dd = end_date+lubridate::days(1:l)-1
  smc.hospD=rbind(smc.hospD,data.table(data.frame(date=dd,N=NA)))
 }
 if(!is.null(test_data))
  smc.testD = smc_convert_test_data(test_data,smc.hospD$date,pi_test)
 else
  smc.testD = NULL

 #Number of days
 ndays = as.numeric(end_date - start_date)+1
 ndays.init = as.numeric(dynR_date - start_date)+1
 ndays.smc = as.numeric(end_date - dynR_date)
 numb.smc = ceiling(ndays.smc / unit.samp)
 nloc = nrow(pop)
 unit.samp.save = unit.samp       #unit.samp changed for last interval
 set.seed(seed)

  ###############
 ##Preparation for dynamic model for R
 cp_probD = NULL
 if(dyn_model=="RW")
 {
  Nsuff = 2
  npardyn = 1
  smc_sim_R = smc_sim_R_RW
  #Transform mean and variance to shape and scale
  param_R = rep(NA,2)
  #param_R[1] = 2 + param_dynamicR$mu * param_dynamicR$mu / param_dynamicR$V
  #param_R[2] = (param_R[1] - 1) * param_dynamicR$mu
  param_R[1] = param_dynamicR$alpha
  param_R[2] = param_dynamicR$beta
 }
 else if(dyn_model=="AR")
 {
  Nsuff = 6
  npardyn = 3
  smc_sim_R = smc_sim_R_AR
  #Transform mean and variance to shape and scale
  param_R = rep(NA,4)
  param_R[1] = param_dynamicR$mu0
  param_R[2] = param_dynamicR$kappa0
  param_R[3] = param_dynamicR$alpha
  param_R[4] = param_dynamicR$beta
 }
 else if(dyn_model=="AR_fix_par")
 {
   Nsuff = 6
   npardyn = 3
   smc_sim_R = smc_sim_R_AR_fix_par
   #Transform mean and variance to shape and scale
   param_R = rep(NA,4)
   param_R[1] = param_dynamicR$a0
   param_R[2] = param_dynamicR$kappa0
   param_R[3] = param_dynamicR$alpha
   param_R[4] = param_dynamicR$beta
 }
 else if(dyn_model=="AR_mu_0")
 {
   Nsuff = 6
   npardyn = 3
   smc_sim_R = smc_sim_R_AR_mu_0
   #Transform mean and variance to shape and scale
   param_R = rep(NA,4)
   param_R[1] = param_dynamicR$a0
   param_R[2] = param_dynamicR$kappa0
   param_R[3] = param_dynamicR$alpha
   param_R[4] = param_dynamicR$beta
 }
 else if(dyn_model=="AR_a8")
 {
   Nsuff = 6
   npardyn = 3
   smc_sim_R = smc_sim_R_AR_a8
   #Transform mean and variance to shape and scale
   param_R = rep(NA,4)
   param_R[1] = param_dynamicR$mu0
   param_R[2] = param_dynamicR$kappa0
   param_R[3] = param_dynamicR$alpha
   param_R[4] = param_dynamicR$beta
 }
 else if(dyn_model=="AR_fix_a")
 {
   Nsuff = 6
   npardyn = 3
   smc_sim_R = smc_sim_R_AR_fix_a
   #Transform mean and variance to shape and scale
   param_R = rep(NA,4)
   param_R[1] = param_dynamicR$mu0
   param_R[2] = param_dynamicR$kappa0
   param_R[3] = param_dynamicR$alpha
   param_R[4] = param_dynamicR$beta
 }
 else if(dyn_model=="CP")
 {
  Nsuff = 3
  npardyn = 2
  smc_sim_R = smc_sim_R_CP
  #Transform mean and variance to shape and scale
  param_R = rep(NA,4)
  param_R[1] = param_dynamicR$phi.a
  param_R[2] = param_dynamicR$phi.b
  param_R[3] = param_dynamicR$alpha
  param_R[4] = param_dynamicR$beta
  cp_probD = rep(NA,ndays)
 }
 else if(dyn_model=="TCP")
 {
  Nsuff = 3
  npardyn = 2
  smc_sim_R = smc_sim_R_TCP
  #Transform mean and variance to shape and scale
  param_R = rep(NA,4)
  param_R[1] = param_dynamicR$phi.a
  param_R[2] = param_dynamicR$phi.b
  param_R[3] = param_dynamicR$alpha
  param_R[4] = param_dynamicR$beta
 }
 else
 {
  cat("Unknown model for dynamic R:",dyn_model,"\n")
 }
 
 
 
 
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
 
 
 logw = rep(NA, B)
 #First run model from start_date to dynR_date
 
 #Simulating static parameters
 R_sim = matrix(nrow = unit.samp + 1, ncol = B)
 if(is.null(fix.R))
  R_sim[1, ] = mu_R0 * exp(rnorm(B, -0.5 * sd_R0*sd_R0, sd_R0))
 else
  R_sim[1,] = fix.R[1,]
 if(DEBUG)
 {
  cat("R_sim for initial period:",R_sim[1,],"\n")
 }
 start = Sys.time()
 AMP_sim = mu_AMP * exp(rnorm(B, -0.5 * sd_AMP*sd_AMP, sd_AMP))
 if (parallel)
 {
  #Simulate initial days with fixed R0
  input = list()
  for (b in 1:B)
   input[[b]] = list(R0 = R_sim[1, b], AMP = AMP_sim[b])
  res  <- parallel::mclapply(input,smc_init_sim_sel,start_date=start_date,end_date=dynR_date,
                             mob_mat=mob_mat[1:(4*ndays.init)],loc.seed=loc.seed,pop=pop,MaxdaysH=MaxdaysH,pbar=pbar,
                             MaxdaysT=MaxdaysT,qbar=qbar,
                             mc.cores = N_cores)
  particle <- sapply(res, '[', 'particle')
  inshist = sapply(res, '[', 'insH')
  instest = sapply(res, '[', 'insT')
 }
 else
 {
  particle = list()
  inshist = list()
  instest = list()
  for (b in 1:B)
  {
   #INITIALIZATION PREVIOUS TO DYNAMIC R!! ALF
   
   input = list(b = b,R0 = R_sim[1, b],AMP = AMP_sim[b])
   res=smc_init_sim_sel(input,start_date=start_date,end_date=dynR_date,
                        mob_mat=mob_mat[1:(4*ndays.init)],loc.seed=loc.seed,pop=pop,MaxdaysH=MaxdaysH,pbar=pbar,MaxdaysT=MaxdaysT,qbar=qbar)
   particle[[b]] = res$particle
   inshist[[b]] = res$insH
   instest[[b]] = res$insT
  }
 }
 end = Sys.time()
 if (VERBOSE)
 {
   cat("Time for running initial particles\n")
   show(end - start)
 }
 #Put all weights equal
 w = rep(1 / B, B)
 logw = log(w)
 #Then start simulations one-day ahead calculating importance weights
 logp.hat = 0
 input_date = dynR_date 
 
 #Store results from initial simulations
 write(t(matrix(rep(R_sim[1, ],ndays.init),nrow=ndays.init,byrow=TRUE)),output.R,ncolumns=B,append=TRUE)
 #Wmatrix = matrix(NA, nrow = ndays, ncol = B)
 #Wmatrix[1:ndays.init, ] = 1
 write(t(matrix(1,nrow=ndays.init,ncol=B)),output.W,ncolumns=B,append=TRUE)
 Rfilt = matrix(NA, nrow = unit.samp, ncol = B)
 Rmatrix_smo = matrix(NA, nrow = ndays, ncol = B)
 #Rmatrix_lag = matrix(NA, nrow = ndays, ncol = B)
 #Rmatrix[1:ndays.init, ] = matrix(rep(R_sim[1, ], ndays.init), nrow = ndays.init, byrow =
 #                                 TRUE)
 Rmatrix_smo[1:ndays.init, ] = matrix(rep(R_sim[1, ], ndays.init), nrow =
                                       ndays.init, byrow = TRUE)
 if(ndays.init>smo.lag)
 {
   write(t(Rmatrix_smo[1:(ndays.init-smo.lag),,drop=FALSE]),output.R_lag,ncolumns=R,append=TRUE)
 }
 #Rmatrix_lag[1:ndays.init, ] = matrix(rep(R_sim[1, ], ndays.init), nrow =
 #                                     ndays.init, byrow = TRUE)
 #AMPmatrix = matrix(NA, nrow = ndays, ncol = B)
 #AMPmatrix[1:ndays.init, ] = matrix(rep(AMP_sim, ndays.init), nrow = ndays.init, byrow =TRUE)
 write(t(matrix(rep(AMP_sim, ndays.init), nrow = ndays.init, byrow =TRUE)),output.AMP,ncolumns=B,append=TRUE)
 #SdRmatrix = array(NA, c(ndays,npardyn,B))
 for(i in 1:npardyn)
 {
   file.create(paste(output.SdR,i,sep=""))
   for(u in 1:ndays.init)
    write(rep(NA,B), file = paste(output.SdR,i,sep=""),ncolumns=B,append=TRUE)
 }
 
 #write(t(matrix(NA,nrow=ndays.init,ncol=B)),output.SdR,ncolumns=B)
 
 #Cmatrix = matrix(NA, nrow = ndays, ncol = B)
 Cmatrix_smo = matrix(NA, nrow = ndays, ncol = B)
 #Cmatrix_lag = matrix(NA, nrow = ndays, ncol = B)
 #Cmatrix[1:ndays.init, ] = smc_extract_c_list(particle)
 #Cmatrix_lag[1:ndays.init, ] = smc_extract_c_list(particle)
 Cmatrix_smo[1:ndays.init, ] = smc_extract_c_list(particle)
 write(t(Cmatrix_smo[1:ndays.init, ]),output.C,ncolumns=B,append=TRUE)

 hospD.sim = NULL
 testD.sim = NULL
 if(sim.data)
 {
  hospD.sim = matrix(nrow=ndays,ncol=B)
  testD.sim = matrix(nrow=ndays,ncol=B)
 }
 
 #Remove initial days from mobility matrix
 #mob_mat = mob_mat[-c(1:(4 * ndays.init))]
 loc.seed$day = loc.seed$day - ndays.init
 loc.seed = loc.seed[day >= 0]
  
 #Initialise for simulations with dynamic R
 ind = 1:B
 #count_days = 0
 Suff_R = matrix(0,nrow=B,ncol=Nsuff)
 ref_change_dates = c()
 if(VERBOSE)
  cat("init_days=", ndays.init, "\n")
 Mdays = list()
 ltest = matrix(nrow=ndays,ncol=Bmem)
 llogw = matrix(nrow=ndays,ncol=Bmem)
 particle.new = particle
 for (d in 1:numb.smc)
 {
  start = Sys.time()
  #Correct for length of last interval
  if ((input_date + unit.samp) > end_date)
  {
   unit.samp = as.numeric(end_date - input_date)
  }
  day.numb = as.numeric(input_date - start_date) + 1
  days = seq(day.numb+1, day.numb + unit.samp, by = 1)
  Mdays[[d]] = days
  if(VERBOSE)
   cat("Processing days=", days, "\n")
  gc() 
  #Check if change date for hospital distr is within these dates
  if(!is.null(lambda_ch) & !is.na(pmatch(lambda_ch,input_date+unit.samp)))
  {
   print("smc_loglik_mem: Changing lambda")
   show(lambda_ch)
   show(input_date)
   show(unit.samp)
   lambda_cur = lambda[2]
   lambda_size_cur = lambda_size[2]
   pbar = smc_calc_lambda(lambda_cur,lambda_size_cur,MaxdaysH)
   show(pbar)
  }
  
  #Sample indices to be used for proposals
  ind.mem = smc_resamp(w,Bmem)
  
  #Simulate R's for next unit.samp days, note first simulate par.R
  Suff_R = Suff_R[ind.mem,,drop=FALSE]
  R_sim = R_sim[,ind.mem]
  AMP_sim = AMP_sim[ind.mem]
  res = smc_sim_R(R_sim,Suff_R,unit.samp,param_R)
  R_sim[1+1:unit.samp,] = res$R_sim[1+1:unit.samp,]
  Suff_R = res$Suff_R
  par.R = res$par.R
  if(!is.null(fix.R))
  {
   R_sim[1+1:unit.samp,] = fix.R[day.numb+1:unit.samp,]
  }
  #Particle simulations unit.samp days ahead
  last.day = max(particle[[1]]$day)
  if (parallel)
  {
   input = list()
  for (b in 1:Bmem)
   {
    b2 = ind.mem[b]
    input.particle = particle[[b2]][day==last.day,.(location_code,b_S,b_E1,b_E2,b_I,b_Ia,b_R)]
    names(input.particle) = names(pop)
    input[[b]] = list(particle=input.particle,R_sim=R_sim[-1,b],AMP_sim=AMP_sim[b],
                      inshist=inshist[[b2]],instest=instest[[b2]])
   }
   res = parallel::mclapply(input,smc_particle_sim_sel,mc.cores = N_cores,
                            input_date=input_date,unit.samp=unit.samp,day.numb=day.numb,loc.seed=loc.seed,
                            mob_mat=mob_mat[4*(day.numb-1)+1:(4*unit.samp)],
                            hospD=smc.hospD[day.numb-1+1:unit.samp,],testD=smc.testD[day.numb-1+1:unit.samp,],
                            pbar=pbar,qbar=qbar,prob_h=prob_h[day.numb-1+1:unit.samp],
                            alpha_hosp=alpha_hosp,alpha_test=alpha_test,pi_test=pi_test,
                            MaxdaysH=MaxdaysH,MaxdaysT=MaxdaysT,sim.data=sim.data)
   particle.new <- sapply(res, '[', 'particle.new')
   inshist = sapply(res, '[', 'insH')
   instest = sapply(res, '[', 'insT')
   numH = sapply(res,'[','numH')
   numT = sapply(res,'[','numT')
   logw.new = unlist(sapply(res, '[', 'logw.new')) + logw[ind.mem]
   ltest[days,] = unlist(sapply(res, '[', 'ltest')) 
   llogw[days,] = logw.new
   if(sim.data)
   {
     D.sim = sapply(res,'[', 'D.sim')
     for(b in 1:B)
     {
       hospD.sim[day.numb+1:unit.samp,b] = D.sim[[b]][,1]
       testD.sim[day.numb+1:unit.samp,b] = D.sim[[b]][,2]
     }
   }
  }
  else
  {
   particle.new = list()
   logw.new = rep(NA, Bmem)
   numH = list()
   numT = list()
   inshist.save = inshist
   instest.save = instest
   for (b in 1:Bmem)
   {
    b2 = ind.mem[b]
    input.particle = particle[[b2]][day==last.day,.(location_code,b_S,b_E1,b_E2,b_I,b_Ia,b_R)] # Stores result from the previous day
    names(input.particle) = names(pop)
    input = list(particle=input.particle,R_sim=R_sim[-1,b],AMP_sim=AMP_sim[b],inshist=inshist.save[[b2]],instest=instest.save[[b2]])
    res = smc_particle_sim_sel(input,input_date=input_date,unit.samp=unit.samp,day.numb=day.numb,
                               loc.seed=loc.seed,mob_mat=mob_mat[4*(day.numb-1)+1:(4*unit.samp)],
                               hospD=smc.hospD[day.numb-1+1:unit.samp,],testD=smc.testD[day.numb-1+1:unit.samp,],
                               pbar=pbar,qbar=qbar,prob_h=prob_h[day.numb-1+1:unit.samp],
                               alpha_hosp=alpha_hosp,alpha_test=alpha_test,pi_test=pi_test,MaxdaysH=MaxdaysH,MaxdaysT=MaxdaysT,
                               sim.data=sim.data)
    particle.new[[b]] = res$particle.new
    inshist[[b]] = res$insH
    instest[[b]] = res$insT
    numH[[b]] = res$numH
    numT[[b]] = res$numT
    logw.new[b] = res$logw.new+logw[b2]
    ltest[days,b] = res$ltest
    llogw[days,b] = logw.new[b]
    if(sim.data)
    {
     hospD.sim[day.numb+1:unit.samp,b2] = res$D.sim[,1]
     testD.sim[day.numb+1:unit.samp,b2] = res$D.sim[,2]
    }
   }
  }
  #show(sort( sapply(ls(),function(x){object.size(get(x))})) )
  #show(lapply(res[[1]],object.size))
  logw = logw.new
  #Marginal log-likelihood estimate based on unnormalized weights
  logp.hat = logp.hat + log(mean(exp(logw)))
  #Normalized weights to be used on next time point
  w = exp(logw - max(logw))
  w = w / sum(w)
  #if(B<15)
  #{
  #  show("w")
  #  show(w)
  #}
  
  if(dyn_model=="CP")
  {
    cp_prob = matrix(NA,nrow=unit.samp,ncol=B)
    for(u in 1:unit.samp)
      cp_prob[u,] = (R_sim[u+1,]-R_sim[u,])>0
  }
  
  #Changing days in particle.new (could be done more efficiently)
  for(b in 1:Bmem)
   particle.new[[b]]$day = particle.new[[b]]$day + day.numb - 1
  
  ind = 1:B
  #Resampling if Bmem>B or Beff<resamp.frac*B,
  if(sum(is.na(w))>0)
  {
   print("smc_particle_sim: Missing values in w")
   show(w)
  }
  Beff = 1 / sum(w * w)
  if ((Bmem>B) | (Beff < (resamp.frac * B)))
  {
   ind = smc_resamp(w, B)
   if(VERBOSE)
    cat("Number of unique particles=",length(unique(ind)),"\n")
   R_sim = R_sim[, ind]
   AMP_sim = AMP_sim[ind]
   Rmatrix_smo = Rmatrix_smo[, ind.mem[ind]]      #Resampling to get smoothed estimate
   Suff_R = Suff_R[ind,,drop=FALSE]
   inshist = inshist[ind]
   instest = instest[ind]
   par.R = par.R[ind,,drop=FALSE]
   if(dyn_model=="CP")
     cp_prob = cp_prob[,ind]
   w = rep(1, B)
   logw = log(w)
  }
  if(dyn_model=="CP")
    cp_probD[day.numb+1:unit.samp] = rowMeans(cp_prob)
  
  
  #Store filtered estimates 
  #Wmatrix[days,]=matrix(rep(w,unit.samp),ncol=unit.samp,byrow=T)
  write(t(R_sim[1+1:unit.samp,]),file=output.R,ncolumns=B,append=TRUE)
  write(t(matrix(rep(w,unit.samp),ncol=unit.samp,byrow=T)),file=output.W,ncolumns=B,append=TRUE)
  #Rmatrix[days,]=R_sim[1+1:unit.samp,]    
  
  #Combine new set of particles with previous
  if(VERBOSE)
   cat("Resampling and combining particles\n")
  particle.new = particle.new[ind]
  #particle = particle[ind.mem[ind]]
  #for (b in 1:B)
  # particle[[b]] = rbind(particle[[b]], particle.new[[b]])
  particle = particle.new
  
  #Store results
  days.lag = days-smo.lag
  days.lag = days.lag[days.lag>0]
  if(VERBOSE)
   cat("Store results\n")
  if(DEBUG)
  {
   show(days)
   show(days.lag)
  }
  days = day.numb + 1:unit.samp
  Rmatrix_smo[days,]=R_sim[1+1:unit.samp,]
  #AMPmatrix[days,]=matrix(rep(AMP_sim,unit.samp),ncol=unit.samp,byrow=T)
  write(t(matrix(rep(AMP_sim,unit.samp),ncol=unit.samp,byrow=T)),file=output.AMP,ncolumns=B,append=TRUE)
  #for(u in 1:unit.samp)
  # SdRmatrix[days[u],,]=t(par.R)
  for(i in 1:npardyn)
    for(u in 1:unit.samp)
      write(par.R[,i], file = paste(output.SdR,i,sep=""),ncolumns=B,append=TRUE)
  #write("Parameter a")
  #show(range(par.R[,1]))
  
  #Rmatrix_lag[days,]=Rmatrix_smo[days,]
  Cmatrix_smo[days,] = smc_extract_c_list(particle.new)
  write(t(Cmatrix_smo[days,]),output.C,ncolumns=B,append=TRUE)
  if(days[unit.samp]>smo.lag)
  {
    days.lag = days-smo.lag
    days.lag = days.lag[days.lag>0]
    write(t(Rmatrix_smo[days.lag,,drop=FALSE]),output.R_lag,ncolumns=B,append=TRUE)
    write(t(Cmatrix_smo[days.lag,,drop=FALSE]),output.C_lag,ncolumns=B,append=TRUE)
  }
  #Cmatrix[days,] = smc_extract_c_list(particle.new)
  #Cmatrix_lag[days,] = Cmatrix_smo[days,]
  #if(length(days.lag)>0)
  #{
   #Rmatrix_lag[days.lag,]=Rmatrix_smo[days.lag,]  #Overwrite earlier to get fixed lag smoother
   #Cmatrix_lag[days.lag,]=Cmatrix_smo[days.lag,]  #Overwrite earlier to get fixed lag smoother
  #}
  #Change to new time-inv
  input_date = input_date + unit.samp
  R_sim[1, ] = R_sim[unit.samp + 1, ]              #Last from current period
  #count_days = count_days + unit.samp
  loc.seed$day = loc.seed$day - unit.samp
  loc.seed = loc.seed[day >= 0]
  #mob_mat = mob_mat[-c(1:(4 * unit.samp))]

  end = Sys.time()
  if (VERBOSE)
  {
    cat("Time for ", unit.samp, "days particle simulation\n")
    show(end - start)
  }
 }
 
 
 write(t(Rmatrix_smo[(ndays-smo.lag+1):ndays,,drop=FALSE]),output.R_lag,ncolumns=B,append=TRUE)
 write(t(Cmatrix_smo[(ndays-smo.lag+1):ndays,,drop=FALSE]),output.C_lag,ncolumns=B,append=TRUE)
 if(VERBOSE)
  cat("Convert to list\n")
 #for (d in (ndays - smo.lag):ndays)
 #  Rmatrix_lag[d, ] = Rmatrix_smo[d, ]
 if(VERBOSE)
  cat("Convert to list2\n")
 #Convert results to list
 #res.param = list()
 #for (b in 1:B)
 #res.param[[b]]=data.table(w=Wmatrix[, b],R=Rmatrix[,b],R_smo=Rmatrix_smo[,b],
 #                            R_lag=Rmatrix_lag[,b],AMP=AMPmatrix[,b],SdR=SdRmatrix[,,b])
 
 if(!is.null(output.R))
 {
  if(VERBOSE)
   cat("Write results to file\n")
  
  #write.table(Wmatrix, file = output.W)
  #write.table(Rmatrix, file = output.R)
  write.table(Rmatrix_smo, file = output.R_smo)
  #write.table(Rmatrix_lag, file = output.R_lag)
  #write.table(AMPmatrix, file = output.AMP)
  #for(i in 1:npardyn)
  # write.table(SdRmatrix[,i,], file = paste(output.SdR,i,sep=""))
  #write.table(Cmatrix, file = output.C)
  #write.table(Cmatrix_lag, file = output.C_lag)
  write.table(Cmatrix_smo, file = output.C_smo)
  write.table(logp.hat, file = output.phat)
 }
 if(VERBOSE)
  cat("Downscale to B2\n")
 if(!is.null(B2))
 {
  ind2 = sample(1:B, B2)
  particle = particle[ind2]
 }
 
 cat("Time before return:\n")
 show(Sys.time())
 
 if (output_latent)
  return(list(latent=particle,logp.hat=logp.hat,ltest=ltest,llogw=llogw,hospD.sim=hospD.sim,testD.sim=testD.sim,cp_prob=cp_probD))
  #return(list(latent=particle,param=res.param,logp.hat=logp.hat,ltest=ltest,llogw=llogw,hospD.sim=hospD.sim,testD.sim=testD.sim))
 else
  return(list(latent=NULL,logp.hat=logp.hat,ltest=ltest,hospD.sim=hospD.sim,testD.sim=testD.sim,cp_prob=cp_probD))
  #return(list(latent=NULL,param=res.param,logp.hat=logp.hat,ltest=ltest,hospD.sim=hospD.sim,testD.sim=testD.sim))
}
