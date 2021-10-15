#' Particle simulation of the spread model and log-likelhood evaluation based on hospital data 
#' Mainly to be used within smc_loglik
#
#' Global variables used:
#' 'mobility_matrices', assumed to start at 'input_date'
#' 'hospD'
#' 'testD'
#'
#' First particles are 'unit.samp' days from 'input_date'
#' Thereafter log-likelhood is calculated based on hospital data
#'
#' There are currently many global parameter used in the underlying routines. No control on that
#'
#' @param 'input' A list with 'particle' to start from, 'Rsim', a vector with R-valyes for the differen days 'AMP' and 'inshist'
#' @param 'input_date' start date for simulation
#' @param 'loc-seed' seeding for call to run_model, should be adjusted for input_date
#' @param 'mob_mat' mobility matrix, should be adjusted for input_date
#' @param 'unit.samp' number of days to simulate
#' @param 'day.numb' The day number for input_date
#' @param 'hospD' hospital data, assumed to start at initial date at whole run (not current input_date)
#' @param 'testD' test data, assumed to start at initial date at whole run (not current input_date)
#' @param 'pbar' Survival probability for delay until hospitalization
#' @param 'qbar' Survival probability for delay until testing
#' @param 'prob_h' Probability for hospitalisation -> NOW A VECTOR WITH LENGTH=NUMBER OF DAYS - NOTE THAT beta_hosp is also a vector then
#' @param 'alpha_hosp' Parameter in beta-binomial distribution for hospital data
#' @param 'alpha_test' Parameter in beta-binomial distribution for test data
#' @param 'pi_test' Parameters describing detection probability
#' @param 'MaxdaysH' Maximum number of days delay until hospitalization
#' @param 'MaxdaysT' Maximum number of days delay until testing
#' @param 'sim.data' If true, observed data are replaced by simulated ones
#' @return A list containing simulated particles of latent variables ('latent'), 'ins' and corresponding log-weights
smc_particle_sim = function(input,input_date,loc.seed,mob_mat,
                            unit.samp,day.numb,hospD,testD,pbar,qbar,
                            prob_h,alpha_hosp,alpha_test,pi_test,MaxdaysH,MaxdaysT,
                            sim.data = FALSE)
{
  D.sim = NULL
  if(sim.data)
    D.sim = matrix(NA,nrow=unit.samp,ncol=2)
 #print("smc_particle_sim")
 #show(input_date+1:unit.samp)
 #If dynamic R, use Reffs to define the different R's
 if(!is.null(alpha_hosp))   #Use beta-binomial distribution
  beta_hosp = alpha_hosp * (1 - prob_h) / prob_h
 if(unit.samp>1)
 {
  ref_change_dates=seq(input_date+1,to=input_date+unit.samp-1,by=1)
  #ref_change_dates=seq(input_date+2,to=input_date+unit.samp,by=1)
  R_eff = input$R_sim[1+1:(unit.samp-1)]
 }
 else
  R_eff = c()
 p.new <- smc_run_model(mob_mat, 
                    start_date = input_date, 
                    end_date = input_date+unit.samp,
                    N_sim = 1,
                    R_0 = input$R_sim[1],
                    reff_change_dates=ref_change_dates,
                    Reffs = R_eff,
                    parallel = F,
                    seeding=loc.seed,
                    se1e2iiar_pop=input$particle,
                    AMP_factor = input$AMP)[[1]]
 #p.new$date = p.new$date+day.numb-1
 #p.new$day = p.new$day+day.numb-1
 #Calculate log-likelihood
 logw.new = 0
 #Simulate days for hospitalisation for all incidences
 ins = p.new[,.(ins=sum(c_symp_incidence+c_asymp_incidence)),by=day]$ins
 insH = input$inshist
 insT = input$instest
 N = length(ins)
 M = nrow(insH)
 M2 = nrow(insT)    #Should be same as M?
 if(abs(M-M2)>0)
 {
  print("smc_particle_sim M,M2")
  print("M and M2 are different!")
  show(M)
  show(M2)
  show(insT)
 }
 for(u in 1:unit.samp)
 {
  insH = rbind(insH,NA)
  insH[M+u,1] = ins[u]
  #for(j in 0:min(MaxdaysH-1,M+u-1))
  for(j in 0:(MaxdaysH-1))
  {
   #        insH[M+u-j,2+j] = rbinom(1,insH[M+u-j,1],pbar[j+1])
   #        insH[M+u-j,1] = insH[M+u-j,1]-insH[M+u-j,2+j]
   insH[M+u,2+j] = rbinom(1,insH[M+u,1],pbar[j+1])
   insH[M+u,1] = insH[M+u,1]-insH[M+u,2+j]
  }
  insT = rbind(insT,NA)
  insT[M+u,1] = ins[u]
  for(j in 0:(MaxdaysT-1))
  {
   insT[M+u,2+j] = rbinom(1,insT[M+u,1],qbar[j+1])
   insT[M+u,1] = insT[M+u,1]-insT[M+u,2+j]
  }
 }
 #Likelihood calculation for hospital data
 N = nrow(insH)
 num = rep(NA,unit.samp)
 for(i in 1:unit.samp)
 {
  if(!is.na(hospD$N[i]))
  {
   NN = min(MaxdaysH,N-unit.samp+i)
   ind = 1+c(1:NN)*(N-1)  #Index to pick out right values of insH
   #num[i] = sum(insH[M+unit.samp+1-i+ind],na.rm=TRUE)
   num[i] = sum(insH[M+i+ind],na.rm=TRUE)
   if(hospD$N[i]>num[i])
    logw.new =  -9999                #Somewhat ad hoc
   else
   {
    if(is.null(alpha_hosp))
     l = dbinom(hospD$N[i],num[i],prob_h[i],log=TRUE)
    else
     l = dbbinom(hospD$N[i],num[i],
                 alpha=alpha_hosp,beta=beta_hosp[i],log=TRUE)
    if(is.na(l))
     show(c(M+i,hospD$N[i],num[i],beta_hosp[i]))
    logw.new = logw.new+l
   }
   if(is.na(logw.new))
    print("smc_particle_sim: Something is wrong")
  }
  if(sim.data)
  {
    NN = min(MaxdaysH,N-unit.samp+i)
    ind = 1+c(1:NN)*(N-1)  #Index to pick out right values of insH
    #num[i] = sum(insH[M+unit.samp+1-i+ind],na.rm=TRUE)
    num[i] = sum(insH[M+i+ind],na.rm=TRUE)
    if(is.null(alpha_hosp))
    D.sim[i,1] = rbinom(1,num[i],prob_h[i])
   else
    D.sim[i,1] = rbbinom(1,num[i],alpha=alpha_hosp,beta=beta_hosp[i])
   
  }
  
 }
 #print("smc_particle_sim: loglik hospD:")
 #show(logw.new)
 
 #likelihood calculation for test data
 l = NA
 ll = rep(NA,unit.samp)
 numT=NA
 if(!is.null(testD))
 {
  N = nrow(insT)
  numT = rep(NA,unit.samp)
  for(i in 1:unit.samp)
  {
   NN = min(MaxdaysT,N-unit.samp+i)
   ind = 1+c(1:NN)*(N-1)
   #show(c(M+i,testD$N[i], ((M+i)%%7)))
   if(!is.na(testD$Positive[i]))
   {       
    #numT[i] = sum(insT[M+unit.samp+1-i+ind],na.rm=TRUE)
    numT[i] = sum(insT[M+i+ind],na.rm=TRUE)
    eta.t = pi_test[1]+pi_test[2]*testD$N[i]
    pi.t = exp(eta.t)/(1+exp(eta.t))
    #print("sim_particle_sim")
    #show(c(M+i,numT[i]))
    #show(testD[i,])
    if(testD$Positive[i]>numT[i])
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
      l = dbinom(testD$Positive[i],numT[i],pi.t,log=TRUE)
     else
      l = dbbinom(testD$Positive[i],numT[i],
                  alpha=alpha_test,beta=alpha_test*(1-pi.t)/pi.t,log=TRUE)
     if(is.na(l))
      show(c(l))
    }
    logw.new = logw.new+l
    if(is.na(logw.new))
     print("smc_particle_sim: Something is wrong")
    ll[i] = l
   }
   if(sim.data & !is.na(testD$N[i]))
   {
     numT[i] = sum(insT[M+i+ind],na.rm=TRUE)
     eta.t = pi_test[1]+pi_test[2]*testD$N[i]
     pi.t = exp(eta.t)/(1+exp(eta.t))
     if(is.null(alpha_hosp))
       D.sim[i,2] = rbinom(1,numT[i],pi.t)
    else
     D.sim[i,2]= rbbinom(1,numT[i], alpha=alpha_test,beta=alpha_test*(1-pi.t)/pi.t)
   }
  }
 }
 return(copy(list(particle.new=p.new,insH=insH,numH=num,numT=numT,insT=insT,logw.new=logw.new,ltest=ll,
        D.sim=D.sim)))
}
