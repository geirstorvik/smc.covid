#' Initial simulation of the spread model 
#' Mainly to be used within smc_loglik
#
#'
#' First initial particles are generated until an initial date
#'
#' There are currently many global parameter used in the underlying run_model routine. No control on that
#'
#' @param 'input' A list with 'R0', 'AMP'
#' @param 'start_date' start date for simulation
#' @param 'end_date' end date for simulation
#' @param 'mob_mat' Mobility matrix, supposed to start at start_date
#' @param 'loc.seed' seeding for call to run_model
#' @param 'pop' population structure
#' @param 'MaxdaysH' maximum delay days for hospitalization
#' @param 'pbar' Survival probability for day of hospitalization
#' @param 'MaxdaysT' maximum delay days for testing
#' @param 'qbar' Survival probability for day of testing
#' @return A list containing 'particle', simulated particles of latent variables, 'insH', simulated numbers of potential Hospitalizations for each day (first index is infection day, second index is delay, starting at zero), 'insT' is delay of potential test with same structure as 'insH'
smc_init_sim = function(input,start_date,end_date,mob_mat,loc.seed,pop,
                        MaxdaysH,pbar,MaxdaysT,qbar){
  particle <- smc_run_model(mob_mat, 
                          start_date = start_date-1, 
                          end_date = end_date,  
                          N_sim = 1,
                          R_0 = input$R0,
                          parallel = F,
                          seeding=loc.seed,
                          se1e2iiar_pop=pop,
                          AMP_factor = input$AMP)[[1]]
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
  return(copy(list(particle=particle,insH=insH,insT=insT)))
}