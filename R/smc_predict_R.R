#' Prediction of R based on output from smc_loglik_mem. 
#' Assumes all inference results are stored in files that are read within the routines.
#' Parameter values from the last inference time points are used (one set for each particle)
#' One routine for each model for R.
#'
#' @param 'dir' Directory (relative to working directory) for files with simulations of R
#' @param 'npred' number of days to do prediction
#' @return A matrix containing both the simulated R's before the prediction date (based on fixed-lag smoothing) and the predicted values
smc_predict_R_rw = function(dir,npred)
{
  R_sim = as.matrix(read.table(paste(dir,"smc_dynamic_R_lag.txt",sep="")))
  sig.R = as.matrix(read.table(paste(dir,"smc_SdR.txt1",sep="")),header=T)
  B = ncol(R_sim)
  #Remove last lines, NA for some reason
  n = nrow(R_sim)
  if(is.na(R_sim[n,1]))    #Error in some earlier runs
  {
   R_sim = R_sim[-n,]
   sig.R = sig.R[-n,]
  }
  n = nrow(R_sim)
  sig.R = sig.R[n,]
  #Predictions
  for(j in 1:npred)
  {
    R_sim = rbind(R_sim,NA)
    R_sim[n+j,] = R_sim[n+j-1,]*exp(sig.R*rnorm(B))
  }
  R_sim
}

smc_predict_R_ar = function(dir,npred)
{
  R_sim = as.matrix(read.table(paste(dir,"smc_dynamic_R_lag.txt",sep="")))
  a.R = as.matrix(read.table(paste(dir,"smc_SdR.txt1",sep="")),header=T)
  sig.R = as.matrix(read.table(paste(dir,"smc_SdR.txt2",sep="")),header=T)
  mu.R = as.matrix(read.table(paste(dir,"smc_SdR.txt3",sep="")),header=T)
  B = ncol(R_sim)
  #Remove last lines, NA for some reason
  n = nrow(R_sim)
  if(is.na(R_sim[n,1]))    #Error in some earlier runs
  {
    R_sim = R_sim[-n,]
    a.R = a.R[-n,]
    sig.R = sig.R[-n,]
    mu.R = mu.R[-n,]
  }
  n = nrow(R_sim)
  a.R = a.R[n,]
  sig.R = sig.R[n,]
  mu.R = mu.R[n,]
  #Predictions
  for(j in 1:npred)
  {
    R_sim = rbind(R_sim,NA)
    R_sim[n+j,] = exp(mu.R+a.R*(log(R_sim[n+j-1,])-mu.R)+sig.R*rnorm(B))
  }
  R_sim
}

smc_predict_R_cp = function(dir,npred)
{
  R_sim = as.matrix(read.table(paste(dir,"smc_dynamic_R_lag.txt",sep="")))
  phi.R = as.matrix(read.table(paste(dir,"smc_SdR.txt1",sep="")),header=T)
  sig.R = as.matrix(read.table(paste(dir,"smc_SdR.txt2",sep="")),header=T)
  B = ncol(R_sim)
  #Remove last lines, NA for some reason
  n = nrow(R_sim)
  if(is.na(R_sim[n,1]))    #Error in some earlier runs
  {
    R_sim = R_sim[-n,]
    phi.R = phi.R[n-1,]
    sig.R = sig.R[n-1,]
  }
  n = nrow(R_sim)
  phi.R = phi.R[n,]
  sig.R = sig.R[n,]
  #Predictions
  for(j in 1:npred)
  {
    R_sim = rbind(R_sim,NA)
    cp = rbinom(B,1,phi.R)
    R_sim[n+j,] = R_sim[n+j-1,]*exp(cp*sig.R*rnorm(B))
  }
  R_sim
}