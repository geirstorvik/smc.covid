#' Resampling using 
#' @param weights Vector proportional to probabilities
#' @param N Number of samples, default to length of weights.
#' @return Vector of samples indices 
smc_resamp = function(weights,N=length(weights))
{
 #Normalize weights
 weights = weights/sum(weights)
 M = length(weights)
 #Calculate deterministic number proportional to probabilities
 N1= as.integer(N*weights)
 ind = rep(1:M,N1)
 #Sampling remaining part
 L = N-sum(N1)
 if(L>0)
 {
  q = N*weights-N1
  ind = c(ind,sample(1:M,L,prob=q,replace=T))
 }
 ind
}
