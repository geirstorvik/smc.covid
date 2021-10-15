#' Modify test data
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
#' @param 'testD' A data.frame/table with columns 'Negative', 'Positive' and 'date'
#' @param 'pi' A vector giving pi_0 and pi_1 for specifying probability for positiv detection
#' @param 'date' A vector for which dates to construct new test data object
#' @return A data.frame with columns date and for every 7 day the sum of the previous positives and totals
smc_convert_test_data = function(testD,date,pi.test=c(-1.175,7.4e-05),week=FALSE)
{
  d = data.frame(date=date,N=NA,Positive=NA,d2=as.Date(NA))
  n = ncol(testD)
  matchdate = pmatch(testD$date,date)
  m2 = matchdate[!is.na(matchdate)]
  if("Tot" %in% names(testD))
    d[m2,]$N = testD$Tot[!is.na(matchdate)]
  else
    d[m2,]$N = testD$Negative[!is.na(matchdate)]+testD$Positive[!is.na(matchdate)]
  d[m2,]$Positive = testD$Positive[!is.na(matchdate)]
  d[m2,]$d2 = testD$date[!is.na(matchdate)]
  scale = 1
  if(week)
   {
      #Take sum of 7 days
      d$N =  rowSums(matrix(unlist(shift(d$N,0:6)),ncol =7))
      d$Positive =  rowSums(matrix(unlist(shift(d$Positive,0:6)),ncol =7))
      scale = 7
    }
  #Calculate pi_t, note that pi_1 is scaled down with 1 or 7 taking into account sum, not mean
  #eta = pi.test[1]+pi.test[2]*d$N/scale
  #d$pi_t = exp(eta)/(1+exp(eta))
  d
}