smc_trans_to_prob = function(l)
{
  l = exp(l-max(l,na.rm=TRUE))
  l/sum(l,na.rm=TRUE)
}
