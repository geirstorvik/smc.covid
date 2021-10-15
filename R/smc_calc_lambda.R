smc_calc_lambda = function(lambda,lambda_size,MaxdaysH=NULL)
{
  if(is.null(lambda_size))
  {
    if(is.null(MaxdaysH))
      MaxdaysH = min(c(1:100)[ppois(1:100,lambda)>0.99])
    pbar = dpois(0:MaxdaysH,lambda)
    pbar = pbar/(sum(pbar))
#    pbar = dpois(0:MaxdaysH,lambda)/(1-ppois(c(0:MaxdaysH)-1,lambda))
  }
  else
  {
    pr = lambda_size/(lambda+lambda_size)
    if(is.null(MaxdaysH))
      MaxdaysH = min(c(1:100)[pnbinom(1:100,lambda_size,pr)>0.99])
    pbar = dnbinom(0:MaxdaysH,lambda_size,pr)
    pbar = pbar/sum(pbar)
#      (1-pnbinom(c(0:MaxdaysH)-1,lambda_size,pr))
  }  
  for(i in 2:length(pbar))
    pbar[i] = pbar[i]/sum(pbar[i:length(pbar)])
  pbar
}