#' @export
smc_smooth_week = function(mat,m=7)
{
  mat.week = mat
  n = nrow(mat)
  B = ncol(mat)
  l = floor(m/2)
  if(l==m/2)
    u = l-1
  else
     u = l
  for(i in 1:n)
  {
    l2 = max(1,i-l)
    u2 = min(n,i+u)
    res = rep(0,B)
    for( v in l2:u2)
      res = res + mat[v,]
    res = res/length(l2:u2)
    mat.week[i,] = res
    #mat.week[i,] = apply(mat[l2:u2,,drop=FALSE],2,mean)
  }
  mat.week  
}