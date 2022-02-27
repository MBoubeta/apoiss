ysim = function(X, beta, phi, nu) {
  
  D = dim(X)[1]
  
  repeat{
    # random effects
    v = rnorm(D, mean=0, sd=1)
    
    # response
    eta = X%*%beta + phi*v
    pd = exp(eta)
    
    if (sum(pd > 1) ==0 ) {break}
  }	
  
  y = rpois(D, nu*pd)
  res = list('y'=y, 'pd'=pd)
  return(res)
}
