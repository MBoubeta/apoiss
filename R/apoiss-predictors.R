#' Calculate the log of p_d.
#' 
#' @param X: Design matrix.
#' @param d: Domain index.
#' @param beta: Beta parameter.
#' @param phi: Variance component.
#' @param vd: Random effects i.i.d. N(0, 1).
#' @return The log of p_d.
#' @export

log_pd = function(X, d, beta, phi, vd) {
  
  log_pd = X[d, ] %*% beta + phi*vd
  return(log_pd)
}


#' Calculate the EBP of p_d.
#'
#' @param y: Response variable.
#' @param X: Design matrix.
#' @param beta: Beta parameter.
#' @param phi: Variance component.
#' @param nu: Offset parameter.
#' @param L: length of the simulated random effects.
#' @return The EBP of p_d.
#' @export
#' @examples
#' EBPp(y, X, beta, phi, nu, L)

EBPp = function(y, X, beta, phi, nu, L) {
  
  D = dim(X)[1]
  p_ebp = numeric(D)
  
  for (d in 1:D) {
    
    # generate random effects
    vd = rnorm(L, mean=0, sd=1)
    vd = c(vd,-vd)
    
    # calculate the components of the EBP pf pd (Nd and Dd) 
    logpd = log_pd(X, d, beta, phi, vd)
    Nd = sum(exp((y[d] + 1) * logpd - nu[d] * exp(logpd))) / (2*L)
    Dd = sum(exp(y[d] * logpd - nu[d] * exp(logpd))) / (2*L)
    p_ebp[d] = Nd / Dd
  }
  
  return(p_ebp)
}	


#' Calculate the EBP of v_d.
#'
#' @param d: Index of the domain.
#' @param y: Response variable.
#' @param X: Design matrix.
#' @param beta: Beta parameter.
#' @param phi: Variance component.
#' @param nu: Offset parameter.
#' @param L: length of the simulated random effects.
#' @return The EBP of v_d.
#' @export
#' @examples
#' EBPvd(d, y, X, beta, phi, nu, L)

EBPvd = function(d, y, X, beta, phi, nu, L) {
  
  vl = rnorm(L, mean=0, sd=1)
  vl = c(vl, -vl)
  N = sum(vl * exp(y[d] * (X[d, ] %*% beta + phi * vl) - nu[d] * exp(X[d, ] %*% beta + phi * vl))) / (2*L)
  D = sum(exp(y[d] * (X[d, ] %*% beta + phi * vl) - nu[d] * exp(X[d, ] %*% beta + phi * vl))) / (2*L)
  vd_ebp = N / D
  return(vd_ebp)
}


#' Calculate the EBP of v.
#'
#' @param y: Response variable.
#' @param X: Design matrix.
#' @param beta: Beta parameter.
#' @param phi: Variance component.
#' @param nu: Offset parameter.
#' @param L: length of the simulated random effects.
#' @return The EBP of v.
#' @export
#' @examples
#' EBPv(y, X, beta, phi, nu, L)

EBPv = function(y, X, beta, phi, nu, L) {
  
  D = dim(X)[1]
  d = 1:D
  v_ebp = sapply(d, EBPvd, y, X, beta, phi, nu, L)
  return(v_ebp)
}
  
  
#' Calculate the plug-in estimator of p.
#'
#' @param X: Design matrix.
#' @param beta: Beta parameter.
#' @param phi: Variance component.
#' @param v: Predicted random effects.
#' @return The plug-in estimator of p.
#' @export
#' @examples
#' PLUGINp(X, beta, phi, v)

PLUGINp = function(X, beta, phi, v) {
  
  p_plugin = exp(X %*% beta + phi * v)
  return(p_plugin)
}


