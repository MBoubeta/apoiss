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
  
  logpd = c(X[d, ] %*% beta) + phi*vd
  return(logpd)
}


#' Calculate the EBP of p_d or v_d.
#'
#' @param d: Index of the domain.
#' @param y: Response variable.
#' @param X: Design matrix.
#' @param beta: Beta parameter.
#' @param phi: Variance component.
#' @param nu: Offset parameter.
#' @param L: length of the simulated random effects.
#' @param parameter: 'target' for p_d and 'random_effects' for v_d.
#' @return The EBP of p_d or v_d.
#' @export
#' @examples
#' ebpd(d, y, X, beta, phi, nu, L, parameter)

ebpd = function(d, y, X, beta, phi, nu, L, parameter) {
  
  vd = rnorm(L, mean=0, sd=1)
  vd = c(vd, -vd)
  
  # calculate the components of the EBP pf pd (Nd and Dd) 
  logpd = log_pd(X, d, beta, phi, vd)
  
  if (parameter == 'target') {
    Nd = sum(exp((y[d] + 1) * logpd - nu[d] * exp(logpd))) / (2*L)
  } else if (parameter == 'random_effects') {
    Nd = sum(vd * exp(y[d] * logpd - nu[d] * exp(logpd))) / (2*L)
  } else {
    stop("parameter not valid!")
  }
  
  Dd = sum(exp(y[d] * logpd - nu[d] * exp(logpd))) / (2*L)
  ebpd = Nd / Dd
  return(ebpd)
}


#' Calculate the EBP of p or v.
#'
#' @param y: Response variable.
#' @param X: Design matrix.
#' @param beta: Beta parameter.
#' @param phi: Variance component.
#' @param nu: Offset parameter.
#' @param L: length of the simulated random effects.
#' @param parameter: 'target' for p and 'random_effects' for v.
#' @return The EBP of p or v.
#' @export
#' @examples
#' ebp(y, X, beta, phi, nu, L, parameter)

ebp = function(y, X, beta, phi, nu, L, parameter) {
  
  D = dim(X)[1]
  d = 1:D
  
  ebp = sapply(d, ebpd, y, X, beta, phi, nu, L, parameter)
  return(ebp)
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
#' plugin(X, beta, phi, v)

plugin = function(X, beta, phi, v) {
  
  plug = exp(X %*% beta + phi * v)
  return(plug)
}


