#' Calculate a natural set of equations to apply the Method of Moments (MM) estimator.
#' 
#' @param theta: Model parameters.
#' @param X: Design matrix.
#' @param nu: Offset parameter.
#' @param y: Response variable.
#' @return A natural set of equations.
#' @export
#' @examples
#' f([-2, .6], 1)

f = function(theta, X, nu, y) {
  
  p = length(theta) - 1
  beta = theta[1:p]
  phi = theta[p+1]
  
  # calculate f_k, k = 1, ..., p
  fk = function(k, y, X, beta, phi, nu) {
    f_k = sum(nu * exp(X %*% beta + 0.5 * phi^2) * X[, k]) - sum(y * X[, k])
    return(f_k)
  }
  
  k = 1:length(beta)
  f_k = sapply(k, fk, y, X, beta, phi, nu)
  
  # calculate f_{p+1}
  f_p1 = sum(nu * exp(X %*% beta + 0.5 * phi^2) + nu^2 * exp(2 * X %*% beta + 2 * phi^2)) - sum(y^2)
  
  f = c(f_k, f_p1)
  return(f)
  
}


#' Calculate the Jacobian matrix elements.
#' 
#' @param X: Design matrix.
#' @param nu: Offset parameter.
#' @param k: index of the function to calculate the partial.
#' @param r: index of the beta parameter to calculate the partial.
#' @param aux1: aux1 = exp(X %*% beta.hat + 0.5 * phi.hat^2).
#' @return \partial f_k / \partial \beta_r.
#' @export

partial_fk_br = function(X, nu, k, r, aux1) {
  
  partial_fk_br = sum(nu * aux1 * X[, r] * X[, k])
  return(partial_fk_br)
  
}


#' Calculate the Jacobian matrix elements.
#' 
#' @param X: Design matrix.
#' @param nu: Offset parameter.
#' @param k: index of the function to calculate the partial.
#' @param phi: variance component.
#' @param aux1: aux1 = exp(X %*% beta.hat + 0.5 * phi.hat^2).
#' @return \partial f_k / \partial \phi.
#' @export

partial_fk_phi = function(X, nu, k, phi, aux1) {
  
  partial_fk_phi = sum(nu * aux1 * phi * X[, k])
  return(partial_fk_phi)
  
}


#' Calculate the Jacobian matrix elements.
#' 
#' @param X: Design matrix.
#' @param nu: Offset parameter.
#' @param r: index of the beta parameter to calculate the partial.
#' @param aux1: aux1 = exp(X %*% beta.hat + 0.5 * phi.hat^2).
#' @param aux2: aux2 = exp(2*X %*% beta.hat + 2*phi.hat^2).
#' @return \partial f_k / \partial \phi.
#' @export

partial_fp1_br = function(X, nu, r, aux1, aux2) {
  
  partial_fp1_br = sum(nu * aux1 * X[, r] + 2 * nu^2 * aux2 * X[, r])    
  return(partial_fp1_br)
  
}


#' Calculate the Jacobian matrix elements.
#' 
#' @param phi: variance component.
#' @param nu: Offset parameter.
#' @param aux1: aux1 = exp(X %*% beta.hat + 0.5 * phi.hat^2).
#' @param aux2: aux2 = exp(2*X %*% beta.hat + 2*phi.hat^2).
#' @return \partial f_p+1 / \partial \phi.
#' @export

partial_fp1_phi = function(phi, nu, aux1, aux2) {
  
  partial_fp1_phi = sum(nu * aux1 * phi + 4*nu^2 * aux2 * phi)
  return(partial_fp1_phi)
  
}


#' Calculate the Jacobian matrix.
#' 
#' @param theta: Model parameters.
#' @param X: Design matrix.
#' @param nu: Offset parameter.
#' @param y: Response variable.
#' @return Jacobian matrix.
#' @export

H = function(theta, X, nu, y) {
  
  p = length(theta) - 1
  beta = theta[1:p]
  phi = theta[p + 1]
  
  H = matrix(0, nrow=p+1, ncol=p+1)
  comb = c(gl(p, 2), combn(p, 2))
  comb = matrix(data=comb, ncol=2, byrow=TRUE)
  
  aux1 = exp(X%*%beta + 0.5*phi^2)
  aux2 = exp(2*X%*%beta + 2*phi^2)
  
  for(i in 1:dim(comb)[1]) {
    H[comb[i,1], comb[i,2]] = partial_fk_br(X, nu, k=comb[i,1], r=comb[i,2], aux1)
    H[comb[i,2], comb[i,1]] = H[comb[i, 1], comb[i,2]]
  }
  
  for(r in 1:p) {
    H[r, p+1] = partial_fk_phi(X, nu, k=r, phi, aux1)
    H[p+1, r] = partial_fp1_br(X, nu, r, aux1, aux2)
  }
  
  H[p+1, p+1] = partial_fp1_phi(phi, nu, aux1, aux2)
  
  return(H)
}


mm_fit = function(theta0, X, nu, y, ...) {
  
  # MM fit
  mm_est = nleqslv(x=theta0, f, jac=H, X, nu, y, ...)
  
  return(mm_est)
}
  
  
#' Calculate the variance-covariance matrix of (\hat{\beta}_0, \hat{\beta}_1, \hat{\phi}) by bootstrap.
#'
#' @param X: Design matrix.
#' @param theta: Model parameters.
#' @param nu: Offset parameter.
#' @param maxiter: Maximum number of iterations for the Newton-Raphson algorithm.
#' @param tol: Tolerance parameter for the Newton-Raphson algorithm.
#' @param B: Number of bootstrap resamples.
#' @return Variance-covariance matrix.
#' @export

varcov = function(X, theta, nu, B, ...) {
  
  # set parameters  
  D = dim(X)[1]
  p = dim(X)[2]
  beta = theta[1:p]
  phi = theta[p + 1]
  theta_boot = matrix(data=NA, ncol=B, nrow=p+1)
  y_boot = numeric(D)
  
  # names of the covariates
  covars = colnames(X)
  
  # fit the model with or without (Intercept)  
  if (covars[1] == "(Intercept)") {
    covars = covars[covars != "(Intercept)"]
    formula = as.formula(paste0("y_boot ~ ", paste0(covars, collapse= "+")))    
  } else {
    formula = as.formula(paste0("y_boot ~ ", paste0(covars, collapse= "+"), "-1"))
  }
  
  # bootstrap resamples
  for (b in 1:B) {
    repeat {
      # TODO: review boot
      yboot_res = yboot(X, beta, phi, nu)
      y_boot = yboot_res$y
      
      # fit the model log(p_{d}) = \bbeta \xx_{d} + \phi v_{d} (d = domain) by MM.
      beta0 = glm(formula=formula, family=poisson, offset=log(nu))$coefficients
      eta_tilde = X %*% beta0
      eta_dir = log((y_boot + 1) / (nu + 1))
      phi0 = sqrt(sum((eta_tilde - eta_dir)^2) / D)
      theta0 = c(beta0, phi0)
      
      mm_est = try(mm_fit(theta0, X, nu, y_boot, ...), TRUE)
      
      if (class(mm_est) != 'try-error') {
        theta_boot[, b] = c(mm_est$x[1:p], mm_est$x[p + 1])
        break
      }
    }
  }
  
  # calculate the mean of theta
  theta_boot_avg = apply(theta_boot, 1, mean)
  
  f = function(x) {x%*%t(x)}
  varcovb = apply(theta_boot - theta_boot_avg, 2, f)
  
  # variance-covariance matrix
  varcov = matrix(rowSums(varcovb), ncol=p+1, nrow=p+1)/B
  return(varcov)
}


#' Add code for the estimation of model parameters \hat{\beta}_k, k=0, ..., p,
#'
#' @param i: Index of beta parameter
#' @param pvalue: pvalue.
#' @return The code for the MM estimations of beta parameter.
#' @export

add_code = function(i, pvalue) {
  
  code = character()
  if (pvalue[i] < 0.001) {code = '***'}
  if (0.001 <= pvalue[i] & pvalue[i] < 0.01) {code = '**'}
  if (0.01 <= pvalue[i] & pvalue[i] < 0.05) {code = '*'}
  if(0.05 <= pvalue[i] & pvalue[i] < 0.1) {code = '.'}
  if(0.1 <= pvalue[i]) {code = ' '}
  return(code)
}

#' Estimation of model parameters (\hat{\beta}_k, k=0, ..., p, \hat{\phi}) using the MM fitting algorithm.
#'
#' @param y: Response variable.
#' @param X: Design matrix.
#' @param beta0: Seed for beta parameter.
#' @param phi0: Seed for variance component.
#' @param nu: Offset parameter.
#' @param maxiter: Maximum number of iterations for the Newton-Raphson algorithm.
#' @param tol: Tolerance parameter for the Newton-Raphson algorithm.
#' @param B: Number of bootstrap resamples.
#' @param add_std_error: If the std error is included.
#' @return The estimations of model parameters using the MM fitting algorithm.
#' @export

mm = function(y, X, beta0, phi0, nu, B, add_std_error=FALSE, ...) {
  
  # parameter definition
  p = length(beta0)
  theta0 = c(beta0, phi0)
  
  # MM fit
  mm_est = mm_fit(theta0, X, nu, y, ...)
  
  iter = mm_est$iter
  beta = mm_est$x[1:p]
  phi = mm_est$x[p + 1]
  theta = c(beta, phi)
  
  # names of beta
  names(beta) = colnames(X)
  
  mm_est = list('coefficients' = beta, 'random_effects' = phi)
  
  if (add_std_error == TRUE) {
    # variance-covariance matrix of (\hat{\beta}_0, \hat{\beta}_1, \hat{\phi})
    var_cov = varcov(X, theta, nu, B, ...)
    std_err = sqrt(diag(var_cov))[1:p]
    
    # results
    zvalue = beta/std_err
    pvalue = 2*(1 - pnorm(abs(zvalue)))
    
    i = 1:p
    code = sapply(i, add_code, pvalue)
  
    cat('Area-level Poisson mixed model by MM \n')
    
    # random effects
    cat('Random effects: \n')
    rand_effs = data.frame('(Intercept)', phi, row.names='')
    names(rand_effs) = c('Random effect', 'phi')
    print(rand_effs)
    
    # fixed effects
    cat('\n Fixed effects: \n')
    
    fixed_effs = data.frame(beta, std_err, zvalue, pvalue, code, row.names=colnames(X))
    names(fixed_effs) = c('Coefficients', 'Std. Error', 'z value', 'Pr(>|z|))', ' ')
    print(fixed_effs)
    
    mm_est_std = list('std_error'=std_err, 'zvalue'=zvalue, 'pvalue'=pvalue, 'iter'=iter)
    mm_est = append(mm_est, mm_est_std)
  }
  
  return(mm_est)  
}


