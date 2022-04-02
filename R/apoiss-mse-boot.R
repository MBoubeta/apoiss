mse_boot = function(X, beta, phi, nu, L, max_iter, tol, B, ...) {

  D = dim(X)[1]
  p = dim(X)[2]
  
  # names of the covariates
  covars = colnames(X)
  
  # fit the model with or without (Intercept)  
  if (covars[1] == "(Intercept)") {
    covars = covars[covars != "(Intercept)"]
    formula = as.formula(paste0("y_boot ~ ", paste0(covars, collapse= "+")))    
  } else {
    formula = as.formula(paste0("y_boot ~ ", paste0(covars, collapse= "+"), "-1"))
  }
  
  mseboot = matrix(data=NA, nrow=B, ncol=D)
  pboot = matrix(data=NA, nrow=B, ncol=D)
  ebpboot = matrix(data=NA, nrow=B, ncol=D)
  
  for (b in 1:B) {
    yboot_res = yboot(X, beta, phi, nu)
    y_boot = yboot_res$y
    p_boot[b, ] = yboot_res$pd
    
    beta0 = glm(formula=formula, family=poisson, offset=log(nu))$coefficients
    eta_tilde = X %*% beta0
    eta_dir = log((y_boot + 1)/(nu + 1))
    phi0 = sqrt(sum((eta_tilde - eta_dir)^2)/D)
    
    mm_fit = mm(y_boot, X, beta0, phi0, nu, add_std_error=FALSE, ...)
    
    beta_boot = mm_fit$x[1:p]
    phi_boot = mm_fit$x[p + 1]
    
    # EBP bootstrap
    ebp_boot[b, ] = EBP(y_boot, X, beta_boot, phi_boot, nu, L, parameter='target')
    mseboot[b, ] = ebpboot[b, ] - pboot[b, ]
  }
  
  # MSE bootstrap
  mseboot = colSums(mseboot ^ 2)/B
  
  # p^{*}_d y \hat{p}^{*}_d.
  pboot = apply(pboot, 2, mean)
  ebpboot = apply(ebpboot, 2, mean)
  
  mseboot_res = list('mse_boot'=mseboot, 'p_boot'=pboot, 'ebp_boot'=ebpboot)
  return(mseboot_res)
}

