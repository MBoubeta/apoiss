Upql = function(xi, X, phi, nu, y) {
  
  D = dim(X)[1]
  p = dim(X)[2]
  beta = xi[1:p]
  v = xi[(p + 1):(p + D)]
  
  # calculating of U_r, r = 1, ..., p
  Ur_pql = function(r, y, X, beta, phi, nu, v) {
    Ur = -sum(nu*exp(X%*%beta + phi*v)*X[, r]) + sum(y*X[, r])
    return(Ur)
  }
  
  r = 1:p
  Ur = sapply(r, Ur_pql, y, X, beta, phi, nu, v)
  
  # calculating of U_{p+d}, d = 1, ..., D
  Upd = -v - nu*exp(X%*%beta + phi*v)*phi + y*phi
  
  U = c(Ur, Upd)
  return(U)
}


Hpql = function(xi, X, phi, nu, y) {
  
  D = dim(X)[1]
  p = dim(X)[2]
  beta = xi[1:p]
  v = xi[(p + 1):(p + D)]
  
  aux = exp(X%*%beta + phi*v)
  
  # Hessian matrix
  Hpql = matrix(0, nrow=p+D, ncol=p+D)  
  r1r2 = combn(p, 2)
  r1r2 = matrix(data=r1r2, ncol=2, byrow=TRUE)
  
  Hr1r2_pql = function(X, nu, aux, r1, r2) {
    Hr1r2 = -sum(nu*aux*X[, r1]*X[, r2])
    return(Hr1r2)
  }
  
  # H_{r1,r2}, r1, r2=1, ..., p
  for(i in 1:dim(r1r2)[1]){
    Hpql[r1r2[i, 1], r1r2[i, 2]] = Hr1r2_pql(X, nu, aux, r1=r1r2[i, 1], r2=r1r2[i, 2])
  }
  
  # H_{r, p+d}, r=1, ..., p and d=1, ..., D
  for (r in 1:p) {
    Hpql[r, (p+1):(p+D)] = -nu*aux*X[, r]*phi
  }
  
  Hpql = Hpql + t(Hpql)
  
  # H_{r,r}, r = 1, ..., p
  Hrr_pql = function(r, X, nu, aux) {
    Hrr = -sum(nu*aux*X[, r]^2)
    return(Hrr)
  }
  
  r = 1:p
  Hrr = sapply(r, Hrr_pql, X, nu, aux)
  
  # H_{p+d, p+d}, d=1, ..., D
  Hpdpd = -1-nu*aux*phi^2
  diag(Hpql) = c(Hrr, Hpdpd)
  return(Hpql)
}


PQL = function(xi0, X, phi0, nu, y, L, tol, maxIter, epsXiPhi, ...) {
  
  # parameters
  D = dim(X)[1]
  p = dim(X)[2]
  xiNew = matrix(data=0, nrow=L, ncol=length(xi0))
  phiNew = numeric(L)
  normXiPhi = numeric(L)
  
  for (ell in 1:L) {
    xiRes = nleqslv(x=xi0, Upql, jac=Hpql, X, phi0, nu, y, ...)
    xiNew[ell, ] = xiRes$x
    vNew = xiNew[ell, (p + 1):(p + D)]
    
    # Newton-Raphson updating equation
    phiNew[ell] = sqrt(mean(vNew^2) * phi0^2)
    
    # convergence
    normXiPhi[ell] = max(abs(c(xiNew[ell, ], phiNew[ell]) - c(xi0, phi0)))
    convXiPhi = (normXiPhi[ell] < epsXiPhi)
    
    if (convXiPhi == TRUE) {
      break
    } else {
      xi0 = xiNew[ell, ]
      phi0 = phiNew[ell]
    }
  }
  
  # PQL results
  PQLres = list('xi'=xiNew[ell, ], 'phi'=phiNew[ell], 'xiRes'=xiNew[1:ell, ], 'phiRes'=phiNew[1:ell],
                'normXiPhi'=normXiPhi[1:ell], 'iter'=ell)
  return(PQLres)
}
