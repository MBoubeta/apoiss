# check the behaviour of the apoiss fitting algorithms
# model: log(\mu_{d}) = \beta\xx_{d} + \phi * v_{d}

# libraries
library(yaml)
library(lme4)
library(nleqslv)


# source code
source('R/metrics.R')
source('R/apoiss-MM-fitting-algorithm.R')
source('R/apoiss-PQL-fitting-algorithm.R')
source('simulations/apoiss/ysim.R')

# load config file
config = read_yaml('simulations/apoiss/config.yml')

# inputs
# D = domains
D = config$D
d = seq(1:D)
dd = as.factor(d)

# p = beta parameters
p = config$p
beta0 = config$beta0
beta1 = config$beta1
beta = rbind(beta0, beta1)

# variance component
phi = config$phi

theta = c(beta, phi)

# offset parameter
nud = config$nud
nu = rep(nud, D)

# design matrix
x0 = rep(1, D)
x1 = d/D 
X = matrix(data=c(x0, x1), nrow=D, ncol=p, byrow=FALSE, dimnames=c(list(NULL, c('(Intercept)', 'x1'))))

# config parameters
K = config$K
tol = config$tol
maxiter = config$maxiter
L = config$L
B = config$B
S = config$S

# define parameters
beta_la = matrix(data=NA, nrow=K, ncol=p)
beta_mm = matrix(data=NA, nrow=K, ncol=p)
beta_pql = matrix(data=NA, nrow=K, ncol=p)

phi_la = numeric(K)
phi_mm = numeric(K)
phi_pql = numeric(K)


for (k in 1:K) {
  
  tit = proc.time()
 
  # response
  y_res = ysim(X, beta, phi, nu)
  pd = y_res$pd
  y = y_res$y 
  
  # laplace
  fit_la = glmer(formula=y~x1+(1|dd), family=poisson, nAGQ=1, offset=log(nu))
  beta_la[k, ] = fit_la@beta
  phi_la[k] = fit_la@theta
  
  # MM
  beta0 = glm(formula=y~x1, family=poisson, offset=log(nu))$coefficients
  eta_tilde = X%*%beta0
  eta_dir = log((y + 1)/(nu + 1))
  phi0 = sqrt(sum((eta_tilde - eta_dir)^2)/D)
  
  fit_mm = mm(y=y, X=X, beta0=beta0, phi0=phi0, nu=nu, add_std_error=FALSE,
              method='Newton', global='dbldog', xscalm="auto", control=list(xtol=tol, ftol=tol, maxit=maxiter))
  
  beta_mm[k, ] = fit_mm$coefficients
  phi_mm[k] = fit_mm$random_effects
  
  # PQL
  v0 = rep(0, times=D)
  xi0 = c(beta0, v0)
  fit_pql = pql(xi0, X, phi0, nu, y, S, epsXiPhi=1e-2,
                method='Newton', global='dbldog', xscalm="auto", control=list(xtol=tol, ftol=tol, maxit=maxiter))
  beta_pql[k, ] = fit_pql$xi[1:p]
  phi_pql[k] = fit_pql$phi
  
  tit = proc.time() - tit
  cat('\n', paste0(rep('-', 25), collapse=""))
  cat(paste('\n Iteracion K=', k, sep='', '\n'))
  print(tit)
}


# accuracy
theta_la = cbind(beta_la, phi_la)
theta_mm = cbind(beta_mm, phi_mm)
theta_pql = cbind(beta_pql, phi_pql)

# bias
bias_la = bias(x_est=theta_la, x_true=theta)
bias_mm = bias(x_est=theta_mm, x_true=theta)
bias_pql = bias(x_est=theta_pql, x_true=theta)

# rbias
rbias_la = rbias(x_est=theta_la, x_true=theta)
rbias_mm = rbias(x_est=theta_mm, x_true=theta)
rbias_pql = rbias(x_est=theta_pql, x_true=theta)

# mse
mse_la = mse(x_est=theta_la, x_true=theta)
mse_mm = mse(x_est=theta_mm, x_true=theta)
mse_pql = mse(x_est=theta_pql, x_true=theta)

# rmse
rmse_la = rmse(x_est=theta_la, x_true=theta)
rmse_mm = rmse(x_est=theta_mm, x_true=theta)
rmse_pql = rmse(x_est=theta_pql, x_true=theta)

# rrmse
rrmse_la = rrmse(x_est=theta_la, x_true=theta)
rrmse_mm = rrmse(x_est=theta_mm, x_true=theta)
rrmse_pql = rrmse(x_est=theta_pql, x_true=theta)

