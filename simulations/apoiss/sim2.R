# check the behaviour of the EBP
# model: log(\mu_{d}) = \beta\xx_{d} + \phi * v_{d}

# libraries
library(yaml)
library(lme4)
library(nleqslv)


# source code
source('R/metrics.R')
source('R/apoiss-MM-fitting-algorithm.R')
source('R/apoiss-PQL-fitting-algorithm.R')
source('R/apoiss-predictors.R')
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

pd = matrix(0, nrow=K, ncol=D)

ebp_la = matrix(0, nrow=K, ncol=D)
ebp_mm = matrix(0, nrow=K, ncol=D)

plug_la = matrix(0, nrow=K, ncol=D)
plug_mm = matrix(0, nrow=K, ncol=D)
plug_pql = matrix(0, nrow=K, ncol=D)


for(k in 1:K) {

  tit = proc.time()
  
  # response
  y_res = ysim(X, beta, phi, nu)
  pd[k, ] = y_res$pd
  y = y_res$y 
  
  # laplace
  fit_la = glmer(formula=y~x1+(1|dd), family=poisson, nAGQ=1,offset=log(nu))
  beta_la[k, ] = fit_la@beta
  phi_la[k] = fit_la@theta
  
  # MM
  beta0 = glm(formula=y~x1, family=poisson, offset=log(nu))$coefficients
  eta_tilde = X%*%beta0
  eta_dir = log((y + 1)/(nu + 1))
  phi0 = sqrt(sum((eta_tilde - eta_dir)^2)/D)
  
  fit_mm = MM(y=y, X=X, beta0=beta0, phi0=phi0, nu=nu, maxiter=maxiter, tol=tol, B=B, add_std_error=FALSE,
              method='Newton', global='dbldog', xscalm="auto", control=list(xtol=tol, ftol=tol, maxit=maxiter))
  
  beta_mm[k, ] = fit_mm$coefficients
  phi_mm[k] = fit_mm$random_effects
  
  # PQL
  v0 = rep(0, times=D)
  xi0 = c(beta0, v0)
  fit_pql = PQL(xi0, X, phi0, nu, y, S, tol, maxiter, epsXiPhi=1e-2,
                method='Newton', global='dbldog', xscalm="auto", control=list(xtol=tol, ftol=tol, maxit=maxiter))
  beta_pql[k, ] = fit_pql$xi[1:p]
  phi_pql[k] = fit_pql$phi
  
  # plug-in estimator
  ebpvd_la = ebp(y, X, beta_la[k, ], phi_la[k], nu, L, parameter='random_effects')
  ebpvd_mm = ebp(y, X, beta_mm[k, ], phi_mm[k], nu, L, parameter='random_effects')
  
  plug_la[k, ] = plugin(X, beta_la[k, ], phi_la[k], ebpvd_la)
  plug_mm[k, ] = plugin(X, beta_mm[k, ], phi_mm[k], ebpvd_mm)
  plug_pql[k, ] = plugin(X, beta_pql[k, ], phi_pql[k], fit_pql$xi[(p + 1) : (p + D)])
  
  # EBP
  ebp_la[k, ] = ebp(y, X, beta_la[k, ], phi_la[k], nu, L, parameter='target') 
  ebp_mm[k, ] = ebp(y, X, beta_mm[k, ], phi_mm[k], nu, L, parameter='target')
  
  tit = proc.time() - tit
  cat('\n', paste0(rep('-', 25), collapse=""))
  cat(paste('\n Iteracion K=', k, sep='', '\n'))
  print(tit)
}


# accuracy
avg_pd = apply(pd, 2, mean)
avg_ebp_la = apply(ebp_la, 2, mean)
avg_abp_mm = apply(ebp_mm, 2, mean)
avg_plug_la = apply(plug_la, 2, mean)
avg_plug_mm = apply(plug_mm, 2, mean)
avg_plug_pql = apply(plug_pql, 2, mean)

# RB_d
rbd_ebp_la = rbias(x_est=ebp_la, x_true=pd)
rbd_ebp_mm = rbias(x_est=ebp_mm, x_true=pd)

rbd_plug_la = rbias(x_est=plug_la, x_true=pd)
rbd_plug_mm = rbias(x_est=plug_mm, x_true=pd)
rbd_plug_pql = rbias(x_est=plug_pql, x_true=pd)

# RRE_d
rred_ebp_la = rrmse(x_est=ebp_la, x_true=pd)
rred_ebp_mm = rrmse(x_est=ebp_mm, x_true=pd)

rred_plug_la = rrmse(x_est=plug_la, x_true=pd)
rred_plug_mm = rrmse(x_est=plug_mm, x_true=pd)
rred_plug_pql = rrmse(x_est=plug_pql, x_true=pd)

# B_d
bd_ebp_la = bias(x_est=ebp_la, x_true=pd)
bd_ebp_mm = bias(x_est=ebp_mm, x_true=pd)

bd_plug_la = bias(x_est=plug_la, x_true=pd)
bd_plug_mm = bias(x_est=plug_mm, x_true=pd)
bd_plug_pql = bias(x_est=plug_pql, x_true=pd)

# E_d
ed_ebp_la = mse(x_est=ebp_la, x_true=pd)
ed_ebp_mm = mse(x_est=ebp_mm, x_true=pd)

ed_plug_la = mse(x_est=plug_la, x_true=pd)
ed_plug_mm = mse(x_est=plug_mm, x_true=pd)
ed_plug_pql = mse(x_est=plug_pql, x_true=pd)

# B
b_ebp_la = mean(abs(bd_ebp_la))
b_ebp_mm = mean(abs(bd_ebp_mm))

b_plug_la = mean(abs(bd_plug_la))
b_plug_mm = mean(abs(bd_plug_mm))
b_plug_pql = mean(abs(bd_plug_pql))

# RB
rb_ebp_la = mean(abs(rbd_ebp_la))
rb_ebp_mm = mean(abs(rbd_ebp_mm))

rb_plug_la = mean(abs(rbd_plug_la))
rb_plug_mm = mean(abs(rbd_plug_mm))
rb_plug_pql = mean(abs(rbd_plug_pql))

# E
e_ebp_la = mean(ed_ebp_la)
e_ebp_mm = mean(ed_ebp_mm)

e_plug_la = mean(ed_plug_la)
e_plug_mm = mean(ed_plug_mm)
e_plug_pql = mean(ed_plug_pql)

# RRE
rre_ebp_la = mean(rred_ebp_la)
rre_ebp_mm = mean(rred_ebp_mm)

rre_plug_la = mean(rred_plug_la)
rre_plug_mm = mean(rred_plug_mm)
rre_plug_pql = mean(rred_plug_pql)
