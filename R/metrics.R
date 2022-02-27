
bias = function(x_est, x_true) {
  
  if (is.null(dim(x_true))) {
    b = colMeans(sweep(x_est, 2, x_true, FUN="-"))
  } else {
    b = colMeans(x_est - x_true)
  }
  
  return(b)
}


rbias = function(x_est, x_true) {
  
  if (is.null(dim(x_true))) {
    x_true = matrix(x_true, nrow=1)
    rb = colMeans(sweep(x_est, 2, x_true, FUN="-")) / colMeans(x_true)
  } else {
    rb = colMeans(x_est - x_true) / colMeans(x_true)
  }
  
  return(rb)
}


mse = function(x_est, x_true) {
  
  if (is.null(dim(x_true))) {
    e = colMeans((sweep(x_est, 2, x_true, FUN="-"))^2)
  } else {
    e = colMeans((x_est - x_true)^2)
  }
  return(e)
}


rmse = function(x_est, x_true) {
  
  if (is.null(dim(x_true))) {
    re = sqrt(colMeans((sweep(x_est, 2, x_true, FUN="-"))^2))
  } else {
    re = sqrt(colMeans((x_est - x_true)^2))
  }
  return(re)
}


rrmse = function(x_est, x_true) {
  
  if (is.null(dim(x_true))) {
    x_true = matrix(x_true, nrow=1)
    rre = sqrt(colMeans((sweep(x_est, 2, x_true, FUN="-"))^2)) / colMeans(x_true)
  } else {
    rre = sqrt(colMeans((x_est - x_true)^2)) / colMeans(x_true)
  }
  
  return(rre)
}
