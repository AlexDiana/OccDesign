
# GP FUNS ------

cov_fun <- function(dist, rho, sigma){
  
  sigma^2 * exp(- dist / rho)
  
}

dist_x <- function(x, y){
  
  sum((x - y)^2)
  
}

Sigma_X <- function(x, y, rho, sigma){
  
  n1 <- dim(x)[1]
  n2 <- dim(y)[1]
  Sigma <- matrix(NA, n1, n2)
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      dist <- dist_x(x[i,,], y[j,,])
      Sigma[i,j] <- cov_fun(dist, rho, sigma)
    }
  }
  
  Sigma
}

SigmaxX <- function(x, X_all, sigma, rho){
  
  dists <- apply(X_all, 1, function(y){
    dist_x(x, y)
  })
  
  cov_fun(dists, rho, sigma)
}

post_meansd <- function(X_all, y_all, rho, sigma){
  
  Sigma_XX <- Sigma_X(X_all, X_all, rho, sigma)
  invSigma_xx <- solve(Sigma_XX)
  
  function(x){
    
    Sigma_x <- SigmaxX(x, X_all, sigma, rho) # x, X_all
    
    mu_x <- Sigma_x %*% invSigma_xx %*% y_all
    sd_x <- sigma - Sigma_x %*% invSigma_xx %*% Sigma_x
    
    list("mu" = mu_x,
         "sd_x" = sd_x)
    
  }
  
}

delta_x <- function(mu_n_x, f_nstar){
  mu_n_x - f_nstar
}

pos_part <- function(x){
  x * (x > 0)
}

phi <- function(x){
  dnorm(x)
}

Phi <- function(x){
  pnorm(x)
}

EI <- function(mu_n_x, sigma_n_x, f_nstar){
  
  delta_nx <- delta_x(mu_n_x, f_nstar)
  
  pos_part(delta_nx) + 
    sigma_n_x * phi(delta_nx / sigma_n_x) - 
    abs(delta_nx) * Phi(delta_nx / sigma_n_x)
  
}
