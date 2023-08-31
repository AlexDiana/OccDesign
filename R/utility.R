

# logistic function
logistic <- function(x){
  1 / (1 + exp(-x))
}

# logit function
logit <- function(x){
  log(x / (1 - x))
}

# simulate a bernoulli random variable using a uniform random variable
generateBernoulliEps <- function(U,  # uniform random variable
                                 p # probability
                                 ){
  as.numeric(U < p)
}

# simulate occupancy datasets using uniform random variables
simDataEps <- function(M,  # number of sample per sites
                       psi, # occupancy probability
                       p,  # detection probability
                       U_n, U_N # auxiliary variables
                       ){
  
  # simulate occupancies
  z <- generateBernoulliEps(U_n, psi)
  
  p_all <- p * rep(z, M)
  
  # simulate detections
  y <- generateBernoulliEps(U_N, p_all)
  
  y
}

rmultinom_singleu <- function(pi, u, n){
  
  vals <- as.numeric(cut(u, breaks = c(0,cumsum(pi))))
  tablevals <- table(vals)  
  samples <- rep(0, n)
  samples[as.numeric(names(tablevals))] <- tablevals
  
  samples
}

# simulate from a multinomial distribution using precomputed uniform variables 
rmultinom_u <- function(pi, # probability
                        u,  # uniform random variables
                        n # sample size
                        ){
  
  samples <- apply(u, 1, function(x){
    rmultinom_singleu(pi, x, n)
  })
  
}


estimateMSE <- function(mu_hat, mu0, Sigma){
  (mu_hat - mu0) %*% (mu_hat - mu0) + sum(diag(Sigma))
}

# compute the maximum posterior covariance
maxVar <- function(Sigma){
  max(diag(as.matrix(Sigma)))
}

# wrapper around the c++ function to compute loglikelihood
# of an occupancy model
loglik_r <- function(y, M, X_psi, X_p){
  
  n <- length(M)
  occ <- as.numeric(sapply(1:n, function(i){
    any(y[1:M[i] + sum(M[seq_len(i-1)])] > 0)
  }))
  occ_all <- rep(occ, M)
  y_notoccupied <- y[occ_all == 0]
  y_occupied <- y[occ_all == 1]
  sumocc <- sum(y_occupied)
  sumMocc <- sum(M[occ == 1])
  Mocc0 <- M[occ == 0]
  
  ncov_psi <- ncol(X_psi)
  ncov_p <- ncol(X_p)
  
  cov_psi <- 1:ncov_psi - 1
  cov_p <- ncov_psi + 1:ncov_p - 1
  
  sumM <- c(0,cumsum(M)[-n])
  
  loglik_data <- function(pars){
    
    # print(pars)
    
    loglik_cpp(pars, y, M, X_psi, X_p, cov_psi, cov_p, 
               y_occupied, occ, occ_all, sumM)
    
  }
  
  loglik_data
  
}

# wrapper around the c++ function to compute gradient of loglikelihood
# of an occupancy model
gr_loglik_r <- function(y, M, X_psi, X_p){
  
  n <- length(M)
  occ <- as.numeric(sapply(1:n, function(i){
    any(y[1:M[i] + sum(M[seq_len(i-1)])] > 0)
  }))
  occ_all <- rep(occ, M)
  y_notoccupied <- y[occ_all == 0]
  y_occupied <- y[occ_all == 1]
  sumocc <- sum(y_occupied)
  sumMocc <- sum(M[occ == 1])
  Mocc0 <- M[occ == 0]
  
  ncov_psi <- ncol(X_psi)
  ncov_p <- ncol(X_p)
  
  cov_psi <- 1:ncov_psi - 1
  cov_p <- ncov_psi + 1:ncov_p - 1
  
  sumM <- c(0,cumsum(M)[-n])
  
  loglik_data <- function(pars){
    
    gr_loglik_cpp(pars, y, M, X_psi, X_p, cov_psi, cov_p,
                  y_occupied, occ, occ_all, sumM)
    
  }
  
  loglik_data
  
}

# find utility function for a given dataset by computing a Laplace approximation
# to the posterior
estimateUtilityData <- function(y,  # data
                                M,  # number of sample per site
                                X_psi, # occupancy covariates matrix
                                X_p, # detection covariates matrix
                                trueParams # true value of the parameters
                                ){
  
  ncov_psi <- ncol(X_psi)
  ncov_p <- ncol(X_p)
  
  startingParams <- rep(0, ncov_psi + ncov_p)
  # startingParams <- trueParams
  
  # maximize the posterior and compute Hessian
  fit_optim <- optim(
    par = startingParams,
    method = c("BFGS"),
    fn = loglik_r(y, M, X_psi, X_p),
    gr = gr_loglik_r(y, M, X_psi, X_p),
    hessian = T
  )
  
  # posterior covariance matrix
  Sigma <- solve(fit_optim$hessian)
  
  list(
    "mu" = fit_optim$par[1:ncov_psi],
    'Sigma' = Sigma[1:ncov_psi, 1:ncov_psi]
  )
}

# convert the parameters gamma to the vector of probabilities pi
mapGammaToPi <- function(gamma){
  
  gamma <- exp(gamma)
  
  gamma / sum(gamma)
  
}

# create the subset of detection auxiliary uniform variables used
idxsubsetUN <- function(M, n, maxM){
  unlist(sapply(1:n, function(i){
    (i - 1) * maxM + seq_len(M[i])
  }))
}

# create design matrix for detection
createDM_P <- function(X){
  cbind(1, X)
}

# compute the utility function 
computeUtility <- function(n,  # number of sites
                           M,  # number of samples per site
                           X_psi, # occupancy covariates matrix
                           X_p, # detection covariates matrix 
                           maxM, # used for assigning maximum dimensions
                           coeffs, # value of the true coefficients
                           U_n, U_N # auxiliary variables to simulate datasets
                           ){
  
  # import true parameters
  beta_psi <- coeffs[["psi"]]
  beta_p <- coeffs[["p"]]
  
  trueParams <- c(beta_psi, beta_p)
  
  # compute true occupancy and detection
  X_p <- X_p[rep(1:n, M),]
  X_p <- createDM_P(X_p)
  
  psi <- logistic(X_psi %*% beta_psi)
  p <- logistic(X_p %*% beta_p)
  
  # subset the detection probability auxiliary variables
  idxUn <- idxsubsetUN(M, n, maxM)
  U_N_subset <- U_N[,idxUn]
  
  # number of simulations
  N_s <- nrow(U_n)
  
  utility_values <- t(sapply(1:N_s, function(i){
    
    # simulate occupancy dateset
    y <- simDataEps(M, psi, p, U_n[i,], U_N_subset[i,])
    
    # estimate utility
    list_score <- estimateUtilityData(y, M, X_psi, X_p,
                                      trueParams)
    mu_hat <- list_score$mu
    Sigma_hat <- list_score$Sigma
    
    maxVar(Sigma_hat)
    
  }))
  
  mean(utility_values)
  
}

# compute gradient of log(f|M)
grad_logfX <- function(x, theta, n_occ){
  x - n_occ * exp(theta) / sum(exp(theta))
}

# compute the gradient using control variates as in formula (4) of the Appendix
computeReducedVarianceGradient <- function(
    H_x, # value of the utility function 
    M_all, # value of the number of samples per site of each iteration
    n_occ, # total number of samples
    theta # current parameter value
) {
  
  N_x <- length(H_x)
  n <- nrow(M_all)
  
  # compute term (1) of formula (4)
  term1 <- sapply(1:N_x, function(i){
    M_current <- M_all[,i]
    H_x[i] * grad_logfX(M_current, theta, n_occ)
  })
  
  # compute term (2) of formula (4)
  term2 <- sapply(1:N_x, function(i){
    M_current <- M_all[,i]
    grad_logfX(M_current, theta, n_occ)
  })
  
  var_h <- apply(term2, 1, var)
  cov_h <- sapply(1:n, function(i){
    cov(term1[i,], term2[i,])
  })
  
  a_i <- cov_h / var_h
  
  # compute the new gradient estimate using formula (4)
  gradient_estimates <- sapply(1:N_x, function(i){
    M_current <- M_all[,i]
    H_x[i] * grad_logfX(M_current, theta, n_occ) -
      grad_logfX(M_current, theta, n_occ) * a_i
  })
  
  apply(gradient_estimates, 1, mean)
}

# compute the utility function without using auxliary variables
computeUtilityBase <- function(n,  # number of sites
                                    M,  # number of samples per site
                                    X_psi, # occupancy covariates matrix
                                    X_p, # detection covariates matrix 
                                    coeffs, # value of the true coefficients
                                    numSims){
  
  beta_psi <- coeffs[["psi"]]
  beta_p <- coeffs[["p"]]
  
  X_p <- X_p[rep(1:n, M),]
  X_p <- createDM_P(X_p)
  
  psi <- logistic(X_psi %*% beta_psi)
  p <- logistic(X_p %*% beta_p)
  
  trueParams <- c(beta_psi, beta_p)
  
  utility_values <- t(sapply(1:numSims, function(i){
    
    y <- simData(M, psi, p)
    
    list_score <- estimateUtilityData(y, M, X_psi, X_p,
                                      trueParams)
    mu_hat <- list_score$mu
    Sigma_hat <- list_score$Sigma
    
    maxVar(Sigma_hat)
    
  }))
  
  mean(utility_values)
  
}
