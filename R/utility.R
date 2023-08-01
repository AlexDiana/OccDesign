

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
    
    maxVar(mu_hat, beta_psi, Sigma_hat)
    
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


# --------


computeM <- function(n, M_avg){
  M_min <- floor(M_avg)
  M_max <- M_min + 1
  n2 <- round((M_avg - M_min) * n)
  n1 <- n - n2
  # (n1 * floor(M) + n2 * (floor(M) + 1)) / (n1 + n2)
  # c(n1, n2)
  c(rep(M_min, each = n1), rep(M_max, n2))
}

simData <- function(M, psi, p){
  
  n <- length(M)
  
  z <- rbinom(n, 1, psi)
  p_all <- p * rep(z, M)
  y <- rbinom(sum(M), 1, p_all)
  
  y
}

loglik <- function(y, M, X_psi, X_p){
  
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
  
  loglik_data <- function(pars){
    
    # print(pars)
    
    beta_psi <- pars[1:ncov_psi]
    beta_p <- pars[ncov_psi + 1:ncov_p]
    
    psi <- logistic(X_psi %*% beta_psi)
    p <- logistic(X_p %*% beta_p)
    
    p_occupied <- p[occ_all == 1]
    psi_occ <- psi[occ == 0]
    
    # occupancies for occupied
    loglik_psi_occ <- 
      sum(log(psi[occ == 1]))
    
    # detections for occupied
    loglik_p_occ <- 
      sum(dbinom(y_occupied, 1, p_occupied, log = T))
    
    p_prod <- sapply(which(occ == 0), function(i){
      prod(1 - p[sum(M[seq_len(i-1)]) + 1:M[i]])
    })
    
    loglik_ppsi_nocc <- 
      sum(
        log(
          (1 - psi_occ) + psi_occ * p_prod
        )
      )
    
    logprior_psi <- dnorm(beta_psi[1], 0, 2, log = T)
    logprior_p <- dnorm(beta_p[1], 0, 2, log = T)
    
    logprior_betapsi <- sum(dnorm(beta_psi[-1], 0, 2, log = T))
    logprior_betap <- sum(dnorm(beta_p[-1], 0, 2, log = T))
    
    - (loglik_psi_occ + loglik_p_occ + loglik_ppsi_nocc +
         logprior_psi + logprior_p + logprior_betapsi + logprior_betap)
    
    # loglik_psi_occ + loglik_p_occ + loglik_ppsi_nocc
    # loglik_ppsi_nocc
  }
  
  loglik_data
  
}


gr_loglik_0 <- function(y, M, X_psi, X_p){
  
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
  
  loglik_data <- function(pars){
    
    # print(pars)
    
    beta_psi <- pars[1:ncov_psi]
    beta_p <- pars[ncov_psi + 1:ncov_p]
    
    psi <- logistic(X_psi %*% beta_psi)
    p <- logistic(X_p %*% beta_p)
    
    p_occupied <- p[occ_all == 1]
    psi_occ <- psi[occ == 0]
    
    grad_l <- rep(0, ncov_psi + ncov_p)
    
    # occupancies for occupied
    loglik_psi_occ <- 
      apply(sapply(which(occ == 1), function(i){
        as.vector(X_psi[i,]) * (exp(- X_psi[i,] %*% beta_psi) / 
                                  (1 + exp(- X_psi[i,] %*% beta_psi)))[1]
      }),1,sum)
    
    grad_l[1:ncov_psi] <- grad_l[1:ncov_psi] + loglik_psi_occ
    
    # detections for occupied
    loglik_p_occ <- 
      apply(sapply(which(occ_all == 1), function(i){
        as.vector(X_p[i,]) *
          ((y[i] * exp(- X_p[i,] %*% beta_p) - 
              (1 - y[i])) / (1 + exp(- X_p[i,] %*% beta_p)))[1]
      }), 1, sum)
    # sum( X_p[occ_all == 1,] *
    # (y_occupied * exp(- X_p[occ_all == 1,] %*% beta_p) + 
    # (1 - y_occupied)) / (1 + exp(- X_p[occ_all == 1,] %*% beta_p))) 
    
    grad_l[ncov_psi + 1:ncov_p] <- grad_l[ncov_psi + 1:ncov_p] + loglik_p_occ
    
    p_prod <- sapply(1:n, function(i){
      prod(1 - p[sum(M[seq_len(i-1)]) + 1:M[i]])
    })
    
    term1 <- 
      1 / ((1 - psi) + psi * p_prod )
    
    grad_betapsi <- 
      apply(
        sapply(which(occ == 0), function(i){
          # (1 / term1[i]) * (p_prod[i] - 1) *
          # as.vector(X_psi[i,]) * (exp(- X_psi[i,] %*% beta_psi) /
          #                           (1 + exp(- X_psi[i,] %*% beta_psi))^2)[1]
          as.vector(X_psi[i,]) *
            ( (p_prod[i] - 1)  /
                ((p_prod[i] + exp(- X_psi[i,] %*% beta_psi)) *
                   (1 + exp( X_psi[i,] %*% beta_psi))))[1]
        }), 1, sum)
    
    grad_betap <- 
      apply(
        sapply(which(occ == 0), function(i){
          
          - term1[i] * psi[i] *
            apply(sapply(1:M[i], function(k){
              
              idx <- sum(M[seq_len(i-1)]) + k
              currentProd <- prod(1 - p[setdiff(sum(M[seq_len(i-1)]) + 1:M[i], idx)])
              term2 <- ( (exp(- X_p[idx,] %*% beta_p)) /
                           (1 + exp(-X_p[idx,] %*% beta_p))^2)[1]
              # print(paste0("idx = ", idx , " - val = ", (exp(- X_p[idx,] %*% beta_p))))
              # print(paste0("idx = ", idx , " - val = ", (1 + exp(-X_p[idx,] %*% beta_p))^2))
              
              as.vector(X_p[idx,]) * currentProd  * term2
              
              # as.vector(X_p[idx,]) *
              #   ( currentProd  *
              #       ( (exp(- X_p[idx,] %*% beta_p)) /
              #          (1 + exp(-X_p[idx,] %*% beta_p))^2))[1]
              
            }),1,sum)
          
        }), 1, sum)
    
    # print(grad_beta_p)
    
    grad_l[1:ncov_psi] <- grad_l[1:ncov_psi] + grad_betapsi
    grad_l[ncov_psi + 1:ncov_p] <- grad_l[ncov_psi + 1:ncov_p] + grad_betap
    
    grad_l <- grad_l - 1 / (2 * 4) * pars
    
    # logprior_psi <- dnorm(beta_psi[1], 0, 2, log = T)
    # logprior_p <- dnorm(beta_p[1], 0, 2, log = T)
    # 
    # logprior_betapsi <- sum(dnorm(beta_psi[-1], 0, 2, log = T))
    # logprior_betap <- sum(dnorm(beta_p[-1], 0, 2, log = T))
    # 
    # - (loglik_psi_occ + loglik_p_occ + loglik_ppsi_nocc + 
    #      logprior_psi + logprior_p + logprior_betapsi + logprior_betap)
    
    - grad_l
  }
  
  loglik_data
  
}






psiModel <- function(coeff_psi, X){
  logistic(X %*% coeff_psi)
}

pModel <- function(coeff_p, X){
  logistic(X %*% coeff_p)
}

findGradientMax <- function(gradient, min = 0, max = 1){
  
  unifGrid <- expand.grid(seq(min, max, length.out = 100),
                          seq(min, max, length.out = 100))
  
  vals <- apply(unifGrid, 1, gradient)
  
  c(min(vals),
    max(vals))
  
}

computeUtilityBase <- function(n, M, X, 
                               gradient_list, coeffs,
                               X_psibreaks,
                               numSims){
  
  coeff_psi <- coeffs[["psi"]]
  coeff_p <- coeffs[["p"]]
  
  X_psip <- computeXgradient(gradient_list, X, M, X_psibreaks)
  X_psi <- X_psip$X_psi
  X_p <- X_psip$X_p
  
  psi <- psiModel(coeff_psi, X_psi)
  p <- pModel(coeff_p, X_p)
  
  trueParams <- c(coeff_psi, coeff_p)
  
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

computeUtilityBaseFinal <- function(n, M, 
                                    X_psi, X_p, maxM,
                                    coeffs,
                                    X_psibreaks,
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

computeUtilityBaseAnthitetic <- function(n, M, X, 
                                         gradient_list, coeffs,
                                         numSims){
  
  coeff_psi <- coeffs[["psi"]]
  coeff_p <- coeffs[["p"]]
  
  X_psip <- computeXgradient(gradient_list, X, M)
  X_psi <- X_psip$X_psi
  X_p <- X_psip$X_p
  
  psi <- psiModel(coeff_psi, X_psi)
  p <- pModel(coeff_p, X_p)
  
  trueParams <- c(coeff_psi, coeff_p)
  
  utility_values <- t(sapply(1:numSims, function(i){
    
    y <- simData(M, psi, p)
    
    list_score <- estimateUtilityData(y, M, 
                                      X_psi, X_p,
                                      trueParams)
    
    list_score$psi
    
  }))
  
  mean(utility_values)
  
}





computeXmatrix <- function(X_psi, 
                           M,
                           X_p){
  
}

computeUtility0 <- function(n, M, X, maxM,
                           gradient_list, coeffs,
                           X_psibreaks,
                           U_n, U_N, N_s){
  
  beta_psi <- coeffs[["psi"]]
  beta_p <- coeffs[["p"]]
  
  X_psip <- computeXgradient(gradient_list, X, M, X_psibreaks)
  X_psi <- X_psip$X_psi
  X_p <- X_psip$X_p
  
  psi <- logistic(X_psi %*% beta_psi)
  p <- logistic(X_p %*% beta_p)
  
  trueParams <- c(beta_psi, beta_p)
  
  idxUn <- idxsubsetUN(M, n, maxM)
  
  U_N_subset <- U_N[,idxUn]
  
  utility_values <- t(sapply(1:N_s, function(i){
    
    y <- simDataEps(M, psi, p, U_n[i,], U_N_subset[i,])
    
    list_score <- estimateUtilityData(y, M, X_psi, X_p,
                                      trueParams)
    mu_hat <- list_score$mu
    Sigma_hat <- list_score$Sigma
    
    maxVar(mu_hat, beta_psi, Sigma_hat)
    
    # score
    # c("var_psi" = list_score$psi,
    #   # "var_p" = list_score$p,
    #   # "var_2" = det(list_score$Sigma)
    #   "var_2" = sum(eigen(list_score$Sigma)$values)
    # )
    
    # estimateMSE
    
  }))
  
  # cor(utility_values)
  
  # util_vals_psi <- utility_values[,1]
  # util_vals_p <- utility_values[,2]
  # 
  # M_psi <- mean(util_vals_psi)
  # M_p <- mean(util_vals_p)
  # 
  # C_psip <- cov(utility_values)[1,2]
  # var_p <- var(util_vals_p)
  # 
  # beta_psip <- - C_psip / var_p
  # 
  # y_star <- util_vals_psi + beta_psip * (util_vals_p - M_p)
  
  # mean(y_star)
  
  # var(util_vals_psi)
  # var(y_star)
  
  # start_M <- 100
  # ggplot() + geom_line(data = NULL, 
  #                      aes(x = start_M:M1, 
  #                          y = (cumsum(util_vals_psi) / 1:M1)[start_M:M1])) +
  #   geom_line(data = NULL, 
  #             aes(x = start_M:M1, 
  #                 y = (cumsum(y_star) / 1:M1)[start_M:M1]), color = "red")
  
  # start_M <- 100
  # qplot(start_M:M1, (cumsum(y_star) / 1:M1)[start_M:M1])
  
  # i <- which.max(utils_vals)
  # y_current <- utility_values[[2 + (i-1)*2]]
  # M_current <- utility_values[[3 + (i-1)*3]]
  # 
  # optim(par = c(0, 0),
  #       fn = loglik(y_current, M_current),
  #       hessian = T)
  
  mean(utility_values)
  # mean(util_vals_psi)
  # - mean(y_star)
  
}

buildStartingX <- function(n){
  
  sqrt_n <- floor(sqrt(n))
  nmsqrt <- n - sqrt_n^2
  
  x_grid <- seq(0, 1, length.out = sqrt_n + 2)[-c(1,sqrt_n+2)]
  
  X <- as.matrix(expand.grid(x_grid, x_grid))
  
  X <- rbind(X, 
             cbind(runif(nmsqrt, 0, 1),
                   runif(nmsqrt, 0, 1))
  )
  
  X
}

createGradient <- function(){
  
  n_c <- rpois(1, lambda = 10)
  m_c <- cbind(
    runif(n_c, 0, 1),
    runif(n_c, 0, 1)
  )
  sd_c <- rgamma(n_c, 1, 3)
  
  gradient <- function(x){
    
    dens_x <- sum(
      sapply(1:n_c, function(k){
        dnorm(x[1], m_c[k,1], sd_c[1]) *   
          dnorm(x[2], m_c[k,2], sd_c[1])
      })
    )
    
    log(dens_x + 1)
  }  
}

createUniformGradient <- function(){
  
  function(x){
    0
  }
  
}

createLinearGradient <- function(){
  
  function(x){
    x[1]
  }
  
}


createDM_Psi <- function(X, breaks){ # c(-1,.25, .5, .75,1)
  # cbind(1, X, X^2)
  # X_psi %>% cut(., breaks = c(-1,.3,.6,2)) %>% 
  # X_psim <- as.character(cut(X, breaks = breaks))
  # as.matrix(model.matrix(~. - 1, as.data.frame(as.matrix(X_psim))))
  X_psim <- cut(X, breaks = breaks)
  as.matrix(model.matrix(
    ~ X_psim - 1,
    data = X_psim))
}

# X_psim <- cut(X, breaks = breaks)
# as.matrix(model.matrix(~ . - 1, X_psim))
# x <- model.matrix(
#    ~ X_psim - 1,
#   data = X_psim)
# as.matrix(model.matrix(~. - 1, as.data.frame(as.matrix(X_psim))))


computeXgradient <- function(gradient_list, X, M, breaks){
  
  n <- length(M)
  
  gradient_list_psi <- gradient_list[["psi"]]
  gradient_list_p <- gradient_list[["p"]]
  
  X_psi <- sapply(gradient_list_psi, function(gradient){
    apply(X, 1, gradient)
  })
  
  X_psi <- createDM_Psi(X_psi, breaks)
  
  X_p <- sapply(gradient_list_p, function(gradient){
    apply(X, 1, gradient)
  })
  
  X_p <- X_p[rep(1:n, M),]
  
  X_p <- createDM_P(X_p)
  
  list("X_psi" = X_psi,
       "X_p" = X_p)
}

moveLocation <- function(n, M, X, sd_x){
  
  i <- sample(1:n, 1)
  
  X[i,] <- X[i,] + rnorm(2, 0, sd_x)
  
  list("n" = n,
       "M" = M,
       "X" = X)
}

proposeLocations <- function(n, M, X, sd_x){
  
  X <- X + matrix(rnorm(n * 2, 0, sd_x), n, 2)
  
  list("n" = n,
       "M" = M,
       "X" = X)
}
