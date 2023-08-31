

T_thetaX <- function(theta, X, n){
  
  apply(X, 2, function(x){
    grad_logfX(x, theta, n)
    # x[-n] - n * exp(theta) / (1 + sum(exp(theta)))  
  })
  
  
}

H_logfX <- function(theta, n_occ){
  
  n <- length(theta)
  
  H_f <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if(i == j){
        H_f[i,j] <- - n_occ * exp(theta[i]) / sum(exp(theta)) + 
          n_occ * exp(theta[i])^2 / sum(exp(theta))^2 
      } else {
        H_f[i,j] <- n_occ * exp(theta[i]) * exp(theta[j]) / (sum(exp(theta)))^2
      }
    }
  }
  
  return(H_f)
  
}




S_theta <- function(hx){
  # (hx - hlb) / (1 + exp(- S0 * (hx - gammatheta)))
  # exp(hx)
  hx
}

computeWeights <- function(x, H){
  
  N_x <- ncol(x)
  
  H_x <- sapply(1:N_x, function(i){
    
    M <- x[,i]
    
    - computeUtilityBase(n, M, X, 
                         gradient_list, coeffs_true)
    
    
  })
  
  Stheta_x <- S_theta(H_x)
  
  w <- Stheta_x / sum(Stheta_x)
  
  Ep_X <- apply(
    apply(x, 1, function(x_i){ x_i * w}),
    2, sum)
  
  Tx <- x
  
  Etheta_X <- apply(x, 1, mean)
  
  Tx_sum <- apply(Tx, 1, sum)
  
  Vartheta_X <- 1 / (N_x - 1) * (Tx %*% t(Tx)) - 
    1 / (N_x^2 - N_x) * Tx_sum %*% t(Tx_sum)
}



generateSamples <- function(theta, N, n){
  
  rmultinom(n, N, theta)
  
}

generateSamplesAntithetic <- function(theta, N, n){
  
  rmultinom(n, N, theta)
  
}


generateNewUN <- function(M1){
  
  U1 <- runif(M1 / 2)
  U2 <- 1 - U1
  
  c(U1, U2)
}

shiftSamplOcc <- function(n, M, X,
                          U_n, U_N){
  
  sumM <- c(0, cumsum(M)[-n])
  dists <- as.matrix(dist(X))
  
  i <- sample(1:n, 1)
  
  closestPoints <- which(abs(dists[i,-i] - min(dists[i,-i])) < 
                           .Machine$double.eps ^ 0.5)
  closestPoint <- sample(closestPoints, 1)
  
  # remove one
  idxOldObs <- sumM[i] + M[i]
  idxNewObs <- sumM[closestPoint] + (M[closestPoint] + 1)
  U_N[,idxOldObs] <- NA
  # insert row
  U_N_first <- U_N[,seq_len(idxNewObs - 1)]
  U_N_new <- matrix(generateNewUN(M1), M1, 1)
  if(idxNewObs <= ncol(U_N)){
    U_N_last <- U_N[,seq(from = idxNewObs, to = ncol(U_N))]
    U_N <- cbind(U_N_first, 
                 U_N_new, 
                 U_N_last)
  } else {
    U_N <- cbind(U_N_first, 
                 U_N_new)
  }
  
  U_N <- U_N[,!is.na(U_N[1,])]
  
  M[i] <- M[i] - 1
  M[closestPoint] <- M[closestPoint] + 1
  
  # resize
  if(M[i] == 0){
    M[i] <- M[n]
    M <- M[-n]
    X[i,] <- X[n,]
    X <- X[-n,]
    U_n[,i] <- U_n[,n]
    U_n <- U_n[,-n]
    n <- n - 1
  }
  
  list("n" = n,
       "M" = M,
       "X" = X,
       "U_n" = U_n,
       "U_N" = U_N)
}


plotgradient <- function(gradient0){
  
  x1 <- seq(0, 1, length.out = 100)
  x2 <- seq(0, 1, length.out = 100)
  loc <- expand.grid(x1, x2)
  
  colorfill <- apply(loc, 1, gradient0)
  
  ggplot() + 
    geom_tile(aes(x = loc[,1], y = loc[,2], fill = colorfill))
  
}


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

# gp funs ---


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

