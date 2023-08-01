loglik_vl <- function(mu){
  
  Xbeta <- X %*% beta
  
  - sum(log(computeloglik(B, y, Xbeta, C)))
  
}


computeUtility <- function(n, X,
                           gradient_list, coeffs,
                           X_psibreaks,
                           U_n, U_N, N_s){
  
  beta_psi <- coeffs[["psi"]]
  beta_p <- coeffs[["p"]]
  # rho <- coeffs[["rho"]]
  
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
  
  
estimateUtilityData <- function(y, M, X_psi, X_p,
                                trueParams){
    
    ncov_psi <- ncol(X_psi)
    ncov_p <- ncol(X_p)
    
    startingParams <- rep(0, ncov_psi + ncov_p)
    # startingParams <- trueParams
    
    fit_optim <- optim(
      par = startingParams,
      method = c("BFGS"),
      fn = loglik_vl(y, M, X_psi, X_p),
      gr = gr_loglik_cpp2(y, M, X_psi, X_p),
      hessian = T
    )
    
    Sigma <- solve(fit_optim$hessian)
    
    list(
      "mu" = fit_optim$par[1:ncov_psi],
      'Sigma' = Sigma[1:ncov_psi, 1:ncov_psi]
    )
  }

