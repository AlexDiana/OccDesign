findOptimalDesign <- function(niter, 
                              n, n_occ,
                              N_x, N_s,
                              coeffs_true,
                              X, X_psibreaks,
                              gradient_list){
  
  utility_vals <- rep(NA, niter)
  pi_vals <- matrix(NA, niter, n)
  
  # variables
  {
    u1 <- matrix(runif(n_occ * N_x / 2), N_x / 2, n_occ)
    u <- rbind(u1, 1 - u1)
  }
  
  # starting values
  {
    # theta <- rep(0, n - 1)
    theta <- seq(0, 0, length.out = n)
    
    alpha_k <- .5
    
    # auxiliary variables
    {
      maxM <- (n_occ / n) * 5 # max M per site
      U_n1 <- matrix(runif(N_s * n / 2), N_s / 2, n)
      U_N1 <- matrix(runif(N_s * (maxM * n) / 2), N_s / 2, maxM * n)
      
      U_n <- rbind(U_n1, 1 - U_n1)
      U_N <- rbind(U_N1, 1 - U_N1)
    }
    
    # utility_val_current <- computeUtility(n, M, X, 
    #                                       gradient_list, coeffs_true,
    #                                       U_n, U_N)
  }
  
  for (i in 1:niter) {
    
    print(paste0("Iter = ",i))
    print(paste0("Theta = ",
                 paste0(round(mapThetaToProb(theta), 3), collapse = "-")))
    
    # sampling values ----------
    
    pi <- mapThetaToProb(theta)
    x <- rmultinom_u(pi, u, n)
    # x <- generateSamples(pi, n_occ, N_x)
    
    # H_x <- sapply(1:N_x, function(i){
    #   print(i)
    #   M <- x[,i]
    #   return(
    #     computeUtility(n, M, X, maxM,
    #                    gradient_list,
    #                    coeffs_true,
    #                    X_psibreaks,
    #                    U_n, U_N, N_s)
    #   )
    # })
    
    H_x <- foreach(i = 1:N_x,
                   .combine = c,
                   .packages = c("OccDesign")) %dopar% {
                     # print(i)
                     M <- x[,i]
                     return( computeUtility(n, M, X, maxM,
                                            gradient_list,
                                            coeffs_true,
                                            X_psibreaks,
                                            U_n, U_N, N_s))
                   }
    
    # newton update
    {
      # Stheta_x <- S_theta(H_x)
      # 
      # grad_est <- sapply(1:N_x, function(i){
      #   x_current <- x[,i]
      #   Stheta_x[i] * grad_logfX(x_current, theta, n_occ)
      # })
      # 
      # grad_mean <- apply(grad_est, 1, mean)
      # 
      # Hessian_list <- lapply(1:N_x, function(i){
      #   x_current <- x[,i]
      #   Stheta_x[i] * (H_logfX(theta, n_occ) + grad_logfX(x_current, theta, n_occ) %*%
      #                    t(grad_logfX(x_current, theta, n_occ)))
      # })
      # 
      # Hessian_mean <- Reduce("+", Hessian_list) / N_x
      # 
      # qplot(1:n, grad_mean)
      # solve(Hessian_mean) %*% grad_mean3
    }
    
    # gradient update
    {
      Stheta_x <- S_theta(H_x)
      
      grad_est <- sapply(1:N_x, function(i){
        x_current <- x[,i]
        Stheta_x[i] * grad_logfX(x_current, theta, n_occ)
      })
      
      grad_est2 <- sapply(1:N_x, function(i){
        x_current <- x[,i]
        grad_logfX(x_current, theta, n_occ)
      })
      
      var_h <- apply(grad_est2, 1, var)
      cov_h <- sapply(1:n, function(i){
        cov(grad_est[i,], grad_est2[i,])
      })
      
      a_i <- cov_h / var_h
      
      grad_est3 <- sapply(1:N_x, function(i){
        x_current <- x[,i]
        Stheta_x[i] * grad_logfX(x_current, theta, n_occ) -
          grad_logfX(x_current, theta, n_occ) * a_i
      })
      
      grad_mean3 <- apply(grad_est3, 1, mean)
      qplot(1:n, grad_mean3)
      theta <- theta - alpha_k * grad_mean3
    }
    
    # old code
    {
      # Stheta_x <- S_theta(H_x)
      # 
      # w <- Stheta_x / sum(Stheta_x)
      # 
      # Ep_X <- apply(
      #   apply(x, 1, function(x_i){ x_i * w}),
      #   2, sum)
      # 
      # Tx <- x
      # 
      # Etheta_X <- apply(x, 1, mean)
      # 
      # Tx_sum <- apply(Tx, 1, sum)
      # 
      # Vartheta_X <- 1 / (N_x - 1) * (Tx %*% t(Tx)) - 
      #   1 / (N_x^2 - N_x) * Tx_sum %*% t(Tx_sum)
      # 
      # theta <- theta + alpha_k * solve(Vartheta_X) %*% (Ep_X - Etheta_X)
    }
    
    utility_vals[i] <- mean(H_x)
    pi_vals[i,] <- pi
  }
  
  list("utility_vals" = utility_vals,
       "pi_vals" = pi_vals)
}



plotDiagnostics <- function(pi_vals){
  
  pi_vals_output <- pi_vals %>% 
    as.data.frame() %>% 
    mutate(Iter = row_number()) %>% 
    pivot_longer(cols = starts_with("V"))
  
  ggplot(data = pi_vals_output, 
         aes(x = Iter, 
             color = name,
             y = value)) + geom_line()
  
  
}

plotObjective <- function(utility_vals){
  
  ggplot(data = NULL, aes(x = seq_along(utility_vals),
                          y = utility_vals)) + 
    geom_line() + geom_point() + 
    theme_bw() + ylab("Objective function") + xlab("Iteration")
  
}

plotDesign <- function(pi_vals, n_occ, X,
                       gradient){
  
  M_estimated <- n_occ * pi_vals[nrow(pi_vals),]
  # M_estimated <- round(M_estimated)
  
  X_M <- data.frame(X, M = M_estimated)
  
  plotgradient(gradient) + 
    geom_point(data = X_M, aes(
      x = X1, y = X2, size = M),
      color = "red") + 
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    theme_bw() + xlab("") + 
    ylab("") + guides(
      fill = F
    )
  
}
