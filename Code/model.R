library(ggplot2); library(here); library(tidyverse)
library(foreach); library(doParallel)
source(here("R","utility.R"))
source(here("R","fun.R"))
numCores <- detectCores()
registerDoParallel(numCores)

{
  library(OccDesign)
  
  # library(Rcpp); library(RcppArmadillo)
  # sourceCpp("src/code.cpp")
}

# COMPUTE OPTIMAL DESIGN -------

# true params
{
  # ppsi0_true <- c(.5, .25)
  # beta_psi_true <- c(1, -1)
  # beta_p_true <- c(-1)
  # 
  # coeffs_true <- list("psi" = c(logit(ppsi0_true[1]),beta_psi_true),
  #                     "p" = c(logit(ppsi0_true[2]),beta_p_true))
  
  beta_psi_true <- seq(0, 0, length.out = 4)
  ppsi0_true <- .5
  beta_p_true <- c(0)
  coeffs_true <- list("psi" = beta_psi_true,
                      "p" = c(logit(ppsi0_true),beta_p_true))
}

# create gradients
{
  set.seed(1)
  gradient1_psi <- createLinearGradient()
  # gradient1_psi <- createGradient()
  plotgradient(gradient1_psi)
  
  gradPsi_minmax <- findGradientMax(gradient1_psi)
  gradPsi_min <- gradPsi_minmax[1]
  gradPsi_max <- gradPsi_minmax[2]
  
  X_psibreaks <- c(gradPsi_min - .05, 
                   seq(gradPsi_min, gradPsi_max, length.out = 5)[-c(1,5)],
                   gradPsi_max + .05)  #c(-1,.25, .5, .75,1)
  
  set.seed(3)
  # gradient1_p <- createLinearGradient()
  gradient1_p <- createUniformGradient()
  plotgradient(gradient1_p)
  
  gradient_psi <- list(
    gradient1_psi
  )
  
  gradient_p <- list(
    gradient1_p
  )
  
  gradient_list <- list("psi" = gradient_psi,
                        "p" = gradient_p)
}

# number of locations
n <- 25

# total number of sampling occasions
n_occ <- 50

# generate sites locations
# X <- expand.grid(seq(0, 1, length.out = sqrt(n)),
#                  seq(0, 1, length.out = sqrt(n)))
# X <- cbind(seq(0, 1, length.out = n), 1)
X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
# X <- cbind(runif(n, 0, 1), runif(n, 0, 1))
# X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)

# starting values
{
  # theta <- rep(0, n - 1)
  theta <- seq(0, 0, length.out = n)
  
  alpha_k <- 2
  
  # auxiliary variables
  {
    maxM <- 15 # max M per site
    U_n1 <- matrix(runif(N_s * n / 2), N_s / 2, n)
    U_N1 <- matrix(runif(N_s * (maxM * n) / 2), N_s / 2, maxM * n)
    
    U_n <- rbind(U_n1, 1 - U_n1)
    U_N <- rbind(U_N1, 1 - U_N1)
  }
  
  # utility_val_current <- computeUtility(n, M, X, 
  #                                       gradient_list, coeffs_true,
  #                                       U_n, U_N)
}

N_x <- 100 # number of simulation replicates
N_s <- 100 # utlity estimation replicates

niter <- 50
utility_vals <- rep(NA, niter)
pi_vals <- matrix(NA, niter, n)

# variables
{
  u1 <- matrix(runif(n_occ * N_x / 2), N_x / 2, n_occ)
  u <- rbind(u1, 1 - u1)
}

for (i in 1:niter) {
  
  print(paste0("Iter = ",i))
  print(paste0("Theta = ",
               paste0(round(mapThetaToProb(theta), 3), collapse = "-")))
  
  # sampling values ----------
  
  pi <- mapThetaToProb(theta)
  x <- rmultinom_u(pi, u, n)
  # x <- generateSamples(pi, n_occ, N_x)
  
  H_x <- sapply(1:N_x, function(i){
    print(i)
    M <- x[,i]
    return(
      computeUtility(n, M, X, 
                     gradient_list, coeffs_true,
                     X_psibreaks,
                     U_n, U_N, N_s)
    )
  })
 
  # H_x <- foreach(i = 1:N_x,
  #                .combine = c,
  #                .packages = c("OccDesign")) %dopar% {
  #                  print(i)
  #                  M <- x[,i]
  #                  return( computeUtilityBase(n, M, X,
  #                                             gradient_list,
  #                                             coeffs_true,
  #                                             X_psibreaks,
  #                                             N_s))
  #                }
  
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

# PI DIAGNOSTICS -----

pi_vals_output <- pi_vals[1:(i-1),] %>% 
  as.data.frame() %>% 
  mutate(Iter = row_number()) %>% 
  pivot_longer(cols = starts_with("V"))

ggplot(data = pi_vals_output, 
       aes(x = Iter, 
           color = name,
           y = value)) + geom_line()

# OUTPUTS --------

data_plot <- pi_vals %>% as.data.frame() %>% 
  mutate(Niter = row_number()) %>% 
  pivot_longer(cols = contains("V"))

ggplot(data = data_plot[data_plot$name %in% c("V6","V23"),], 
       aes(x = Niter, color = name,
                             y = value)) + geom_line() + theme_bw()

qplot(1:niter, utility_vals)

# OUTPUT GRID ------

M_estimated <- n_occ * mapThetaToProb(theta)
# M_estimated <- round(M_estimated)

X_M <- data.frame(X, M_estimated)

plotgradient(gradient1_psi) + 
  geom_point(data = NULL, aes(
  x = X[,1], y = X[,2], size = M_estimated,
  color = "red")) + 
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  theme_bw() + xlab("") + 
  ylab("") 

computeUtilityBase(n, M, X,
                   gradient_list,
                   coeffs_true,
                   X_psibreaks,
                   N_s)

computeUtilityBase(n, M_estimated, X, 
                   gradient_list, 
                   coeffs_true,
                   X_psibreaks,
                   N_s)

M_base <- round(rep(n_occ / n, n))
M_base[1:2] <- M_base[1:2] + 1
computeUtilityBase(n, round(rep(n_occ / n, n)), 
                   X, 
                   gradient_list, 
                   coeffs_true,
                   X_psibreaks,
                   N_s)

# CHECK -----

computeUtilityBase(n, M, X, 
                   gradient_list, coeffs_true)

n_new <- 25
M_new <- rep(4, n_new)
X_new <- cbind(seq(0, 1, length.out = n_new), 1)

computeUtilityBase(n_new, M_new, X_new, 
                   gradient_list, coeffs_true)

computeUtility(n, M, X, 
               gradient_list, coeffs_true,
               U_n, U_N)

computeUtility(n_new, M_new, X_new, 
               gradient_list, coeffs_true,
               U_n, U_N)

#

X_new <- cbind(seq(0.2, .75, length.out = n), 1)

computeUtility(n_new, M_new, X_new, 
               gradient_list, coeffs_true,
               U_n, U_N)

X_new <- cbind(seq(0, .5, length.out = n), 1)

computeUtility(n_new, M_new, X_new, 
               gradient_list, coeffs_true,
               U_n, U_N)

X_new <- cbind(seq(.5, 1, length.out = n), 1)

computeUtility(n_new, M_new, X_new, 
               gradient_list, coeffs_true,
               U_n, U_N)

X_new <- cbind(seq(0, 1, length.out = n), 1)

computeUtility(n_new, M_new, X_new, 
               gradient_list, coeffs_true,
               U_n, U_N)

# EMPIRICAL CHECKS --------

n <- 20
X <- cbind(seq(0, 1, length.out = n),
           1)

# left skewed
{
  M <- rep(n_occ / n, n)
  M[1] <- M[1] + 2; M[n] <- M[n] - 2
  M[2] <- M[2] + 2; M[n - 1] <- M[n - 1] - 2
  M[3] <- M[3] + 2; M[n - 2] <- M[n - 2] - 2
  M[4] <- M[4] + 2; M[n - 3] <- M[n - 3] - 2
  
  computeUtilityBase(n, M, X, 
                     gradient_list, coeffs_true,
                     numSims = 10000)  
}

# right skewed
{
  M <- rep(n_occ / n, n)
  M[1] <- M[1] - 2; M[n] <- M[n] + 2
  M[2] <- M[2] - 2; M[n - 1] <- M[n - 1] + 2
  M[3] <- M[3] - 2; M[n - 2] <- M[n - 2] + 2
  M[4] <- M[4] - 2; M[n - 3] <- M[n - 3] + 2
  
  computeUtilityBase(n, M, X, 
                     gradient_list, coeffs_true,
                     numSims = 10000)  
}

# uniform
{
  M <- rep(n_occ / n, n)
  
  computeUtilityBase(n, M, X, 
                     gradient_list, coeffs_true,
                     numSims = 10000)
  
}


