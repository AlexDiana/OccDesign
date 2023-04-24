library(ggplot2); library(here); library(tidyverse)
library(foreach); library(doParallel)
source(here("R","utility.R"))
numCores <- detectCores()
registerDoParallel(numCores)

{
  # library(OccDesign)
  
  library(Rcpp); library(RcppArmadillo)
  sourceCpp("src/code.cpp")
  source("~/OccDesign/R/utility.R", echo=TRUE)
}

# COMPUTE OPTIMAL DESIGN -------

# true params
{
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

N_x <- 200 # number of simulation replicates
N_s <- 100 # utlity estimation replicates

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

niter <- 50
utility_vals <- rep(NA, niter)
pi_vals <- matrix(NA, niter, n)

# variables
{
  u1 <- matrix(runif(n_occ * N_x / 2), N_x / 2, n_occ)
  u <- rbind(u1, 1 - u1)
}

#

pi <- mapThetaToProb(theta)
M <- rmultinom(1, n_occ, pi)

M[1] <- M[1] - 1
M[2] <- M[2] + 1

computeUtility(n, M, X, 
               gradient_list, coeffs,
               U_n, U_N)

coeffs <- coeffs_true

coeff_psi <- coeffs[["psi"]]
coeff_p <- coeffs[["p"]]

X_psip <- computeXgradient(gradient_list, X, M, X_psibreaks)
X_psi <- X_psip$X_psi
X_p <- X_psip$X_p

psi <- psiModel(coeff_psi, X_psi)
p <- pModel(coeff_p, X_p)

trueParams <- c(coeff_psi, coeff_p)

idxsubsetUN <- function(M, n, maxM){
  unlist(sapply(1:n, function(i){
    (i - 1) * maxM + seq_len(M[i])
  }))
}

idxUn <- idxsubsetUN(M, n, maxM)

U_N_subset <- U_N[,idxUn]

utility_values2 <- t(sapply(1:N_s, function(i){
  
  y <- simDataEps(M, psi, p, U_n[i,], U_N_subset[i,])
  # y <- simData(M, psi, p)
  
  list_score <- estimateUtilityData(y, M, 
                                    X_psi, X_p,
                                    trueParams)
  
  list_score$psi
  
}))

mean(utility_values)

qplot(utility_values[1,])

qplot(1:N_s, cumsum(utility_values1 - utility_values2) / 1:N_s)
ggplot() + 
  geom_line(data = NULL, aes(x = 1:N_s, y = cumsum(utility_values1) / 1:N_s)) + 
  geom_line(data = NULL, aes(x = 1:N_s, y = cumsum(utility_values2) / 1:N_s),
            color = "red")


# EXPERIMENTS 2 ---------

x <- rmultinom_u(pi, u, n)

H_x2 <- sapply(1:N_x, function(i){
  print(i)
  M <- x[,i]
  return(
    computeUtility(n, M, X, maxM,
                   gradient_list, coeffs_true,
                   X_psibreaks,
                   U_n, U_N, N_s)
  )
})

{
  qplot(1:N_x, cumsum(H_x1 - H_x2) / 1:N_x)
  ggplot() + 
    geom_line(data = NULL, aes(x = 1:N_x, y = cumsum(H_x1) / 1:N_x)) + 
    geom_line(data = NULL, aes(x = 1:N_x, y = cumsum(H_x2) / 1:N_x),
              color = "red")
}

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

