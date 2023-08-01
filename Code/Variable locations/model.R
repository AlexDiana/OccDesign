library(ggplot2); library(MASS)
library(Rcpp); library(RcppArmadillo)
library(here)
sourceCpp(here("Code","Variable Locations","ghk.cpp"))

n <- 100

# spatial
{
  sigma_gp <- 1
  l_gp <- .05
  tau <- .01
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

# simulate data
{
  ncov_psi <- 1; ncov_p <- 1
  beta_psi <- rnorm(ncov_psi)
  beta_p <- rnorm(ncov_p)
  
  X_psi <- matrix(rnorm(n * ncov_psi), n, ncov_psi)
  X_p <- matrix(rnorm(n * ncov_p), n, ncov_p)
}

# number of locations
n <- 100

# starting locations
{
  X_S <- cbind(runif(n), runif(n))
  X_S <- as.matrix(X_S)
  
  ggplot(data = NULL, aes(x = X_S[,1], y = X_S[,2])) + geom_point(size = .1) + theme_bw()
  
  Sigma_X <- K2(X_S, X_S, sigma_gp^2, l_gp)
  max(Sigma_X[1,-1])
}

