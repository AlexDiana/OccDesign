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

simulateCovariatesLocation <- function(n){
  
  X <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  
  # create gradients
  set.seed(1)
  # gradient1_psi <- createLinearGradient()
  gradient1_psi <- createGradient()
  plotgradient(gradient1_psi)
  
  gradPsi_minmax <- findGradientMax(gradient1_psi)
  gradPsi_min <- gradPsi_minmax[1]
  gradPsi_max <- gradPsi_minmax[2]
  
  X_psirange <- gradPsi_minmax
  
  # X_psibreaks <- c(gradPsi_min - .05, 
  #                  seq(gradPsi_min, gradPsi_max, length.out = 5)[-c(1,5)],
  #                  gradPsi_max + .05)  #c(-1,.25, .5, .75,1)
  
  set.seed(3)
  # gradient1_p <- createLinearGradient()
  # gradient1_p <- createUniformGradient()
  gradient1_p <- createGradient()
  plotgradient(gradient1_p)
  
  gradient_psi <- list(
    gradient1_psi
  )
  
  gradient_p <- list(
    gradient1_p
  )
  
  X_psi <- sapply(gradient_psi, function(gradient){
    apply(X, 1, gradient)
  })
  
  X_p <- sapply(gradient_p, function(gradient){
    apply(X, 1, gradient)
  })
  
  list("X_psi" = X_psi,
       "X_p" = X_p,
       "X_psirange" = X_psirange)
   
}

# SETTINGS ------

# fixed design parameter

designSettings <- list(
  n = 25, # number of sites
  n_occ = 50 # number of sampling occasions
)

# true parameters (theta_0)

beta_psi_true <- seq(0, 0, length.out = 4)
ppsi0_true <- .5
beta_p_true <- c(0)

trueParamSettings <-  list(
  "psi" = beta_psi_true,
  "p" = c(logit(ppsi0_true),beta_p_true)
  )


# simulate covariates locations
{
  list_X <- simulateCovariatesLocation(designSettings$n)
  X_psi <- list_X$X_psi
  X_p <- list_X$X_p
  X_psirange <- list_X$X_psirange
}

# covariates values
covariatesValues <- list(
  X_psi = X_psi, # covariates for occupancy
  X_psirange = X_psirange, # 
  X_p = X_p # covariates for detection
  )

algoParams <- list(
  eps = .01,
  N_x = 100, # number of simulation replicates
  N_s = 100 # utlity estimation replicates
)

## ALGORITHM -------

optimalConfiguration <- 
  findOptimalDesign(
    designSettings,
    trueParamSettings,
    covariatesValues,
    algoParams
  )

# CHECKS -----

# loss function for uniform setting

n <- designSettings$n
n_occ <- designSettings$n_occ

M <- rep(n_occ / n, n)

computeUtilityBaseFinal(n, M, 
                   X_psi, X_p, maxM,
                   trueParamSettings,
                   X_psibreaks,
                   numSims = 10000)

# loss function for optimal setting

n <- designSettings$n
n_occ <- designSettings$n_occ

M <- round(optimalConfiguration)

computeUtilityBaseFinal(n, M, 
                   X_psi, X_p, maxM,
                   trueParamSettings,
                   X_psibreaks,
                   numSims = 10000)
