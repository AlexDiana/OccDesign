
library(ggplot2); library(here); library(tidyverse)
library(foreach); library(doParallel)
library(OccDesign)
# source(here("R","utility.R"))
# source(here("R","fun.R"))
numCores <- detectCores()
registerDoParallel(numCores)

N_x <- 500 # number of simulation replicates
N_s <- 100 # utlity estimation replicates

niter <- 100

setwd(here("Output"))

# scenario 1 (no effects) ----------

# true params
{
  beta_psi_true <- seq(0, 0, length.out = 4)
  p0_true <- 0.2
  beta_p_true <- 0
  coeffs_true <- list("psi" = beta_psi_true,
                      "p" = c(logit(p0_true),beta_p_true))
}

# create gradients
{
  set.seed(1)
  gradient1_psi <- createLinearGradient()
  plotgradient(gradient1_psi)
  
  gradPsi_minmax <- findGradientMax(gradient1_psi)
  gradPsi_min <- gradPsi_minmax[1]
  gradPsi_max <- gradPsi_minmax[2]
  
  X_psibreaks <- c(gradPsi_min - .05, 
                   seq(gradPsi_min, gradPsi_max, length.out = 5)[-c(1,5)],
                   gradPsi_max + .05)  #c(-1,.25, .5, .75,1)
  
  set.seed(3)
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
n <- 20

# total number of sampling occasions
n_occ <- 200

# generate sites locations
{
  # X <- expand.grid(seq(0, 1, length.out = sqrt(n)),
  # seq(0, 1, length.out = sqrt(n)))
  X1 <- cbind(seq(0, 1, length.out = n), .5)
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  # X <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  
}

list_design_scenario1 <- findOptimalDesign(niter, 
                                           n, n_occ,
                                           N_x, N_s,
                                           coeffs_true,
                                           X1, X_psibreaks,
                                           gradient_list)
utility_vals <- list_design_scenario1$utility_vals
pi_vals <- list_design_scenario1$pi_vals

diag1 <- plotDiagnostics(pi_vals)
ggsave("diag_1.jpeg", diag1)
obj1 <- plotObjective(utility_vals)
ggsave("obj_1.jpeg", obj1)
design1 <- plotDesign(pi_vals, n_occ, X1, gradient1_psi)
ggsave("design_1.jpeg", design1)

# scenario 2 (positive effect on occ) ----------

# true params
{
  beta_psi_true <- seq(-1, 1, length.out = 4)
  p0_true <- .2
  beta_p_true <- c(0)
  coeffs_true <- list("psi" = beta_psi_true,
                      "p" = c(logit(p0_true),beta_p_true))
}

# create gradients
{
  set.seed(1)
  gradient1_psi <- createLinearGradient()
  plotgradient(gradient1_psi)
  
  gradPsi_minmax <- findGradientMax(gradient1_psi)
  gradPsi_min <- gradPsi_minmax[1]
  gradPsi_max <- gradPsi_minmax[2]
  
  X_psibreaks <- c(gradPsi_min - .05, 
                   seq(gradPsi_min, gradPsi_max, length.out = 5)[-c(1,5)],
                   gradPsi_max + .05)  #c(-1,.25, .5, .75,1)
  
  set.seed(3)
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
n <- 20

# total number of sampling occasions
n_occ <- 200

# generate sites locations
{
  # X <- expand.grid(seq(0, 1, length.out = sqrt(n)),
  # seq(0, 1, length.out = sqrt(n)))
  X <- cbind(seq(0, 1, length.out = n), .5)
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  # X <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  
}

list_design_scenario2 <- findOptimalDesign(niter, 
                                           n, n_occ, 
                                           N_x, N_s,
                                           coeffs_true,
                                           X, X_psibreaks,
                                           gradient_list)

utility_vals2 <- list_design_scenario2$utility_vals
pi_vals2 <- list_design_scenario2$pi_vals

diag2 <- plotDiagnostics(pi_vals2)
ggsave("diag_2.jpeg", diag2)
obj2 <- plotObjective(utility_vals2)
ggsave("obj_2.jpeg", obj2)
design2 <- plotDesign(pi_vals2, n_occ, X, gradient1_psi)
ggsave("design_2.jpeg", design2)

# scenario 3 (positive effect on det) ----------

# true params
{
  beta_psi_true <- seq(0, 0, length.out = 4)
  p0_true <- .2
  beta_p_true <- c(3)
  coeffs_true <- list("psi" = beta_psi_true,
                      "p" = c(logit(p0_true),beta_p_true))
}

# create gradients
{
  set.seed(1)
  gradient1_psi <- createLinearGradient()
  plotgradient(gradient1_psi)
  
  gradPsi_minmax <- findGradientMax(gradient1_psi)
  gradPsi_min <- gradPsi_minmax[1]
  gradPsi_max <- gradPsi_minmax[2]
  
  X_psibreaks <- c(gradPsi_min - .05, 
                   seq(gradPsi_min, gradPsi_max, length.out = 5)[-c(1,5)],
                   gradPsi_max + .05)  #c(-1,.25, .5, .75,1)
  
  set.seed(3)
  gradient1_p <- createLinearGradient()
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
n <- 20

# total number of sampling occasions
n_occ <- 200

# generate sites locations
{
  # X <- expand.grid(seq(0, 1, length.out = sqrt(n)),
  # seq(0, 1, length.out = sqrt(n)))
  X <- cbind(seq(0, 1, length.out = n), .5)
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  # X <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  
}

list_design_scenario3 <- findOptimalDesign(niter,
                                           n, n_occ, 
                                           N_x, N_s,
                                           coeffs_true,
                                           X, X_psibreaks,
                                           gradient_list)

utility_vals3 <- list_design_scenario3$utility_vals
pi_vals3 <- list_design_scenario3$pi_vals

diag3 <- plotDiagnostics(pi_vals3)
ggsave("diag_3.jpeg", diag3)
obj3 <- plotObjective(utility_vals3)
ggsave("obj_3.jpeg", obj3)
design3 <- plotDesign(pi_vals3, n_occ, X, gradient1_psi)
ggsave("design_3.jpeg", design3)

# scenario 4 (no effect, non uniform design) ----------

# true params
{
  beta_psi_true <- seq(0, 0, length.out = 4)
  p0_true <- .2
  beta_p_true <- c(0)
  coeffs_true <- list("psi" = beta_psi_true,
                      "p" = c(logit(p0_true),beta_p_true))
}

# create gradients
{
  set.seed(1)
  gradient1_psi <- createLinearGradient()
  plotgradient(gradient1_psi)
  
  gradPsi_minmax <- findGradientMax(gradient1_psi)
  gradPsi_min <- gradPsi_minmax[1]
  gradPsi_max <- gradPsi_minmax[2]
  
  X_psibreaks <- c(gradPsi_min - .05, 
                   seq(gradPsi_min, gradPsi_max, length.out = 5)[-c(1,5)],
                   gradPsi_max + .05)  #c(-1,.25, .5, .75,1)
  
  set.seed(3)
  gradient1_p <- createLinearGradient()
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
n <- 20

# total number of sampling occasions
n_occ <- 200

# generate sites locations
{
  # X <- expand.grid(seq(0, 1, length.out = sqrt(n)),
  # seq(0, 1, length.out = sqrt(n)))
  # X <- cbind(seq(0, 1, length.out = n), .5)
  X <- cbind(sqrt(seq(0, 1, length.out = n)), .5)
  # X <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  
}

list_design_scenario4 <- findOptimalDesign(niter, 
                                           n, n_occ, 
                                           N_x, N_s,
                                           coeffs_true,
                                           X, X_psibreaks,
                                           gradient_list)

utility_vals4 <- list_design_scenario4$utility_vals
pi_vals4 <- list_design_scenario4$pi_vals

diag4 <- plotDiagnostics(pi_vals4)
ggsave("diag_4.jpeg", diag4)
obj4 <- plotObjective(utility_vals4)
ggsave("obj_4.jpeg", obj4)
design4 <- plotDesign(pi_vals4, n_occ, X, gradient1_psi)
ggsave("design_4.jpeg", design4)

# scenario 5 (negative effect on prob, non uniform design) ----------

# true params
{
  beta_psi_true <- seq(0, 0, length.out = 4)
  p0_true <- .2
  beta_p_true <- c(-3)
  coeffs_true <- list("psi" = beta_psi_true,
                      "p" = c(logit(p0_true),beta_p_true))
}

# create gradients
{
  set.seed(1)
  gradient1_psi <- createLinearGradient()
  plotgradient(gradient1_psi)
  
  gradPsi_minmax <- findGradientMax(gradient1_psi)
  gradPsi_min <- gradPsi_minmax[1]
  gradPsi_max <- gradPsi_minmax[2]
  
  X_psibreaks <- c(gradPsi_min - .05, 
                   seq(gradPsi_min, gradPsi_max, length.out = 5)[-c(1,5)],
                   gradPsi_max + .05)  #c(-1,.25, .5, .75,1)
  
  set.seed(3)
  gradient1_p <- createLinearGradient()
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
n <- 20

# total number of sampling occasions
n_occ <- 200

# generate sites locations
{
  # X <- expand.grid(seq(0, 1, length.out = sqrt(n)),
  # seq(0, 1, length.out = sqrt(n)))
  # X <- cbind(seq(0, 1, length.out = n), .5)
  X <- cbind(sqrt(seq(0, 1, length.out = n)), .5)
  # X <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  # X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
  
}

list_design_scenario5 <- findOptimalDesign(niter, 
                                           n, n_occ, 
                                           N_x, N_s,
                                           coeffs_true,
                                           X, X_psibreaks,
                                           gradient_list)

utility_vals5 <- list_design_scenario5$utility_vals
pi_vals5 <- list_design_scenario5$pi_vals

diag5 <- plotDiagnostics(pi_vals5)
ggsave("diag_5.jpeg", diag5)
obj5 <- plotObjective(utility_vals5)
ggsave("obj_5.jpeg", obj5)
design5 <- plotDesign(pi_vals5, n_occ, X, gradient1_psi)
ggsave("design_5.jpeg", design5)
