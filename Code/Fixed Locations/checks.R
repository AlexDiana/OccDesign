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
  plotgradient(gradient1_psi)
  
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
n <- 50

# total number of sampling occasions
n_occ <- 100

# generate sites locations
X <- cbind(sqrt(sqrt(seq(0, 1, length.out = n))), 1)
ggplot(data = NULL, aes(
  x = X[,1], y = X[,2], size = 1)) + geom_point()

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
