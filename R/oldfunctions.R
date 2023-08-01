

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
