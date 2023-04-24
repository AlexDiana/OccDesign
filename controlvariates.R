cummean <- function(x){
  cumsum(x) / 1:length(x)
}

# grad_cummean <- apply(grad_est, 1, cummean)
# 
# data_plot <- grad_cummean %>% as.data.frame() %>% 
#   mutate(Iter = row_number()) %>% 
#   pivot_longer(cols = starts_with("V")) 
# 
# ggplot(data_plot[data_plot$name %in% c("V1","V2","V3"),], aes(x = Iter, y = value, color = name)) + 
#   geom_point()

#  -----

Stheta_x <- S_theta(H_x)

grad_est <- sapply(1:N_x, function(i){
  x_current <- x[,i]
  Stheta_x[i] * grad_logfX(x_current, theta, n_occ)
})

apply(grad_est, 1, var)

grad_est2 <- sapply(1:N_x, function(i){
  x_current <- x[,i]
  grad_logfX(x_current, theta, n_occ)
})

apply(grad_est2, 1, mean)

var_h <- apply(grad_est2, 1, var)
cov_h <- sapply(1:n, function(i){
  cov(grad_est[i,], grad_est2[i,])
})

cor_h <- sapply(1:n, function(i){
  cor(grad_est[i,], grad_est2[i,])
})

a_i <- cov_h / var_h

grad_est3 <- sapply(1:N_x, function(i){
  x_current <- x[,i]
  Stheta_x[i] * grad_logfX(x_current, theta, n_occ) - 
    grad_logfX(x_current, theta, n_occ) * a_i
})

# multivariate
{
  # VarM <- cov(t(grad_est2))
  # CovMT <- cov(t(grad_est), t(grad_est2))
  # C <- - CovMT %*% solve(VarM)
  # 
  # grad_est3 <- sapply(1:N_x, function(i){
  #   x_current <- x[,i]
  #   Stheta_x[i] * grad_logfX(x_current, theta, n_occ) +
  #     t(C) %*% grad_logfX(x_current, theta, n_occ) 
  # })
}

apply(grad_est, 1, var)
apply(grad_est3, 1, var)

apply(grad_est, 1, mean)
apply(grad_est3, 1, mean)

grad_mean <- apply(grad_est, 1, mean)

grad_mean3 <- apply(grad_est3, 1, mean)

grad_cummean <- apply(grad_est, 1, cummean)

data_plot <- grad_cummean %>% as.data.frame() %>% 
  mutate(Iter = row_number()) %>% 
  pivot_longer(cols = starts_with("V")) 

ggplot(data_plot[data_plot$name %in%
                   # c("V1","V2","V3","V4","V5","V6")
                 c("V46","V47","V48","V49","V50")
                 & data_plot$Iter > 50,], aes(x = Iter, y = value, color = name)) + 
  geom_point()

qplot(grad_est[2,])

qplot(1:(n-1), grad_mean)
qplot(1:(n-1), grad_mean3)

# HESSIAN ----------

H_logfX(theta, n_occ)

Hessian_list <- lapply(1:N_x, function(i){
  x_current <- x[,i]
  Stheta_x[i] * (H_logfX(theta, n_occ) + grad_logfX(x_current, theta, n_occ) %*%
                   t(grad_logfX(x_current, theta, n_occ)))
})

Hessian_list2 <- lapply(1:N_x, function(i){
  x_current <- x[,i]
  grad_logfX(x_current, theta, n_occ) %*%
                   t(grad_logfX(x_current, theta, n_occ))
})

Hessian_vals <- array(NA, dim = c(N_x, n, n))
Hessian_vals2 <- array(NA, dim = c(N_x, n, n))
for (i in 1:N_x) {
  Hessian_vals[i,,] <- Hessian_list[[i]]
  Hessian_vals2[i,,] <- Hessian_list2[[i]]
}

Hessian_vals_var <- apply(Hessian_vals, c(2,3), var)
Hessian_vals2_var <- apply(Hessian_vals2, c(2,3), var)

# var_h <- apply(grad_est2, 1, var)
cov_h <- matrix(NA, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    cov_h[i,j] <- cov(Hessian_vals[,i,j],
                      Hessian_vals2[,i,j])
  }
}

a_i <- cov_h / Hessian_vals2_var

Hessian_list3 <- lapply(1:N_x, function(i){
  x_current <- x[,i]
  Stheta_x[i] * (H_logfX(theta, n_occ) + grad_logfX(x_current, theta, n_occ) %*%
                   t(grad_logfX(x_current, theta, n_occ))) - 
    a_i * grad_logfX(x_current, theta, n_occ) %*%
    t(grad_logfX(x_current, theta, n_occ))
})

Hessian_vals3 <- array(NA, dim = c(N_x, n, n))
for (i in 1:N_x) {
  Hessian_vals3[i,,] <- Hessian_list3[[i]]
}

Hessian_val_mean <- apply(Hessian_vals3, c(2,3), mean)

apply(Hessian_vals, c(2,3), var)[1,1]
apply(Hessian_vals3, c(2,3), var)[1,1]

qplot(1:n, diag(Hessian_val_mean))

qplot(1:n, solve(Hessian_val_mean) %*% grad_mean3)

qplot(1:n, theta) + ylim(c(-1.5, .2))
qplot(1:n, theta + solve(Hessian_val_mean) %*% grad_mean3) + 
  ylim(c(-1.5, .2))

qplot(1:n, theta)
