library(MASS) ; library(invgamma) ; library(base)


generate_sim_tensor = function(V,N,R){
  # Generate Random Factor Matricies (Since A = B, we just generate A)
  A <- abs(matrix(rnorm(V * R), nrow = V)) ; C <- abs(matrix(rnorm(N * R), nrow = N))
  
  # Generate Tensor Data:
  lambda <- runif(R, 0.5, 1.5)  # Random weights for the components
  X <- array(0, dim = c(V, V, N))
  error <- array(rnorm(V*V*N, mean = 0, sd = 0.1), dim = c(V,V,N))
  for (r in 1:R) {
    X = X + lambda[r] * (A[,r] %o% A[,r] %o% C[,r])
  }
  data = X + error
  return(list(dataTensor=data,A=A,C=C))
}

sim_data <- generate_sim_tensor(V = 10, N = 4, R = 3)
X = sim_data$dataTensor

# -------------------- BAYESIAN INFERENCE (GIBBS SAMPLING) --------------------- 
# (i) ----------------------- Initialization -------------------- 
V = 10; N =4; R = 3; n_iter = 1000
# Hyperparameters
alpha_a <- 2 ; gamma_a <- 2
alpha_c <- 2 ; gamma_c <- 2
alpha_lambda <- 2 ; gamma_lambda <- 2
alpha_tau <- 2 ;gamma_tau <- 2

# Initialize parameters
A <- abs(matrix(rnorm(V * R), nrow = V))
C <- abs(matrix(rnorm(N * R), nrow = N))
lambda <- runif(R, 0.5, 1.5)
tau <- rgamma(1, alpha_tau, gamma_tau) #Precision of error
beta_a <- rgamma(1, alpha_a, gamma_a) #Precision of A
beta_c <- rgamma(1, alpha_c, gamma_c) #Precision of C
beta_lambda <- rgamma(1, alpha_lambda, gamma_lambda) #Precision of lambda

# Storage for samples
samples_A <- array(0, dim = c(n_iter, V, R))
samples_C <- array(0, dim = c(n_iter, N, R))
samples_lambda <- matrix(0, nrow = n_iter, ncol = R)
samples_tau <- numeric(n_iter)
samples_beta_a <- numeric(n_iter)
samples_beta_c <- numeric(n_iter)
samples_beta_lambda <- numeric(n_iter)

# (ii) ----------------------- GIbbs Sampler (Updates) --------------------
for (iter in 1:n_iter) {
  print(iter)
  for (r in 1:R) {
    for (v in 1:V) {
      precision_A <- tau * lambda[r]^2 * sum(C[, r]^2) + beta_a
      sigma_A_cond <- sqrt(1 / precision_A)
      mu_A <- sigma_A_cond^2 * tau * lambda[r] * sum(X[v, , ] * A[, r] * C[, r])
      A[v, r] <- rnorm(1, mu_A, sigma_A_cond)
    }
  }
  
  # Sample C
  for (r in 1:R) {
    for (n in 1:N) {
      precision_C <- tau * lambda[r]^2 * sum(A[, r]^2) + beta_c
      sigma_C_cond <- sqrt(1 / precision_C)
      mu_C <- sigma_C_cond^2 * tau * lambda[r] * sum(X[, , n] * A[, r] * A[, r])
      C[n, r] <- rnorm(1, mu_C, sigma_C_cond)
    }
  }
  
  #Sample Lambda
  for (r in 1:R) {
    precision_lambda <- tau * sum(A[, r]^2 %o% C[, r]^2) + beta_lambda
    sigma_lambda_cond <- sqrt(1 / precision_lambda)
    mu_lambda <- sigma_lambda_cond^2 * tau * sum(outer(X, A[, r]) * C[, r])
    lambda[r] <- rnorm(1, mu_lambda, sigma_lambda_cond)
  }
  
  # Sample beta_a, beta_c, beta_lambda
  beta_a <- rgamma(1, alpha_a + V * R / 2, gamma_a + sum(A^2) / 2)
  beta_c <- rgamma(1, alpha_c + N * R / 2, gamma_c + sum(C^2) / 2)
  beta_lambda <- rgamma(1, alpha_lambda + R / 2, gamma_lambda + sum(lambda^2) / 2)
  
  #Sample Tau
  residuals <- X
  for (r in 1:R) {
    contribution <- lambda[r] * (A[, r] %o% A[, r] %o% C[, r])
    residuals <- residuals - contribution
  }
  residual_sum_squares <- sum(residuals^2)
  alpha_tau_prime <- alpha_tau + (V * N) / 2
  beta_tau_prime <- beta_tau + residual_sum_squares / 2
  tau <- rgamma(1, shape = alpha_tau_prime, rate = beta_tau_prime)
  
  
  a# Store samples
  samples_A[iter,,] <- A
  samples_C[iter,,] <- C
  samples_lambda[iter,] <- lambda
  samples_tau[iter] <- tau
  samples_beta_a[iter] <- beta_a
  samples_beta_c[iter] <- beta_c
  samples_beta_lambda[iter] <- beta_lambda
}



x = rgamma(100,shape = 10,rate = 10) ; hist(x)
