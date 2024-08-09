# -------------------- BAYESIAN INFERENCE (GIBBS SAMPLING) --------------------- 
# (i) ----------------------- Initialization -------------------- 
V = 10; N = 4; R = 3; n_iter = 1000

# Hyperparameters
alpha_a <- 2 ; gamma_a <- 2
alpha_c <- 2 ; gamma_c <- 2
alpha_lambda <- 2 ; gamma_lambda <- 2
alpha_tau <- 2 ; gamma_tau <- 2

# Initialize parameters (in log scale)
log_A <- matrix(rnorm(V * R), nrow = V)  # log(A)
log_C <- matrix(rnorm(N * R), nrow = N)  # log(C)
log_lambda <- runif(R, log(0.5), log(1.5))  # log(lambda)
tau <- rgamma(1, alpha_tau, gamma_tau)
beta_a <- rgamma(1, alpha_a, gamma_a)
beta_c <- rgamma(1, alpha_c, gamma_c)
beta_lambda <- rgamma(1, alpha_lambda, gamma_lambda)

# Storage for samples
samples_A <- array(0, dim = c(n_iter, V, R))
samples_C <- array(0, dim = c(n_iter, N, R))
samples_lambda <- matrix(0, nrow = n_iter, ncol = R)
samples_tau <- numeric(n_iter)
samples_beta_a <- numeric(n_iter)
samples_beta_c <- numeric(n_iter)
samples_beta_lambda <- numeric(n_iter)

# (ii) ----------------------- Gibbs Sampler --------------------
for (iter in 1:n_iter) {
  print(iter)
  
  # Sample log(A)
  for (r in 1:R) {
    for (v in 1:V) {
      precision_log_A <- tau * exp(2 * log_lambda[r]) * sum(exp(2 * log_C[, r])) + beta_a
      sigma_log_A_cond <- sqrt(1 / precision_log_A)
      mu_log_A <- sigma_log_A_cond^2 * tau * exp(log_lambda[r]) * sum(X[v, , ] * exp(log_A[, r]) * exp(log_C[, r]))
      log_A[v, r] <- rnorm(1, mu_log_A, sigma_log_A_cond)
    }
  }
  
  # Sample log(C)
  for (r in 1:R) {
    for (n in 1:N) {
      precision_log_C <- tau * exp(2 * log_lambda[r]) * sum(exp(2 * log_A[, r])) + beta_c
      sigma_log_C_cond <- sqrt(1 / precision_log_C)
      mu_log_C <- sigma_log_C_cond^2 * tau * exp(log_lambda[r]) * sum(X[, , n] * exp(log_A[, r]) * exp(log_A[, r]))
      log_C[n, r] <- rnorm(1, mu_log_C, sigma_log_C_cond)
    }
  }
  
  # Sample log(lambda)
  for (r in 1:R) {
    precision_log_lambda <- tau * sum(exp(4 * log_A[, r]) * exp(2 * log_C[, r])) + beta_lambda
    sigma_log_lambda_cond <- sqrt(1 / precision_log_lambda)
    mu_log_lambda <- sigma_log_lambda_cond^2 * tau * sum(X * exp(log_A[, r]) * exp(log_A[, r]) * exp(log_C[, r]))
    log_lambda[r] <- rnorm(1, mu_log_lambda, sigma_log_lambda_cond)
  }
  
  # Sample tau
  alpha_tau_post <- alpha_tau + V * V * N / 2
  fitted <- array(0, dim = c(V, V, N))
  for (r in 1:R) {
    fitted <- fitted + exp(log_lambda[r]) * exp(log_A[, r] %o% log_A[, r] %o% log_C[, r])
  }
  beta_tau_post <- gamma_tau + 0.5 * sum((X - fitted)^2)
  tau <- rgamma(1, alpha_tau_post, beta_tau_post)
  
  # Sample beta_a, beta_c, beta_lambda
  beta_a <- rgamma(1, alpha_a + V * R / 2, gamma_a + sum(exp(2 * log_A)) / 2)
  beta_c <- rgamma(1, alpha_c + N * R / 2, gamma_c + sum(exp(2 * log_C)) / 2)
  beta_lambda <- rgamma(1, alpha_lambda + R / 2, gamma_lambda + sum(exp(2 * log_lambda)) / 2)
  
  # Store samples (on the original scale)
  samples_A[iter,,] <- exp(log_A)
  samples_C[iter,,] <- exp(log_C)
  samples_lambda[iter,] <- exp(log_lambda)
  samples_tau[iter] <- tau
  samples_beta_a[iter] <- beta_a
  samples_beta_c[iter] <- beta_c
  samples_beta_lambda[iter] <- beta_lambda
}
