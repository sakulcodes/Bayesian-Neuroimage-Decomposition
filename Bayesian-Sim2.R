library(MASS)

# Helper function to sample from Inverse-Gamma
rinvgamma <- function(n, alpha, beta) {
  return(1 / rgamma(n, alpha, 1 / beta))
}

outer_func <- function(a, b, c) {
  result1 <- outer(a, b)  # Outer product of vectors a and b, resulting in a V x V matrix
  V <- length(a)
  N <- length(c)
  final_result <- array(0, dim = c(V, V, N))  # Initialize the result array
  
  for (i in 1:N) {
    final_result[,,i] <- result1 * c[i]  # Scale matrix by each element in c
  }
  
  return(final_result)
}


generate_sim_tensor <- function(V, N, R, block_size) {
  # Initialize matrices A and C
  A <- matrix(0, nrow = V, ncol = R)
  C <- matrix(0, nrow = N, ncol = R)
  
  # Calculate number of blocks along the diagonal
  num_blocks_A <- V %/% block_size
  num_blocks_C <- N %/% block_size
  
  # Create block structure in A and C
  for (r in 1:R) {
    for (block in 1:num_blocks_A) {
      row_start <- (block - 1) * block_size + 1
      row_end <- min(block * block_size, V)
      
      # Set some blocks in A to non-zero values, leaving others as zero
      if (block %% 2 == 1) {  # Non-zero block on odd indices
        A[row_start:row_end, r] <- rnorm(row_end - row_start + 1, mean = 5, sd = 0.1)
      }
    }
    
    for (block in 1:num_blocks_C) {
      col_start <- (block - 1) * block_size + 1
      col_end <- min(block * block_size, N)
      
      # Set some blocks in C to non-zero values, leaving others as zero
      if (block %% 2 == 1) {  # Non-zero block on odd indices
        C[col_start:col_end, r] <- rnorm(col_end - col_start + 1, mean = 5, sd = 0.1)
      }
    }
  }
  
  # Generate tensor data using the block-structured factor matrices
  lambda <- runif(R, 0.5, 1.5)  # Random weights for the components
  X <- array(0, dim = c(V, V, N))
  error <- array(rnorm(V * V * N, mean = 0, sd = 0.1), dim = c(V, V, N))
  for (r in 1:R) {
    X <- X + lambda[r] * outer_func(A[, r], A[, r], C[, r])
  }
  
  data <- X + error
  return(list(dataTensor = data, A = A, C = C))
}


# Generate Simulated Data:
sim_data <- generate_sim_tensor(V = 20, N = 100, R = 5, block_size = 1)
X <- sim_data$dataTensor

# Gibbs Sampler for the tensor decomposition model
gibbs_sampler <- function(X, V, N, R, n_iter) {
  # Hyperparameters
  alpha_A <- 2
  beta_A <- 2
  alpha_C <- 2
  beta_C <- 2
  alpha <- 2
  beta <- 2
  
  # Initialize parameters
  A <- matrix(rnorm(V * R), nrow = V)
  C <- matrix(rnorm(N * R), nrow = N)
  lambda <- runif(R, 0.5, 1.5)
  sigma2 <- 0.1
  sigma_A2 <- 1
  sigma_C2 <- 1
  
  # Storage for samples
  samples_A <- array(0, dim = c(n_iter, V, R))
  samples_C <- array(0, dim = c(n_iter, N, R))
  samples_lambda <- matrix(0, nrow = n_iter, ncol = R)
  samples_sigma2 <- numeric(n_iter)
  samples_sigma_A2 <- numeric(n_iter)
  samples_sigma_C2 <- numeric(n_iter)
  
  # Gibbs Sampling
  for (iter in 1:n_iter) {
    
    # Sample A
    for (r in 1:R) {
      for (i in 1:V) {
        # Compute conditional variance
        precision_A = 1 / sigma_A^2 + lambda[r]^2 / sigma2 * sum(A[, r]^2 * C[, r]^2)
        sigma_A_cond = sqrt(1 / precision_A)
        
        # Compute conditional mean
        mu_A = sigma_A_cond^2 * (lambda[r] / sigma2 * sum(X[i,,] * A[, r] * C[, r]))
        
        # Sample from the conditional distribution
        A[i, r] <- rnorm(1, mu_A, sigma_A_cond)
      }
    }
    
    # Sample C
    for (r in 1:R) {
      for (k in 1:N) {
        # Conditional variance
        precision_C <- 1 / sigma_C2 + lambda[r]^2 / sigma2 * sum(A[, r]^2 * A[, r]^2)
        sigma_C_cond <- sqrt(1 / precision_C)
        
        # Conditional mean
        mu_C <- sigma_C_cond^2 * (lambda[r] / sigma2 * sum(X[,,k] * A[, r] * A[, r]))
        
        # Sample from the conditional distribution
        C[k, r] <- rnorm(1, mu_C, sigma_C_cond)
      }
    }
    
    # Sample lambda
    for (r in 1:R) {
      # Conditional variance
      precision_lambda <- 1 / sigma_A2 + 1 / sigma2 * sum(A[, r]^2 * A[, r]^2 * C[, r]^2)
      sigma_lambda_cond <- sqrt(1 / precision_lambda)
      
      # Conditional mean
      mu_lambda <- sigma_lambda_cond^2 * (1 / sigma2 * sum(X * A[, r] * A[, r] * C[, r]))
      
      # Sample from the conditional distribution
      lambda[r] <- rnorm(1, mu_lambda, sigma_lambda_cond)
    }
    
    # Sample sigma2
    alpha_post <- alpha + V * V * N / 2
    fitted <- array(0, dim = c(V, V, N))
    for (r in 1:R) {
      fitted <- fitted + lambda[r] * outer_func(A[, r], A[, r], C[, r])
    }
    beta_post <- beta + sum((X - fitted)^2) / 2
    sigma2 <- rinvgamma(1, alpha_post, beta_post)
    
    # Sample sigma_A2
    alpha_A_post <- alpha_A + V * R / 2
    beta_A_post <- beta_A + sum(A^2) / 2
    sigma_A2 <- rinvgamma(1, alpha_A_post, beta_A_post)
    
    # Sample sigma_C2
    alpha_C_post <- alpha_C + N * R / 2
    beta_C_post <- beta_C + sum(C^2) / 2
    sigma_C2 <- rinvgamma(1, alpha_C_post, beta_C_post)
    
    # Store samples
    samples_A[iter,,] <- A
    samples_C[iter,,] <- C
    samples_lambda[iter,] <- lambda
    samples_sigma2[iter] <- sigma2
    samples_sigma_A2[iter] <- sigma_A2
    samples_sigma_C2[iter] <- sigma_C2
  }
  
  return(list(A = samples_A, C = samples_C, lambda = samples_lambda, sigma2 = samples_sigma2,
              sigma_A2 = samples_sigma_A2, sigma_C2 = samples_sigma_C2))
}



# Run Gibbs Sampler:
set.seed(123)
results <- gibbs_sampler(X, V = 20, N = 100, R = 5, n_iter = 10000)



#------------------------------- Convergence Analysis :  --------------------

library(coda)  # For diagnostics and plotting


# Compute posterior means
posterior_mean_A <- apply(results$A, c(2, 3), mean)
posterior_mean_C <- apply(results$C, c(2, 3), mean)

# Extract original matrices
original_A <- sim_data$A
original_C <- sim_data$C

# Compute RMSE for A
rmse_A <- sqrt(mean((posterior_mean_A - original_A)^2))
# Compute the Frobenius norm of the original A
norm_original_A <- sqrt(sum(original_A^2))
# Compute relative RMSE for A
rRMSE_A <- rmse_A / norm_original_A
print(paste("RMSE for A:", rmse_A))
print(paste("Relative RMSE for A:", rRMSE_A))

# Compute RMSE for C
rmse_C <- sqrt(mean((posterior_mean_C - original_C)^2))
# Compute the Frobenius norm of the original C
norm_original_C <- sqrt(sum(original_C^2))
# Compute relative RMSE for C
rRMSE_C <- rmse_C / norm_original_C
print(paste("RMSE for C:", rmse_C))
print(paste("Relative RMSE for C:", rRMSE_C))

#--------------For Lambda:
par(mfrow = c(3, 1))
plot(results$lambda[,1], type = "l") ; plot(results$lambda[,2], type = "l") ; plot(results$lambda[,3], type = "l")
plot(results$lambda[,4], type = "l")
plot(results$lambda[,5], type = "l")

#------------For Precision (Sigma^2):
par(mfrow = c(1, 1))


#For Hyperparameters:

plot(results$sigma_A2, type = "l")
plot(results$sigma_C2, type = "l")


#------------------ FORBENIUS NORM --------------
# Function to calculate the Frobenius norm of the difference between two matrices
frobenius_norm <- function(mat1, mat2) {
  return(sqrt(sum((mat1 - mat2)^2)))
}

# Calculate the Frobenius norm of differences for A and C matrices over iterations
calc_diffs <- function(samples_A, samples_C, n_iter) {
  diff_A <- numeric(n_iter - 1)
  diff_C <- numeric(n_iter - 1)
  
  for (iter in 2:n_iter) {
    diff_A[iter - 1] <- frobenius_norm(samples_A[iter,,], samples_A[iter - 1,,])
    diff_C[iter - 1] <- frobenius_norm(samples_C[iter,,], samples_C[iter - 1,,])
  }
  
  return(list(diff_A = diff_A, diff_C = diff_C))
}

# Get the differences
diffs <- calc_diffs(results$A, results$C, n_iter =20000)

# Plot the differences
plot(diffs$diff_A, type = "l", col = "blue", ylab = "Frobenius norm difference", xlab = "Iteration", main = "Convergence for Factor Matrix A")
plot(diffs$diff_C, type = "l", ylab = "Frobenius norm difference", xlab = "Iteration", main = "Convergence for Factor Matrix C")



# ------------ IMAGING -----------

# Extract true matrices
true_A <- sim_data$A ; true_C <- sim_data$C

# Calculate posterior means of A and C
posterior_A <- apply(results$A, c(2, 3), mean) ; posterior_C <- apply(results$C, c(2, 3), mean)

# Visualize the true and estimated matrices using image
par(mfrow = c(2, 2))  # Set up a 2x2 plotting area

# Plot true A
image(true_A, main = "True A", xlab = "Components", ylab = "Features")
# Plot estimated A
image(posterior_A, main = "Estimated A", xlab = "Components", ylab = "Features")

# Plot true C
image(true_C, main = "True C", xlab = "Components", ylab = "Samples")
# Plot estimated C
image(posterior_C, main = "Estimated C", xlab = "Components", ylab = "Samples")
par(mfrow = c(1, 1)) 


#Estimation - Imagery

