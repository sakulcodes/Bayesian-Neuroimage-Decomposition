# Introduce a significant amount of error to true_A to create estimated_A
set.seed(123)  # For reproducibility
estimated_A <- true_A + matrix(rnorm(length(true_A), mean = 0, sd = 1.2), nrow = nrow(true_A), ncol = ncol(true_A))

# Introduce a small amount of error to true_C to create estimated_C
estimated_C <- true_C + matrix(rnorm(length(true_C), mean = 0, sd = 0.4), nrow = nrow(true_C), ncol = ncol(true_C))

# Visualize the matrices
par(mfrow = c(2, 2))

# Plot true A matrix
image(true_A, main = "True A Matrix", xlab = "Components", ylab = "Features")

# Plot estimated A matrix
image(estimated_A, main = "Estimated A Matrix", xlab = "Components", ylab = "Features")

# Plot true C matrix
image(true_C, main = "True C Matrix", xlab = "Components", ylab = "Samples")

# Plot estimated C matrix
image(estimated_C, main = "Estimated C Matrix", xlab = "Components", ylab = "Samples")
