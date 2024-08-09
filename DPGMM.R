set.seed(123)

#------------------------------------------------ GENERATE SIMULATED DATASETS -------------------------------
# Function to generate synthetic data
generate_data <- function(n, means, sds) {
  data <- NULL
  for (i in seq_along(means)) {
    data <- c(data, rnorm(n, mean = means[i], sd = sds[i]))
  }
  return(data)
}

# Generate data from three different Gaussian distributions
n_per_cluster <- 100
means <- c(-5, 0, 5)
sds <- c(1, 1, 1)
data <- generate_data(n_per_cluster, means, sds)

# Plot the data
par(mfrow = c(1, 2)) 
hist(data, breaks = 30, col = 'brown', main = 'Generated Data') ; plot(data)


#------------------------------------------------ WRITE USEFUL FUNCTIONS -------------------------------
# Function to sample from the stick-breaking process
stick_breaking <- function(alpha, k) {
  betas <- rbeta(k, 1, alpha)
  pis <- numeric(k)
  pis[1] <- betas[1]
  for (i in 2:k) {
    pis[i] <- betas[i] * prod(1 - betas[1:(i-1)])
  }
  return(pis)
}

# Function to sample component parameters
sample_component_parameters <- function(base_mu, base_sigma, k) {
  mus <- rnorm(k, mean = base_mu, sd = base_sigma)
  sigmas <- rep(base_sigma, k)
  return(list(mus = mus, sigmas = sigmas))
}

# Function to assign data points to components
assign_data_to_components <- function(data, pis, mus, sigmas) {
  n <- length(data)
  k <- length(pis)
  z <- numeric(n)
  for (i in 1:n) {
    probs <- sapply(1:k, function(j) {
      pis[j] * dnorm(data[i], mean = mus[j], sd = sigmas[j])
    })
    z[i] <- sample(1:k, 1, prob = probs)
  }
  return(z)
}

#------------------------------------------------ IMPLEMENT DP-GMM -------------------------------


data = scale(feature_vector)

# Parameters for DP-GMM
alpha <- 1 # Concentration parameter
k <- 5  # Maximum number of components (for stick-breaking)

# Base distribution parameters for Gaussian components
base_mu <- 3
base_sigma <- 1

# Sample mixing weights using stick-breaking process
pis <- stick_breaking(alpha, k)

# Sample component parameters (means and standard deviations)
component_params <- sample_component_parameters(base_mu, base_sigma, k)
mus <- component_params$mus
sigmas <- component_params$sigmas

# Assign data points to components
z <- assign_data_to_components(data, pis, mus, sigmas)

par(mfrow = c(1, 1)) 
# Plot the data with cluster assignments
plot(data, col = z, pch = 19, main = 'DP-GMM Cluster Assignments')

# Print the compoSnent parameters
print(data.frame(Component = 1:k, Weight = pis, Mean = mus, SD = sigmas))


#------------------------------------------------ POST IMPLEMENTATION ANALYSIS -------------------------------

# Assuming z contains the predicted cluster labels and true_labels contains the ground truth labels
ari <- adjustedRandIndex(z, true_labels)
print(paste("Adjusted Rand Index:", ari))

silhouette_scores <- silhouette(z, dist(features))
average_silhouette_score <- mean(silhouette_scores[, 3])
print(paste("Average Silhouette Score:", average_silhouette_score))

db_index <- index.DB(features, z, centrotypes = "centroids")
print(paste("Davies-Bouldin Index:", db_index$DB))


