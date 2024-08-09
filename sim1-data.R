library(dirichletprocess) ; library(stats) ; library(cluster); library(clusterSim); library(multiway); library(mclust)

# Example list of VxV brain connectivity matricies for N subjects
V <- 200 # Number of nodes
N <- 300  # Number of subjects
set.seed(123)

# Function to create a symmetric matrix with cluster-specific patterns
create_cluster_symmetric_matrix <- function(V, cluster_mean) {
  mat <- matrix(rnorm(V * V, mean = cluster_mean), nrow = V)
  sym_mat <- (mat + t(mat)) / 2  # Making the matrix symmetric
  diag(sym_mat) <- abs(diag(sym_mat))  # Ensuring diagonal elements are non-negative
  return(sym_mat)
}

# Establishing Ground Truth
K = 5 # no. of clusters
cluster_means = seq(1:K)  
true_labels = sample(1:K,N, replace = TRUE)

brain_connectivity_matrices <- lapply(1:N, function(i) {
  cluster <- true_labels[i]
  create_cluster_symmetric_matrix(V, cluster_means[cluster])
})

# Add Gaussian noise to each matrix
noise_sd <- 0.01  # Standard deviation of the noise
brain_connectivity_matrices_noisy <- lapply(brain_connectivity_matrices, function(mat) {
  mat + matrix(rnorm(V * V, mean = 0, sd = noise_sd), nrow = V)
})


brain_connectivity_tensor <- array(unlist(brain_connectivity_matrices_noisy), dim = c(V, V, N))

#Parafac Decomposition

parafac_result <- parafac(brain_connectivity_tensor, nfac = 13) #Alternating Least Squares (ALS) algorithm. 
features = parafac_result$C

#--------------------- ELBOW METHOD FOR K MEANS --------------------------

# Function to calculate total within-cluster sum of squares
wcss <- function(features, max_k) {
  wcss_values <- numeric(max_k)
  for (k in 1:max_k) {
    kmeans_result <- kmeans(features, centers = k, nstart = 25)
    wcss_values[k] <- sum(kmeans_result$tot.withinss)
  }
  return(wcss_values)
}

# Set maximum number of clusters to check
max_k <- 15

# Calculate WCSS for each k from 1 to max_k
wcss_values <- wcss(features, max_k)

# Plot the WCSS values to visualize the Elbow Method
plot(1:max_k, wcss_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters K",
     ylab = "Total Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Optimal K")


#-------------------- K MEANS ALGORITHM WITH K = 5 --------------------------
set.seed(123)
feature_vec  = features[,1]
kmeans_result_1 <- kmeans(feature_vec, centers = K)
kmeans_result_2 = kmeans(features, centers = K)

# Predicted labels from K-means
predicted_labels_1 <- kmeans_result_1$cluster; predicted_labels_2 <- kmeans_result_2$cluster

comparison_table_1 <- table(Predicted = predicted_labels_1, True = true_labels)
print(comparison_table) #THE RESULT IS OKAYISH

comparison_table_2 <- table(Predicted = predicted_labels_2, True = true_labels)
print(comparison_table) #THE RESULT IS OKAYISH


#--------------------- DP CLUSTERING -----------------------


# dp <- DirichletProcessMvnormal(scale(features))

feature_vector = features[,1]
dp <- DirichletProcessGaussian(scale(feature_vector))
# Fit the model
dp <- Fit(dp, 2000)

# Print the results
print(dp)

# Retrieve the cluster assignments
cluster_assignments <- dp$clusterLabels
print(cluster_assignments)

# Compare the predicted labels with the true labels
comparison_table <- table(Predicted = cluster_assignments, True = true_labels)
print(comparison_table)

   #--------------------- BEST RANK? -----------------------

# Function to calculate explained variance
explained_variance <- function(tensor, model) {
  # Total sum of squares of the original tensor
  total_ss <- sum(tensor^2)
  
  # Reconstructed tensor from the model
  reconstructed_tensor <- array(0, dim = dim(tensor))
  for (r in 1:length(model$A[1, ])) {
    reconstructed_tensor <- reconstructed_tensor + outer(model$A[, r], model$B[, r]) %o% model$C[, r]
  }
  
  # Sum of squares of the reconstructed tensor
  model_ss <- sum(reconstructed_tensor^2)
  
  # Explained variance
  explained_variance <- model_ss / total_ss
  return(explained_variance)
}

# Trying different ranks and calculating explained variance
ranks <- 1:15
explained_var <- sapply(ranks, function(r) {
  model <- parafac(brain_connectivity_tensor, nfac = r)
  explained_variance(brain_connectivity_tensor, model)
})

# Plot explained variance against rank
plot(ranks, explained_var, type = "b", xlab = "Rank", ylab = "Explained Variance", main = "Explained Variance vs. Rank")
abline(v = which.max(diff(diff(explained_var))), col = "red", lty = 2)


#--------------------- SHILHOUETTE AND DAVIES-BOULDIN SCORES -----------------------'

#------------------------------------(I) K MEANS:
# Calculate silhouette scores
silhouette_scores <- silhouette(kmeans_result$cluster, dist(features))
cat("Average silhouette width:", mean(silhouette_scores[, 'sil_width']), "\n")

# Calculate Davies-Bouldin Index
dbi <- index.DB(features, kmeans_result$cluster, centrotypes = "centroids")
cat("Davies-Bouldin Index:", dbi$DB, "\n")

ari <- adjustedRandIndex(kmeans_result$cluster, true_labels)
print(paste("Adjusted Rand Index:", ari))

#------------------------------------(II) DIRICHLET PROCESS:

# Calculate silhouette scores for DP clustering
silhouette_scores_dp <- silhouette(cluster_assignments, dist(feature_vector))
cat("Average silhouette width for DP clustering:", mean(silhouette_scores_dp[, 'sil_width']), "\n")

dbi_dp <- index.DB(features, cluster_assignments, centrotypes = "centroids")
cat("Davies-Bouldin Index for DP clustering:", dbi_dp$DB, "\n")

ari <- adjustedRandIndex(cluster_assignments, true_labels)
print(paste("Adjusted Rand Index:", ari))


#--------------------- PLOTS -----------------------'
centers =  kmeans_result$centers
plot(features[,1], features[,2], xlab = "1st Rank", ylab = "2nd Rank", main = "" ,col = "black")
# Overlay the cluster centers
points(centers[,1], centers[,2], pch = 4, col = "blue", cex = 2, lwd = 2)

# Add a legend
legend("topleft", legend = c("Data Points", "Cluster Centers"), col = c("black", "blue"), pch = c(19, 4), cex = 0.8)






