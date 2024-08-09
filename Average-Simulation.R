

library(dirichletprocess) ; library(stats) ; library(cluster); library(clusterSim); library(multiway); library(mclust)
# Number of replications
num_replications <- 10

# Initialize lists to store results
kmeans_results <- list()
dp_results <- list()

# Function to calculate clustering accuracy
calculate_accuracy <- function(true_labels, predicted_labels) {
  comparison_table <- table(Predicted = predicted_labels, True = true_labels)
  correct_predictions <- sum(diag(comparison_table))
  total_predictions <- sum(comparison_table)
  accuracy <- correct_predictions / total_predictions
  return(accuracy)
}

for (rep in 1:num_replications) {
  # Simulate data
  true_labels <- sample(1:K, N, replace = TRUE)
  brain_connectivity_matrices <- lapply(1:N, function(i) {
    cluster <- true_labels[i]
    create_cluster_symmetric_matrix(V, cluster_means[cluster])
  })
  
  brain_connectivity_matrices_noisy <- lapply(brain_connectivity_matrices, function(mat) {
    mat + matrix(rnorm(V * V, mean = 0, sd = noise_sd), nrow = V)
  })
  
  brain_connectivity_tensor <- array(unlist(brain_connectivity_matrices_noisy), dim = c(V, V, N))
  
  # PARAFAC decomposition
  parafac_result <- parafac(brain_connectivity_tensor, nfac = 13)
  features <- parafac_result$C
  
  # K-means clustering
  set.seed(123)
  feature_vec <- features[, 1]
  kmeans_result <- kmeans(feature_vec, centers = K)
  predicted_labels_kmeans <- kmeans_result$cluster
  kmeans_accuracy <- calculate_accuracy(true_labels, predicted_labels_kmeans)
  kmeans_results[[rep]] <- kmeans_accuracy
  
  # DP clustering
  feature_vector <- scale(features[, 1])
  dp <- DirichletProcessGaussian(feature_vector)
  dp <- Fit(dp, 2000)
  cluster_assignments <- dp$clusterLabels
  dp_accuracy <- calculate_accuracy(true_labels, cluster_assignments)
  dp_results[[rep]] <- dp_accuracy
}

# Calculate average results
avg_kmeans_accuracy <- mean(unlist(kmeans_results))
avg_dp_accuracy <- mean(unlist(dp_results))

# Print the average accuracies
cat("Average K-means Accuracy:", avg_kmeans_accuracy, "\n")
cat("Average DP Clustering Accuracy:", avg_dp_accuracy, "\n")
