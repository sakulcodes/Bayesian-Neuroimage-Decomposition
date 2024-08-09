

#-------------------- SIMULATED BRAIN CONNECTOME DATA --------------------------
library(dirichletprocess) ; library(stats) ; library(cluster); library(clusterSim)

# Example list of VxV brain connectivity matricies for N subjects
V <- 200 # Number of nodes
N <- 100  # Number of subjects
set.seed(123)

# Function to create a symmetric matrix
create_symmetric_matrix <- function(V) {
  mat <- matrix(rnorm(V * V), nrow = V)
  sym_mat <- (mat + t(mat)) / 2  # Making the matrix symmetric
  diag(sym_mat) <- abs(diag(sym_mat))  # Ensuring diagonal elements are non-negative
  return(sym_mat)
}

# Generate the list of symmetric brain connectivity matrices
brain_connectivity_matrices <- replicate(N, create_symmetric_matrix(V), simplify = FALSE)

# Function to perform SVD and extract singular values
extract_singular_values <- function(matrix) {
  svd_result <- svd(matrix)
  return(svd_result$d)
}

# Extract features (singular values) for each subject
features <- t(sapply(brain_connectivity_matrices, extract_singular_values))

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

#--------------------- K MEANS CLUSTERING  --------------------------


# Perform K-means clustering
num_clusters <- 4
kmeans_result <- kmeans(features, centers = num_clusters)

# Print clustering result
print(kmeans_result$cluster)

# Calculate silhouette scores
silhouette_scores <- silhouette(kmeans_result$cluster, dist(features))

# Print average silhouette width
cat("Average silhouette width:", mean(silhouette_scores[, 'sil_width']), "\n")

# Plot silhouette scores
#plot(silhouette_scores)

# Calculate Davies-Bouldin Index
dbi <- index.DB(features, kmeans_result$cluster, centrotypes = "centroids")
cat("Davies-Bouldin Index:", dbi$DB, "\n")



#--------------------- DIRICHLET PROCESS CLUSTERING -----------------------


dp <- DirichletProcessMvnormal(as.matrix(scale(features)))#, alphaPrior = c(40, 0.001)

# Fit the model
dp <- Fit(dp, 1000)

# Print the results
print(dp)

# Retrieve the cluster assignments
cluster_assignments <- dp$clusterLabels
print(cluster_assignments)


#--------------------- SILHOUETTE SCORE CALCULATION FOR DP CLUSTERING --------------------------

# Calculate silhouette scores for DP clustering
silhouette_scores_dp <- silhouette(cluster_assignments, dist(features))

# Print average silhouette width for DP clustering
cat("Average silhouette width for DP clustering:", mean(silhouette_scores_dp[, 'sil_width']), "\n")

# Plot silhouette scores for DP clustering
plot(silhouette_scores_dp, main = "Silhouette Plot for DP Clustering")

#--------------------- DAVIES-BOULDIN INDEX CALCULATION FOR DP CLUSTERING --------------------------

# Calculate Davies-Bouldin Index for DP clustering
dbi_dp <- index.DB(features, cluster_assignments, centrotypes = "centroids")
cat("Davies-Bouldin Index for DP clustering:", dbi_dp$DB, "\n")

#----------------------------IMAGES--------------------------
x = brain_connectivity_matrices[[1]]

# Visualize the matrix using the image function
image(x, main = "Brain Connectivity Matrix",
      xlab = "Brain Nodes", ylab = "Brain Nodes")


#------------------ DEBUGGING -------------
V = 5; N =3
brain_connectivity_matrices <- replicate(N, create_symmetric_matrix(V), simplify = FALSE)



mat = matrix(data = 1:25, nrow = 5, byrow = TRUE)
