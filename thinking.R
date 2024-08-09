# Load the required library
library(multiway)

# Create a simple symmetric matrix (representing a brain connectome)
set.seed(123)
V <- 5 ; N <- 4
tensor_data <- array(0, dim = c(V, V, N))
for (i in 1:N) {
  connectome_matrix <- matrix(sample(1:10, V*V, replace=TRUE), nrow=V)
  connectome_matrix[lower.tri(connectome_matrix)] <- t(connectome_matrix)[lower.tri(connectome_matrix)]
  tensor_data[,,i] <- connectome_matrix 
}


parafac_decomp <- parafac(tensor_data, nfac = 3)

# Print the decomposition result
parafac_decomp$A; parafac_decomp$B
