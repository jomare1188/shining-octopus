library(Matrix)

# Create the ortholog matrix (3 sorghum genes × 4 sugarcane genes)
ortholog_matrix <- Matrix(c(
    1, 1, 0, 0,  # Sb1 -> Sc1A, Sc1B
    0, 0, 1, 0,  # Sb2 -> Sc2
    0, 0, 0, 1   # Sb3 -> Sc3
), nrow=3, sparse=TRUE)

# Create sugarcane adjacency matrix (4×4)
sugarcane_matrix <- Matrix(c(
    0, 0, 1, 0,  # Sc1A -> Sc2
    0, 0, 1, 1,  # Sc1B -> Sc2, Sc3
    1, 1, 0, 0,  # Sc2 -> Sc1A, Sc1B
    0, 1, 0, 0   # Sc3 -> Sc1B
), nrow=4, sparse=TRUE)

# Create sorghum adjacency matrix (3×3)
sorghum_matrix <- Matrix(c(
    0, 1, 0,  # Sb1 -> Sb2
    1, 0, 1,  # Sb2 -> Sb1, Sb3
    0, 1, 0   # Sb3 -> Sb2
), nrow=3, sparse=TRUE)

# Add row and column names for better visualization
rownames(ortholog_matrix) <- c("Sb1", "Sb2", "Sb3")
colnames(ortholog_matrix) <- c("Sc1A", "Sc1B", "Sc2", "Sc3")

rownames(sugarcane_matrix) <- c("Sc1A", "Sc1B", "Sc2", "Sc3")
colnames(sugarcane_matrix) <- c("Sc1A", "Sc1B", "Sc2", "Sc3")

rownames(sorghum_matrix) <- c("Sb1", "Sb2", "Sb3")
colnames(sorghum_matrix) <- c("Sb1", "Sb2", "Sb3")

# Print the matrices to verify
cat("Ortholog Matrix:\n")
print(ortholog_matrix)

cat("\nSugarcane Network Matrix:\n")
print(sugarcane_matrix)

cat("\nSorghum Network Matrix:\n")
print(sorghum_matrix)

# Perform the matrix multiplication step by step
step1 <- ortholog_matrix %*% sugarcane_matrix
cat("\nStep 1 (ortholog_matrix %*% sugarcane_matrix):\n")
print(step1)

# Final multiplication with transposed ortholog matrix
mapped_edges <- (step1 %*% t(ortholog_matrix)) > 0
cat("\nMapped Edges (after converting to logical):\n")
print(mapped_edges)

# Compare with original sorghum network
cat("\nOriginal Sorghum Network:\n")
print(sorghum_matrix > 0)

# Show which edges are conserved
cat("\nConserved Edges (edges present in both networks):\n")
conserved <- (mapped_edges & (sorghum_matrix > 0))
print(conserved)
