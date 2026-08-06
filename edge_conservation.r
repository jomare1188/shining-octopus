# Core of the matrix-based edge-conservation method.
#
# Depends on Matrix only, so it can be sourced and tested without the
# tidyverse/OrthoFinder stack that the full pipeline in shared_edges.r needs.

library(Matrix)

# Build a symmetric sparse adjacency matrix from an undirected edge list.
#
# Edge lists store each undirected edge once (gene1, gene2), so filling only
# those (i, j) positions yields a triangular, non-symmetric matrix. Every step
# downstream then depends on the arbitrary gene order of the input file: a
# network compared against itself scores 1.0 or 0.0 depending only on whether
# the two files happen to list each pair the same way round. Filling both
# orientations restores the symmetry an undirected network requires.
#
# use.last.ij = TRUE keeps duplicated (i, j) entries at their last value rather
# than summing them, so a repeated edge - or a self-loop, which is emitted in
# both orientations - cannot inflate a weight.
build_adjacency <- function(edges, gene_to_idx, n_genes) {
    i <- gene_to_idx[edges$gene1]
    j <- gene_to_idx[edges$gene2]

    adjacency <- sparseMatrix(
        i = c(i, j),
        j = c(j, i),
        x = c(edges$weight, edges$weight),
        dims = c(n_genes, n_genes),
        use.last.ij = TRUE
    )

    stopifnot(isSymmetric(adjacency))
    adjacency
}

# Project a network through an orthology relation.
#
# ortholog_matrix is the n_A x n_B indicator matrix of orthologous pairs and
# adjacency_b the adjacency matrix of species B. Entry [i, k] of the triple
# product counts the edges of B running between any ortholog of gene i and any
# ortholog of gene k, so a positive entry means at least one pair of their
# orthologs is connected. Thresholding at > 0 turns that count into the
# many-to-many OR the method is after.
#
# The result is symmetric whenever adjacency_b is, which is what makes the
# per-edge lookup in process_networks_matrix() independent of edge orientation.
map_edges_through_orthologs <- function(ortholog_matrix, adjacency_b) {
    (ortholog_matrix %*% adjacency_b %*% t(ortholog_matrix)) > 0
}
