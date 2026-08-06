# Regression tests for the matrix-based edge-conservation method.
#
# Depends on Matrix only, so it runs without the tidyverse/OrthoFinder stack:
#     Rscript test_edge_conservation.r

source("edge_conservation.r")

failures <- 0

check <- function(label, ok) {
    cat(if (isTRUE(ok)) "PASS  " else "FAIL  ", label, "\n", sep = "")
    if (!isTRUE(ok)) failures <<- failures + 1
}

index_of <- function(genes) setNames(seq_along(genes), genes)

# Fraction of a network's edges whose orthologs are connected in the other one.
conservation_rate <- function(edges, gene_to_idx, mapped) {
    flags <- mapped[cbind(gene_to_idx[edges$gene1], gene_to_idx[edges$gene2])]
    sum(flags) / nrow(edges)
}

###########################################
# build_adjacency
###########################################

genes <- c("g1", "g2", "g3")
idx <- index_of(genes)
edges <- data.frame(gene1 = c("g1", "g2"), gene2 = c("g2", "g3"), weight = c(1, 1))

# build_adjacency leaves the matrix unnamed, as the pipeline does - dimnames on
# a 100k-gene matrix are pure overhead there. Name it here for readable tests.
named_adjacency <- function(edges, gene_to_idx, genes) {
    adjacency <- build_adjacency(edges, gene_to_idx, length(genes))
    dimnames(adjacency) <- list(genes, genes)
    adjacency
}

adjacency <- named_adjacency(edges, idx, genes)

check("adjacency of an undirected edge list is symmetric", isSymmetric(adjacency))
check("both orientations of each edge are filled",
      adjacency["g1", "g2"] == 1 && adjacency["g2", "g1"] == 1)

# A repeated edge must not have its weight summed into something larger.
duplicated_edges <- data.frame(gene1 = c("g1", "g1"), gene2 = c("g2", "g2"),
                               weight = c(1, 1))
check("a duplicated edge does not inflate its weight",
      named_adjacency(duplicated_edges, idx, genes)["g1", "g2"] == 1)

###########################################
# Orientation independence
#
# This is the regression guard: building only the (gene1, gene2) triangle made
# the conservation score depend on the order genes happened to appear in the
# input file, so a network compared with itself could score 0.
###########################################

identity_orthology <- Diagonal(length(genes))
mapped_self <- map_edges_through_orthologs(identity_orthology, adjacency)

check("a network compared against itself is fully conserved",
      conservation_rate(edges, idx, mapped_self) == 1)

flipped <- data.frame(gene1 = edges$gene2, gene2 = edges$gene1, weight = edges$weight)
check("conservation is unchanged when the edge list is written the other way round",
      conservation_rate(flipped, idx, mapped_self) == 1)

check("projecting a flipped edge list gives the same result",
      all(as.matrix(map_edges_through_orthologs(
              identity_orthology, build_adjacency(flipped, idx, length(genes)))) ==
          as.matrix(mapped_self)))

###########################################
# Worked example from the README
#
# Sb1 has two sugarcane orthologs (Sc1A, Sc1B); Sb2 -> Sc2; Sb3 -> Sc3.
###########################################

sorghum_genes <- c("Sb1", "Sb2", "Sb3")
sugarcane_genes <- c("Sc1A", "Sc1B", "Sc2", "Sc3")
sb_idx <- index_of(sorghum_genes)
sc_idx <- index_of(sugarcane_genes)

ortholog_matrix <- sparseMatrix(
    i = sb_idx[c("Sb1", "Sb1", "Sb2", "Sb3")],
    j = sc_idx[c("Sc1A", "Sc1B", "Sc2", "Sc3")],
    x = 1,
    dims = c(length(sorghum_genes), length(sugarcane_genes)),
    dimnames = list(sorghum_genes, sugarcane_genes)
)

sorghum_edges <- data.frame(gene1 = c("Sb1", "Sb2"), gene2 = c("Sb2", "Sb3"),
                            weight = c(1, 1))
sugarcane_edges <- data.frame(gene1 = c("Sc1A", "Sc1B", "Sc1B"),
                              gene2 = c("Sc2", "Sc2", "Sc3"),
                              weight = c(1, 1, 1))

sorghum_adjacency <- named_adjacency(sorghum_edges, sb_idx, sorghum_genes)
sugarcane_adjacency <- named_adjacency(sugarcane_edges, sc_idx, sugarcane_genes)

counts <- ortholog_matrix %*% sugarcane_adjacency %*% t(ortholog_matrix)
mapped <- map_edges_through_orthologs(ortholog_matrix, sugarcane_adjacency)

# Sb1-Sb2 is reachable through two ortholog pairs (Sc1A-Sc2 and Sc1B-Sc2),
# Sb1-Sb3 through one (Sc1B-Sc3), and Sb2-Sb3 through none.
expected_counts <- matrix(c(0, 2, 1,
                            2, 0, 0,
                            1, 0, 0), nrow = 3, byrow = TRUE,
                          dimnames = list(sorghum_genes, sorghum_genes))
check("ortholog-pair counts match the worked example",
      all(as.matrix(counts) == expected_counts))

conserved <- mapped & (sorghum_adjacency > 0)
check("Sb1-Sb2 is conserved", conserved["Sb1", "Sb2"])
check("Sb2-Sb3 is not conserved (Sc2 and Sc3 are unconnected)",
      !conserved["Sb2", "Sb3"])
check("Sb1-Sb3 is ortholog-implied, not conserved",
      mapped["Sb1", "Sb3"] && !conserved["Sb1", "Sb3"])
check("conservation rate of the worked example is 1/2",
      conservation_rate(sorghum_edges, sb_idx, mapped) == 0.5)

###########################################

cat("\n", if (failures == 0) "all checks passed\n" else
    paste0(failures, " check(s) failed\n"), sep = "")
quit(status = if (failures == 0) 0 else 1)
