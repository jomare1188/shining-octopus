library(readr)
library(igraph)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(cogeqc)

# Read and preprocess data as in original code
orthologues <- read_orthogroups("/home/dmpachon/jorge/comparative_cane/results/orthofinder/Results_Aug30/Orthogroups/Orthogroups.tsv")
orthologues$Gene <- sub("\\.p[0-9]+$", "", orthologues$Gene)
orthologues$Species <- recode(orthologues$Species,
    maize_GCF_902167145.1_Zm.B73.REFERENCE.NAM.5.0_protein = "maize",
    rice_GCF_034140825.1_ASM3414082v1_protein = "rice",
    sorghum_protein_of_longest_cds_per_orthogroup = "sorghum",
    sugarcane_proteins_of_longest_cds_per_OG = "sugarcane"
)

edges_sorghum <- read_tsv("../results/networks/abs_triplets_p60_sorghum_fancy.tsv",
                             col_names = c("gene1", "gene2", "weight"))
edges_sugarcane <- read_tsv("../results/networks/abs_triplets_p60_cane_fancy.tsv",
                               col_names = c("gene1", "gene2", "weight"))

# Create ortholog pairs as in original code
sorghum_orthologues <- orthologues %>%
    filter(Species == "sorghum") %>%
    select(Orthogroup, sorghum_gene = Gene)

sugarcane_orthologues <- orthologues %>%
    filter(Species == "sugarcane") %>%
    select(Orthogroup, sugarcane_gene = Gene)

ortholog_pairs <- sorghum_orthologues %>%
    inner_join(sugarcane_orthologues, by = "Orthogroup")

###########################################
# Matrix-Based Approach Implementation
###########################################

process_networks_matrix <- function(sorghum_edges, sugarcane_edges, ortholog_pairs) {
    # Get unique gene lists from both edges and orthologs
    sorghum_genes <- unique(c(
        sorghum_edges$gene1, 
        sorghum_edges$gene2,
        ortholog_pairs$sorghum_gene  # Include genes from ortholog pairs
    ))
    sugarcane_genes <- unique(c(
        sugarcane_edges$gene1, 
        sugarcane_edges$gene2,
        ortholog_pairs$sugarcane_gene  # Include genes from ortholog pairs
    ))
    
    # Add diagnostic prints
    cat("Number of unique sorghum genes in edges:", 
        length(unique(c(sorghum_edges$gene1, sorghum_edges$gene2))), "\n")
    cat("Number of unique sorghum genes in orthologs:", 
        length(unique(ortholog_pairs$sorghum_gene)), "\n")
    cat("Total unique sorghum genes after combining:", 
        length(sorghum_genes), "\n")
    
    # Create gene to index mappings
    sorghum_gene_to_idx <- setNames(seq_along(sorghum_genes), sorghum_genes)
    sugarcane_gene_to_idx <- setNames(seq_along(sugarcane_genes), sugarcane_genes)
    
    # Verify no NAs in the mapping
    if(any(is.na(sorghum_gene_to_idx[ortholog_pairs$sorghum_gene]))) {
        missing_genes <- ortholog_pairs$sorghum_gene[
            is.na(sorghum_gene_to_idx[ortholog_pairs$sorghum_gene])
        ]
        stop("Some sorghum genes in ortholog pairs are missing from the mapping: ",
             paste(head(missing_genes), collapse=", "))
    }
    
    if(any(is.na(sugarcane_gene_to_idx[ortholog_pairs$sugarcane_gene]))) {
        missing_genes <- ortholog_pairs$sugarcane_gene[
            is.na(sugarcane_gene_to_idx[ortholog_pairs$sugarcane_gene])
        ]
        stop("Some sugarcane genes in ortholog pairs are missing from the mapping: ",
             paste(head(missing_genes), collapse=", "))
    }
    
    # Create sparse adjacency matrices for both networks
    sorghum_matrix <- sparseMatrix(
        i = sorghum_gene_to_idx[sorghum_edges$gene1],
        j = sorghum_gene_to_idx[sorghum_edges$gene2],
        x = sorghum_edges$weight,
        dims = c(length(sorghum_genes), length(sorghum_genes))
    )
    
    sugarcane_matrix <- sparseMatrix(
        i = sugarcane_gene_to_idx[sugarcane_edges$gene1],
        j = sugarcane_gene_to_idx[sugarcane_edges$gene2],
        x = sugarcane_edges$weight,
        dims = c(length(sugarcane_genes), length(sugarcane_genes))
    )
    
    # Create ortholog mapping matrix
    ortholog_matrix <- sparseMatrix(
        i = sorghum_gene_to_idx[ortholog_pairs$sorghum_gene],
        j = sugarcane_gene_to_idx[ortholog_pairs$sugarcane_gene],
        x = 1,
        dims = c(length(sorghum_genes), length(sugarcane_genes))
    )        
    # Find conserved edges through matrix multiplication
    # This maps sorghum edges to sugarcane space and checks for existence
    mapped_edges <- (ortholog_matrix %*% sugarcane_matrix %*% t(ortholog_matrix)) > 0
    
    # Convert results back to edge format
    conserved_edges <- sorghum_edges %>%
        mutate(
            edge_idx = row_number(),
            conserved = mapped_edges[cbind(
                sorghum_gene_to_idx[gene1],
                sorghum_gene_to_idx[gene2]
            )]
        )
    
    return(conserved_edges)
}

# Calculate shared edges using either implementation
calculate_jaccard <- function(conserved_edges) {
    intersection_size <- sum(conserved_edges$conserved)
    total_edges <- nrow(conserved_edges)
    jaccard_index <- intersection_size / total_edges
    return(jaccard_index)
}

# Example usage:
# Matrix approach
#matrix_results <- process_networks_matrix(adjacency_sorghum, adjacency_sugarcane, ortholog_pairs)
#matrix_jaccard <- calculate_jaccard(matrix_results)

# Functiion to calculate sorghum to cane and cane to sorghum
process_networks_bidirectional <- function(sorghum_edges, sugarcane_edges, ortholog_pairs) {
    # Original direction (sorghum to sugarcane)
    sorghum_to_sugarcane <- process_networks_matrix(
        sorghum_edges, sugarcane_edges, ortholog_pairs)
    
    # Reverse direction (sugarcane to sorghum)
    # Need to swap the columns in ortholog_pairs
    reversed_pairs <- ortholog_pairs %>%
        select(sorghum_gene = sugarcane_gene, 
               sugarcane_gene = sorghum_gene)
    
    sugarcane_to_sorghum <- process_networks_matrix(
        sugarcane_edges, sorghum_edges, reversed_pairs)
    
    # Compare results
    sorghum_to_sugarcane_jaccard <- calculate_jaccard(sorghum_to_sugarcane)
    sugarcane_to_sorghum_jaccard <- calculate_jaccard(sugarcane_to_sorghum)
    
    return(list(
        sorghum_to_sugarcane = sorghum_to_sugarcane_jaccard,
        sugarcane_to_sorghum = sugarcane_to_sorghum_jaccard,
	edges_sorghum = sorghum_to_sugarcane,
	edges_sugarcane = sugarcane_to_sorghum
    ))
}

# example run
results <- process_networks_bidirectional(edges_sorghum, edges_sugarcane, ortholog_pairs)

# write in files the classified edges 
edges_sorghum <- results$edges_sorghum
edges_sugarcane <- results$edges_sugarcane

write.csv(edges_sorghum, "classification_edges_sorghum.csv", quote = F, row.names = F)
write.csv(edges_sugarcane, "classification_edges_sugarcane.csv", quote = F, row.names = F)
