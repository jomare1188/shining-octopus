suppressPackageStartupMessages({
    library(readr)
    library(igraph)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(cogeqc)
    library(optparse)
})

# Matrix construction and the ortholog projection live here; sourced relative
# to the working directory, so run this script from the repository root.
source("edge_conservation.r")

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

    # Add diagnostic prints. The gap between the first two numbers is the
    # orthology coverage of the network: edges whose endpoints have no ortholog
    # can never be counted as conserved.
    cat("  genes in the query network:      ",
        length(unique(c(sorghum_edges$gene1, sorghum_edges$gene2))), "\n")
    cat("  query genes with an ortholog:    ",
        length(unique(ortholog_pairs$sorghum_gene)), "\n")
    cat("  union used for the matrix:       ",
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

    # Create symmetric sparse adjacency matrices for both networks
    sorghum_matrix <- build_adjacency(
        sorghum_edges, sorghum_gene_to_idx, length(sorghum_genes))

    sugarcane_matrix <- build_adjacency(
        sugarcane_edges, sugarcane_gene_to_idx, length(sugarcane_genes))

    # Create ortholog mapping matrix
    ortholog_matrix <- sparseMatrix(
        i = sorghum_gene_to_idx[ortholog_pairs$sorghum_gene],
        j = sugarcane_gene_to_idx[ortholog_pairs$sugarcane_gene],
        x = 1,
        dims = c(length(sorghum_genes), length(sugarcane_genes))
    )
    # Find conserved edges through matrix multiplication
    # This maps sugarcane edges into sorghum space and checks for existence
    mapped_edges <- map_edges_through_orthologs(ortholog_matrix, sugarcane_matrix)

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

# Fraction of this species' edges whose orthologs are connected in the other.
#
# NOTE: despite the name this is a containment (overlap) coefficient, not a
# Jaccard index - the denominator is the edge count of the query network, not
# the size of the union. It is therefore asymmetric, which is why
# process_networks_bidirectional() reports it in both directions.
calculate_jaccard <- function(conserved_edges) {
    intersection_size <- sum(conserved_edges$conserved)
    total_edges <- nrow(conserved_edges)
    jaccard_index <- intersection_size / total_edges
    return(jaccard_index)
}

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

###########################################
# Inputs
###########################################

read_network <- function(path) {
    read_tsv(path, col_names = c("gene1", "gene2", "weight"),
             col_types = cols(gene1 = col_character(),
                              gene2 = col_character(),
                              weight = col_double()))
}

# Ortholog pairs between two species of the orthogroups table.
#
# NOTE: the sorghum_gene / sugarcane_gene column names are positional slots
# inherited from the original sorghum-vs-sugarcane script, not a claim about
# which species is which - sorghum_gene holds species_a and sugarcane_gene
# holds species_b whatever they are. Renaming them throughout would be clearer
# and is worth doing separately.
build_ortholog_pairs <- function(orthologues, species_a, species_b) {
    genes_a <- orthologues %>%
        filter(Species == species_a) %>%
        select(Orthogroup, sorghum_gene = Gene)

    genes_b <- orthologues %>%
        filter(Species == species_b) %>%
        select(Orthogroup, sugarcane_gene = Gene)

    inner_join(genes_a, genes_b, by = "Orthogroup")
}

###########################################
# Command line interface
###########################################

parse_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
    option_list <- list(
        make_option("--orthogroups", type = "character", default = NULL,
                    help = "OrthoFinder Orthogroups.tsv [required]"),
        make_option("--network-a", type = "character", default = NULL,
                    dest = "network_a",
                    help = paste("Edge list of the first species, tab separated,",
                                 "columns gene1 gene2 weight, no header [required]")),
        make_option("--network-b", type = "character", default = NULL,
                    dest = "network_b",
                    help = "Edge list of the second species [required]"),
        make_option("--species-a", type = "character", default = NULL,
                    dest = "species_a",
                    help = paste("Name of the first species exactly as it appears in the",
                                 "orthogroups file; see --list-species [required]")),
        make_option("--species-b", type = "character", default = NULL,
                    dest = "species_b",
                    help = "Name of the second species [required]"),
        make_option("--outdir", type = "character", default = ".",
                    help = "Directory for the classified edge tables [default %default]"),
        make_option("--gene-suffix", type = "character", default = "\\.p[0-9]+$",
                    dest = "gene_suffix",
                    help = paste("Regex stripped from gene IDs in the orthogroups file, to",
                                 "reconcile protein IDs with the gene IDs used in the",
                                 "networks. Pass '' to disable. [default %default]")),
        make_option("--list-species", action = "store_true", default = FALSE,
                    dest = "list_species",
                    help = "Print the species found in the orthogroups file and exit")
    )

    parse_args(OptionParser(
        option_list = option_list,
        usage = paste("usage: %prog --orthogroups FILE --network-a FILE --network-b FILE",
                      "--species-a NAME --species-b NAME [options]")),
        args = args)
}

require_options <- function(opt, names) {
    missing <- names[vapply(names, function(n) is.null(opt[[n]]), logical(1))]
    if (length(missing) > 0) {
        stop("missing required argument(s): ",
             paste0("--", gsub("_", "-", missing), collapse = ", "),
             "\nrun with --help for usage", call. = FALSE)
    }
}

# Species names come from the orthogroups column headers and are often long
# filenames, so they cannot be dropped into a path unescaped.
safe_filename <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Report how much of a network the orthogroups file actually covers, and refuse
# to run if the two use different identifier conventions.
#
# This is the failure mode worth being loud about: OrthoFinder runs on proteins
# and the networks are built on genes, so the IDs routinely differ by an isoform
# suffix. Nothing downstream errors when they do not match - every edge simply
# fails to find an ortholog and the conservation score comes out a clean, silent
# zero, which reads as a biological result rather than a configuration mistake.
check_id_overlap <- function(edges, ortholog_genes, species, gene_suffix) {
    network_genes <- unique(c(edges$gene1, edges$gene2))
    shared <- intersect(network_genes, unique(ortholog_genes))

    cat(sprintf("  %s: %d of %d network genes have an ortholog (%.1f%%)\n",
                species, length(shared), length(network_genes),
                100 * length(shared) / length(network_genes)))

    if (length(shared) == 0) {
        stop("no gene IDs are shared between the ", species, " network and the ",
             "orthogroups file.\n",
             "  network IDs look like:    ",
             paste(utils::head(network_genes, 3), collapse = ", "), "\n",
             "  orthogroup IDs look like: ",
             paste(utils::head(unique(ortholog_genes), 3), collapse = ", "), "\n",
             "  --gene-suffix is currently '", gene_suffix, "'; it is stripped from ",
             "the orthogroup IDs to reconcile the two.", call. = FALSE)
    }

    invisible(length(shared))
}

###########################################
# Entry point
###########################################

main <- function(opt = parse_arguments()) {
    require_options(opt, "orthogroups")

    orthologues <- read_orthogroups(opt$orthogroups)

    # read_orthogroups() returns Species as a factor. Left as one, cat() prints
    # its integer codes rather than the species names, and any comparison
    # against a name not in the levels is a warning-free NA.
    orthologues$Species <- as.character(orthologues$Species)

    if (nzchar(opt$gene_suffix)) {
        orthologues$Gene <- sub(opt$gene_suffix, "", orthologues$Gene)
    }

    available_species <- sort(unique(orthologues$Species))

    if (opt$list_species) {
        cat(available_species, sep = "\n")
        return(invisible(available_species))
    }

    require_options(opt, c("network_a", "network_b", "species_a", "species_b"))

    unknown <- setdiff(c(opt$species_a, opt$species_b), available_species)
    if (length(unknown) > 0) {
        stop("species not found in ", opt$orthogroups, ": ",
             paste(unknown, collapse = ", "),
             "\navailable: ", paste(available_species, collapse = ", "), call. = FALSE)
    }

    ortholog_pairs <- build_ortholog_pairs(orthologues, opt$species_a, opt$species_b)
    if (nrow(ortholog_pairs) == 0) {
        stop("no orthogroup contains both ", opt$species_a, " and ", opt$species_b,
             call. = FALSE)
    }
    cat("Ortholog pairs:", nrow(ortholog_pairs), "\n")

    edges_a <- read_network(opt$network_a)
    edges_b <- read_network(opt$network_b)

    cat("Orthology coverage\n")
    check_id_overlap(edges_a, ortholog_pairs$sorghum_gene, opt$species_a, opt$gene_suffix)
    check_id_overlap(edges_b, ortholog_pairs$sugarcane_gene, opt$species_b, opt$gene_suffix)
    cat("\n")

    results <- process_networks_bidirectional(edges_a, edges_b, ortholog_pairs)

    cat("\nConserved fraction",
        "\n  ", opt$species_a, " -> ", opt$species_b, ": ",
        format(results$sorghum_to_sugarcane, digits = 4),
        "\n  ", opt$species_b, " -> ", opt$species_a, ": ",
        format(results$sugarcane_to_sorghum, digits = 4), "\n", sep = "")

    dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
    out_a <- file.path(opt$outdir,
                       paste0("classification_edges_", safe_filename(opt$species_a), ".csv"))
    out_b <- file.path(opt$outdir,
                       paste0("classification_edges_", safe_filename(opt$species_b), ".csv"))

    write.csv(results$edges_sorghum, out_a, quote = FALSE, row.names = FALSE)
    write.csv(results$edges_sugarcane, out_b, quote = FALSE, row.names = FALSE)
    cat("\nWrote ", out_a, "\n      ", out_b, "\n", sep = "")

    invisible(results)
}

# Runs only under `Rscript shared_edges.r`; sourcing the file does not trigger it.
if (sys.nframe() == 0L) {
    main()
}
