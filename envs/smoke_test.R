# Smoke test for the matrix2-expression environment.
#
# Loading a package is not the same as it working - bioconductor packages built
# against the wrong R version import fine and then fail inside compiled code.
# This runs each tool on tiny synthetic data instead.
#
#     conda activate matrix2-expression
#     Rscript envs/smoke_test.R

suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(topGO)
    library(tximport)
    library(Matrix)
})

set.seed(42)

n_genes <- 200
n_samples <- 6
counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 1/0.2),
                 nrow = n_genes,
                 dimnames = list(paste0("gene", seq_len(n_genes)),
                                 paste0("sample", seq_len(n_samples))))
condition <- factor(rep(c("control", "treated"), each = 3))

# Plant a real effect in the first 20 genes so the tools have something to find.
counts[1:20, condition == "treated"] <-
    counts[1:20, condition == "treated"] * 8

###########################################
# DESeq2
###########################################

dds <- DESeqDataSetFromMatrix(counts, DataFrame(condition), ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)
cat("DESeq2      ", sum(res$padj < 0.05, na.rm = TRUE), "DE genes of 20 planted\n")

# lfcShrink exercises apeglm, which is compiled and a common breakage point.
shrunk <- lfcShrink(dds, coef = "condition_treated_vs_control",
                    type = "apeglm", quiet = TRUE)
cat("apeglm      shrunk LFC range",
    sprintf("%.2f to %.2f", min(shrunk$log2FoldChange), max(shrunk$log2FoldChange)), "\n")

###########################################
# edgeR
###########################################

dge <- DGEList(counts = counts, group = condition)
dge <- calcNormFactors(dge)
design <- model.matrix(~ condition)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
cat("edgeR       ", sum(topTags(qlf, n = Inf)$table$FDR < 0.05), "DE genes of 20 planted\n")

###########################################
# topGO
###########################################

# Three fake GO terms; the first is assigned to the genes we made DE.
gene2GO <- setNames(
    lapply(seq_len(n_genes), function(i)
        if (i <= 20) c("GO:0008150", "GO:0003674") else "GO:0005575"),
    rownames(counts))

gene_scores <- setNames(as.numeric(res$padj), rownames(counts))
gene_scores[is.na(gene_scores)] <- 1

go_data <- new("topGOdata", ontology = "BP",
               allGenes = gene_scores,
               geneSelectionFun = function(p) p < 0.05,
               annot = annFUN.gene2GO, gene2GO = gene2GO,
               nodeSize = 1)
go_res <- runTest(go_data, algorithm = "classic", statistic = "fisher")
cat("topGO        ran Fisher test,", length(score(go_res)), "GO term(s) scored\n")

###########################################
# tximport
###########################################

# tximport needs quantification files on disk, so exercise it on a minimal
# fake salmon quant.sf rather than just checking that the function exists.
quant_dir <- file.path(tempdir(), "salmon_sample1")
dir.create(quant_dir, showWarnings = FALSE)
quant_file <- file.path(quant_dir, "quant.sf")
write.table(
    data.frame(Name = paste0("tx", 1:10), Length = 1000L, EffectiveLength = 800,
               TPM = runif(10, 0, 100), NumReads = rpois(10, 50)),
    quant_file, sep = "\t", quote = FALSE, row.names = FALSE)

tx2gene <- data.frame(TXNAME = paste0("tx", 1:10),
                      GENEID = rep(paste0("gene", 1:5), each = 2))
txi <- tximport(setNames(quant_file, "sample1"), type = "salmon", tx2gene = tx2gene)
cat("tximport     collapsed 10 transcripts to", nrow(txi$counts), "genes\n")

###########################################
# The method itself
###########################################

source("edge_conservation.r")
genes <- paste0("g", 1:3)
idx <- setNames(seq_along(genes), genes)
adj <- build_adjacency(
    data.frame(gene1 = c("g1", "g2"), gene2 = c("g2", "g3"), weight = c(1, 1)),
    idx, 3)
cat("edge_conservation  adjacency symmetric:", isSymmetric(adj), "\n")

cat("\nsmoke test passed\n")
