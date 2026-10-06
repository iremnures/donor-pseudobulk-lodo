input <- readRDS(
  "data/processed/cd14_monocyte_pseudobulk_unfiltered.rds"
)

counts <- input$counts
samples <- input$samples
genes <- input$genes

stopifnot(
  identical(colnames(counts), rownames(samples)),
  identical(rownames(counts), genes$source_gene_name),
  nrow(counts) == nrow(genes)
)

gene_ids <- as.character(genes$ENSEMBL)
stopifnot(
  !anyNA(gene_ids),
  all(nzchar(gene_ids)),
  !anyDuplicated(gene_ids)
)
rownames(counts) <- gene_ids

# Use the prespecified primary filtering rule.
keep <- Matrix::rowSums(counts >= 10) >= 8

filtered_counts <- counts[keep, , drop = FALSE]
filtered_genes <- genes[keep, , drop = FALSE]
rownames(filtered_genes) <- gene_ids[keep]

samples$donor_id <- factor(samples$donor_id)
samples$condition <- factor(
  samples$condition,
  levels = c("ctrl", "stim")
)

stopifnot(!anyNA(samples$condition))

design <- model.matrix(
  ~ donor_id + condition,
  data = samples
)

stopifnot(
  qr(design)$rank == ncol(design),
  nrow(design) > ncol(design)
)

saveRDS(
  list(
    counts = filtered_counts,
    samples = samples,
    genes = filtered_genes,
    design_matrix = design,
    gene_filter = "count >=10 in at least 8 primary samples"
  ),
  "data/processed/cd14_monocyte_pseudobulk_filtered.rds"
)

write.csv(
  data.frame(
    ENSEMBL = gene_ids,
    SYMBOL = genes$SYMBOL,
    samples_count_at_least_10 = as.integer(
      Matrix::rowSums(counts >= 10)
    ),
    retained = as.logical(keep)
  ),
  "results/tables/cd14_monocyte_gene_filter.csv",
  row.names = FALSE
)

cat("\n--- Gene filtering ---\n")
cat("Before:", nrow(counts), "\n")
cat("Retained:", nrow(filtered_counts), "\n")
cat("Removed:", sum(!keep), "\n")

cat("\n--- Filtered dimensions: genes x samples ---\n")
print(dim(filtered_counts))

cat("\n--- Model design ---\n")
print(design)

cat("\nDesign columns:", ncol(design), "\n")
cat("Design rank:", qr(design)$rank, "\n")
cat("Full rank: TRUE\n")
cat("Reference condition: ctrl\n")
cat("Contrast direction: stim versus ctrl\n")

cat("\nFiltered input saved. No differential-expression testing yet.\n")
