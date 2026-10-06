suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

sce <- readRDS("data/raw/kang18_8vs8.rds")

meta <- read.csv(
  "results/tables/kang_eligible_cell_metadata.csv",
  stringsAsFactors = FALSE
)

meta <- meta[
  meta$cell_type == "CD14+ Monocytes", , drop = FALSE
]

stopifnot(
  nrow(meta) > 0,
  !anyDuplicated(meta$source_column),
  all(meta$source_column %in% seq_len(ncol(sce)))
)

source_meta <- as.data.frame(colData(sce))
stopifnot(
  all(as.character(source_meta$ind[meta$source_column]) == meta$donor_id),
  all(as.character(source_meta$stim[meta$source_column]) == meta$condition),
  all(as.character(source_meta$cell[meta$source_column]) == meta$cell_type),
  all(as.character(source_meta$multiplets[meta$source_column]) == "singlet")
)

samples <- unique(
  meta[, c("sample_id", "donor_id", "condition", "cell_type")]
)
stopifnot(!anyDuplicated(samples$sample_id))

samples <- samples[
  order(samples$donor_id, samples$condition), , drop = FALSE
]
rownames(samples) <- samples$sample_id

group_index <- match(meta$sample_id, samples$sample_id)
stopifnot(!anyNA(group_index))

# Each selected cell contributes to exactly one sample.
membership <- sparseMatrix(
  i = seq_len(nrow(meta)),
  j = group_index,
  x = 1,
  dims = c(nrow(meta), nrow(samples))
)

cell_counts <- assay(sce, "counts")[
  , meta$source_column, drop = FALSE
]

pseudobulk <- cell_counts %*% membership
colnames(pseudobulk) <- samples$sample_id

samples$n_cells <- tabulate(group_index, nbins = nrow(samples))
samples$total_counts <- as.numeric(Matrix::colSums(pseudobulk))

stopifnot(
  all(samples$n_cells >= 50),
  nrow(samples) == 16,
  identical(colnames(pseudobulk), rownames(samples)),
  isTRUE(all.equal(
    as.numeric(sum(cell_counts)),
    as.numeric(sum(pseudobulk)),
    tolerance = 0
  ))
)

pairs <- table(samples$donor_id, samples$condition)
stopifnot(
  identical(sort(colnames(pairs)), c("ctrl", "stim")),
  all(pairs == 1)
)

annotation <- as.data.frame(rowData(sce))
annotation$source_gene_name <- rownames(sce)
annotation$source_row <- seq_len(nrow(sce))

saveRDS(
  list(
    counts = pseudobulk,
    samples = samples,
    genes = annotation
  ),
  "data/processed/cd14_monocyte_pseudobulk_unfiltered.rds"
)

write.csv(
  samples,
  "results/tables/cd14_monocyte_sample_metadata.csv",
  row.names = FALSE
)

write.csv(
  annotation,
  "results/tables/kang_gene_annotation.csv",
  row.names = FALSE
)

cat("\n--- Pseudobulk dimensions: genes x samples ---\n")
print(dim(pseudobulk))

cat("\n--- Sample metadata ---\n")
print(samples, row.names = FALSE)

cat("\n--- Paired design ---\n")
print(pairs)

cat("\n--- Gene annotation: first six rows ---\n")
print(head(annotation))

cat("\nTotal counts preserved: TRUE\n")
cat("Saved unfiltered pseudobulk. No DE testing performed.\n")
