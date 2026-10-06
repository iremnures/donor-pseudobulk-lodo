suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

sce <- readRDS("data/raw/kang18_8vs8.rds")
meta <- as.data.frame(colData(sce))
counts <- assay(sce, "counts")

stopifnot(
  inherits(counts, "sparseMatrix"),
  all(is.finite(counts@x)),
  all(counts@x >= 0),
  all(counts@x == round(counts@x))
)

meta$source_column <- seq_len(ncol(sce))
meta$donor_id <- as.character(meta$ind)
meta$condition <- as.character(meta$stim)
meta$cell_type <- as.character(meta$cell)
meta$sample_id <- paste(meta$donor_id, meta$condition, sep = "_")

singlet <- as.character(meta$multiplets) == "singlet"
named_type <- !is.na(meta$cell_type) & nzchar(meta$cell_type)
keep <- singlet & named_type

cat("\n--- Cell eligibility ---\n")
cat("Starting cells:", nrow(meta), "\n")
cat("Non-singlet labels:", sum(!singlet), "\n")
cat("Missing/empty cell type:", sum(!named_type), "\n")
cat("Excluded unique cells:", sum(!keep), "\n")
cat("Eligible cells:", sum(keep), "\n")

selected <- meta[keep, , drop = FALSE]
selected$total_counts <- Matrix::colSums(counts)[keep]

write.csv(
  selected,
  "results/tables/kang_eligible_cell_metadata.csv",
  row.names = FALSE
)

coverage <- aggregate(
  source_column ~ donor_id + sample_id + condition + cell_type,
  data = selected,
  FUN = length
)
names(coverage)[names(coverage) == "source_column"] <- "n_cells"

write.csv(
  coverage,
  "results/tables/kang_cell_type_sample_coverage.csv",
  row.names = FALSE
)

cat("\n--- CD14+ monocytes: donor x condition ---\n")
monocytes <- coverage[
  coverage$cell_type == "CD14+ Monocytes", , drop = FALSE
]
print(xtabs(n_cells ~ donor_id + condition, data = monocytes))

cat("\n--- Coverage summary by cell type ---\n")
summary <- do.call(rbind, lapply(
  split(coverage, coverage$cell_type),
  function(tab) {
    data.frame(
      cell_type = tab$cell_type[1],
      observed_samples = nrow(tab),
      minimum_cells = min(tab$n_cells),
      median_cells = median(tab$n_cells),
      maximum_cells = max(tab$n_cells)
    )
  }
))
rownames(summary) <- NULL
print(summary, row.names = FALSE)

write.csv(
  summary,
  "results/tables/kang_coverage_summary.csv",
  row.names = FALSE
)

cat("\nCoverage audit saved. No pseudobulk or statistical tests yet.\n")
