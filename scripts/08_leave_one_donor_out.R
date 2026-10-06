suppressPackageStartupMessages(library(DESeq2))

primary <- readRDS(
  "data/processed/cd14_monocyte_primary_dds.rds"
)
primary_tab <- read.csv(
  "results/tables/cd14_monocyte_primary_deseq2.csv",
  stringsAsFactors = FALSE
)
primary_tab <- primary_tab[
  match(rownames(primary), primary_tab$ENSEMBL), , drop = FALSE
]
stopifnot(identical(primary_tab$ENSEMBL, rownames(primary)))

raw_counts <- counts(primary, normalized = FALSE)
meta <- as.data.frame(colData(primary))
# Remove inherited normalization before creating fresh models.
meta$sizeFactor <- NULL
donors <- levels(meta$donor_id)

dir.create("data/processed/lodo", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/lodo", recursive = TRUE, showWarnings = FALSE)

all_results <- list()
run_summary <- list()

for (donor in donors) {
  cat("\n--- Excluding donor", donor, "---\n", flush = TRUE)

  keep <- as.character(meta$donor_id) != donor
  reduced_meta <- droplevels(meta[keep, , drop = FALSE])

  pairs <- table(reduced_meta$donor_id, reduced_meta$condition)
  stopifnot(all(pairs == 1), nrow(pairs) == 7)

  # Fresh object: no inherited normalization or dispersion estimates.
  reduced <- DESeqDataSetFromMatrix(
    countData = raw_counts[, keep, drop = FALSE],
    colData = reduced_meta,
    design = ~ donor_id + condition
  )
  stopifnot(is.null(sizeFactors(reduced)))

  design_matrix <- model.matrix(
    ~ donor_id + condition,
    data = reduced_meta
  )
  stopifnot(qr(design_matrix)$rank == ncol(design_matrix))

  reduced <- DESeq(reduced, parallel = FALSE)
  res <- results(
    reduced,
    contrast = c("condition", "stim", "ctrl"),
    alpha = 0.05,
    pAdjustMethod = "BH"
  )

  stopifnot(identical(rownames(reduced), rownames(primary)))

  tab <- as.data.frame(res)
  tab$ENSEMBL <- rownames(tab)
  tab$SYMBOL <- primary_tab$SYMBOL
  tab$excluded_donor <- donor
  tab$primary_log2FC <- primary_tab$log2FoldChange
  tab$primary_lfcSE <- primary_tab$lfcSE
  tab$primary_padj <- primary_tab$padj
  tab$log2FC_difference <- tab$log2FoldChange - tab$primary_log2FC

  finite <- is.finite(tab$log2FoldChange) &
    is.finite(tab$primary_log2FC)
  tab$same_direction <- NA
  tab$same_direction[finite] <- (
    sign(tab$log2FoldChange[finite]) ==
      sign(tab$primary_log2FC[finite])
  )

  significant <- is.finite(tab$padj) & tab$padj < 0.05
  primary_sig <- is.finite(tab$primary_padj) &
    tab$primary_padj < 0.05

  write.csv(
    tab,
    paste0("results/tables/lodo/no_", donor, "_results.csv"),
    row.names = FALSE
  )
  saveRDS(
    reduced,
    paste0("data/processed/lodo/no_", donor, "_dds.rds")
  )

  write.csv(
    data.frame(
      sample_id = colnames(reduced),
      size_factor = sizeFactors(reduced)
    ),
    paste0("results/tables/lodo/no_", donor, "_size_factors.csv"),
    row.names = FALSE
  )

  run_summary[[donor]] <- data.frame(
    excluded_donor = donor,
    donors = nlevels(reduced$donor_id),
    samples = ncol(reduced),
    genes = nrow(reduced),
    genes_FDR_below_0.05 = sum(significant),
    primary_sig_same_direction = sum(
      primary_sig & tab$same_direction, na.rm = TRUE
    ),
    primary_sig_retaining_FDR = sum(primary_sig & significant),
    primary_sig_direction_and_FDR = sum(
      primary_sig & significant & tab$same_direction,
      na.rm = TRUE
    ),
    median_absolute_LFC_change_primary_sig = median(
      abs(tab$log2FC_difference[primary_sig]), na.rm = TRUE
    ),
    median_SE_ratio_primary_sig = median(
      tab$lfcSE[primary_sig] / tab$primary_lfcSE[primary_sig],
      na.rm = TRUE
    ),
    missing_LFC = sum(!is.finite(tab$log2FoldChange)),
    missing_padj = sum(is.na(tab$padj))
  )

  all_results[[donor]] <- tab
  cat("Completed donor", donor, "\n", flush = TRUE)
}

summary <- do.call(rbind, run_summary)
rownames(summary) <- NULL

write.csv(
  summary,
  "results/tables/lodo/run_summary.csv",
  row.names = FALSE
)
write.csv(
  do.call(rbind, all_results),
  "results/tables/lodo/all_gene_results.csv",
  row.names = FALSE
)

cat("\n--- LODO run summary ---\n")
print(summary, row.names = FALSE)
cat("\nAll eight LODO models completed.\n")
