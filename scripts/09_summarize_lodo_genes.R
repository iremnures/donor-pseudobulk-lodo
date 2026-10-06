primary <- read.csv(
  "results/tables/cd14_monocyte_primary_deseq2.csv",
  stringsAsFactors = FALSE
)
long <- read.csv(
  "results/tables/lodo/all_gene_results.csv",
  stringsAsFactors = FALSE
)

donors <- unique(long$excluded_donor)
stopifnot(length(donors) == 8, !anyDuplicated(primary$ENSEMBL))

tables <- lapply(donors, function(donor) {
  tab <- long[long$excluded_donor == donor, , drop = FALSE]
  stopifnot(
    nrow(tab) == nrow(primary),
    !anyDuplicated(tab$ENSEMBL),
    setequal(tab$ENSEMBL, primary$ENSEMBL)
  )
  tab[match(primary$ENSEMBL, tab$ENSEMBL), , drop = FALSE]
})

lfc <- do.call(cbind, lapply(tables, function(tab) tab$log2FoldChange))
se <- do.call(cbind, lapply(tables, function(tab) tab$lfcSE))
padj <- do.call(cbind, lapply(tables, function(tab) tab$padj))
colnames(lfc) <- colnames(se) <- colnames(padj) <- as.character(donors)

difference <- sweep(lfc, 1, primary$log2FoldChange, "-")
same_direction <- sweep(sign(lfc), 1, sign(primary$log2FoldChange), "==")
significant <- is.finite(padj) & padj < 0.05

summary <- primary
summary$LODO_same_direction <- rowSums(same_direction, na.rm = TRUE)
summary$LODO_FDR_below_0.05 <- rowSums(significant)
summary$LODO_direction_and_FDR <- rowSums(
  same_direction & significant, na.rm = TRUE
)
summary$LODO_missing_LFC <- rowSums(!is.finite(lfc))
summary$LODO_missing_padj <- rowSums(!is.finite(padj))

finite_range <- function(values, fun) {
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  fun(values)
}
summary$LODO_min_log2FC <- apply(lfc, 1, finite_range, fun = min)
summary$LODO_max_log2FC <- apply(lfc, 1, finite_range, fun = max)
summary$maximum_absolute_LFC_change <- apply(
  abs(difference), 1, finite_range, fun = max
)
summary$maximum_lfcSE <- apply(se, 1, finite_range, fun = max)

summary$donor_for_largest_LFC_change <- apply(
  abs(difference), 1, function(values) {
    if (!any(is.finite(values))) return(NA_character_)
    values[!is.finite(values)] <- -Inf
    colnames(difference)[which.max(values)]
  }
)

write.csv(
  summary,
  "results/tables/lodo/gene_stability_summary.csv",
  row.names = FALSE
)

primary_sig <- is.finite(summary$padj) & summary$padj < 0.05
selected <- summary[primary_sig, , drop = FALSE]

cat("\n--- Primary-significant genes: LODO FDR retention ---\n")
print(table(selected$LODO_FDR_below_0.05))

cat("\n--- Primary-significant genes: direction retention ---\n")
print(table(selected$LODO_same_direction))

columns <- c(
  "SYMBOL", "log2FoldChange", "lfcSE",
  "LODO_min_log2FC", "LODO_max_log2FC",
  "maximum_absolute_LFC_change",
  "donor_for_largest_LFC_change", "LODO_FDR_below_0.05"
)

cat("\n--- Largest LFC changes among primary-significant genes ---\n")
ordered <- selected[
  order(selected$maximum_absolute_LFC_change, decreasing = TRUE), ]
print(head(ordered[, columns], 15), row.names = FALSE)

cat("\n--- Previously plotted genes ---\n")
print(
  summary[
    summary$SYMBOL %in% c("CCL8", "ISG20", "IRF7", "RPL6"),
    columns
  ],
  row.names = FALSE
)

cat("\nGene-level sensitivity summary saved.\n")
