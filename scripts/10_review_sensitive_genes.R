suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

stability <- read.csv(
  "results/tables/lodo/gene_stability_summary.csv",
  stringsAsFactors = FALSE
)
long <- read.csv(
  "results/tables/lodo/all_gene_results.csv",
  stringsAsFactors = FALSE
)
dds <- readRDS("data/processed/cd14_monocyte_primary_dds.rds")

fragile <- stability[
  is.finite(stability$padj) &
    stability$padj < 0.05 &
    stability$LODO_FDR_below_0.05 == 0,
  , drop = FALSE
]
stopifnot(nrow(fragile) == 1)

selected <- stability[
  stability$SYMBOL %in% c("JUP", "AXL"), , drop = FALSE
]
stopifnot(nrow(selected) == 2)
selected <- rbind(selected, fragile)

selected$display_label <- ifelse(
  is.na(selected$SYMBOL) | selected$SYMBOL == "",
  selected$ENSEMBL,
  selected$SYMBOL
)

cat("\n--- Selected sensitivity examples ---\n")
print(selected[, c(
  "ENSEMBL", "SYMBOL", "log2FoldChange", "lfcSE",
  "padj", "LODO_FDR_below_0.05",
  "maximum_absolute_LFC_change"
)], row.names = FALSE)

write.csv(
  selected,
  "results/tables/lodo/selected_sensitivity_examples.csv",
  row.names = FALSE
)

normalized <- counts(dds, normalized = TRUE)
raw <- counts(dds, normalized = FALSE)
meta <- as.data.frame(colData(dds))

expression <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
  id <- selected$ENSEMBL[i]
  data.frame(
    gene_id = id,
    gene = selected$display_label[i],
    donor = as.character(meta$donor_id),
    condition = factor(meta$condition, levels = c("ctrl", "stim")),
    raw_count = as.numeric(raw[id, ]),
    normalized_count = as.numeric(normalized[id, ])
  )
}))

write.csv(
  expression,
  "results/tables/lodo/sensitivity_examples_paired_counts.csv",
  row.names = FALSE
)

p <- ggplot(
  expression,
  aes(condition, log2(normalized_count + 1), group = donor)
) +
  geom_line(aes(color = donor)) +
  geom_point(aes(color = donor), size = 2) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  labs(
    title = "Donor-level expression in sensitivity examples",
    x = "Condition", y = "log2(normalized count + 1)"
  ) +
  theme_bw()

ggsave(
  "results/figures/cd14_sensitive_genes_paired.png",
  p, width = 11, height = 5, dpi = 180
)

estimates <- long[
  long$ENSEMBL %in% selected$ENSEMBL,
  c("ENSEMBL", "excluded_donor", "log2FoldChange", "lfcSE", "padj")
]
estimates$model <- paste0("Without ", estimates$excluded_donor)

primary_rows <- data.frame(
  ENSEMBL = selected$ENSEMBL,
  excluded_donor = NA,
  log2FoldChange = selected$log2FoldChange,
  lfcSE = selected$lfcSE,
  padj = selected$padj,
  model = "Primary"
)
estimates <- rbind(primary_rows, estimates)

estimates$gene <- selected$display_label[
  match(estimates$ENSEMBL, selected$ENSEMBL)
]
model_order <- c(
  "Primary", paste0("Without ", unique(long$excluded_donor))
)
estimates$model <- factor(
  estimates$model, levels = rev(model_order)
)
estimates$lower <- estimates$log2FoldChange - 1.96 * estimates$lfcSE
estimates$upper <- estimates$log2FoldChange + 1.96 * estimates$lfcSE
estimates$FDR_status <- ifelse(
  !is.finite(estimates$padj), "Unavailable",
  ifelse(estimates$padj < 0.05, "FDR < 0.05", "FDR >= 0.05")
)

write.csv(
  estimates,
  "results/tables/lodo/sensitivity_examples_model_estimates.csv",
  row.names = FALSE
)

q <- ggplot(estimates, aes(log2FoldChange, model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_segment(aes(x = lower, xend = upper, yend = model)) +
  geom_point(aes(color = FDR_status), size = 2.5) +
  facet_wrap(~ gene, scales = "free_x", ncol = 3) +
  labs(
    title = "Primary and leave-one-donor-out estimates",
    subtitle = "Intervals: approximate 95% Wald confidence intervals",
    x = "Model log2 fold change: stim versus ctrl",
    y = NULL
  ) +
  theme_bw()

ggsave(
  "results/figures/cd14_sensitive_genes_lodo_estimates.png",
  q, width = 12, height = 6, dpi = 180
)

cat("\nPaired counts and model-estimate figures saved.\n")
