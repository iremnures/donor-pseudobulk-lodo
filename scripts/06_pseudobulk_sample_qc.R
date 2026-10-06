suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- readRDS("data/processed/cd14_monocyte_primary_dds.rds")
meta <- as.data.frame(colData(dds))

# Transformation for visualization only; DE used raw counts.
transformed <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(transformed)

variances <- apply(mat, 1, var)
selected <- order(variances, decreasing = TRUE)[
  seq_len(min(500, nrow(mat)))
]

pca <- prcomp(t(mat[selected, , drop = FALSE]))
variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)

plot_data <- meta
plot_data$PC1 <- pca$x[rownames(meta), 1]
plot_data$PC2 <- pca$x[rownames(meta), 2]

write.csv(
  plot_data,
  "results/tables/cd14_monocyte_pca_sample_coordinates.csv",
  row.names = FALSE
)

p <- ggplot(plot_data, aes(PC1, PC2)) +
  geom_line(aes(group = donor_id), color = "grey65") +
  geom_point(aes(color = condition), size = 3) +
  geom_text(
    aes(label = donor_id),
    vjust = -0.7, size = 3, check_overlap = TRUE
  ) +
  labs(
    title = "CD14+ monocyte pseudobulk PCA",
    subtitle = "Lines connect paired donor samples",
    x = sprintf("PC1 (%.1f%%)", variance_percent[1]),
    y = sprintf("PC2 (%.1f%%)", variance_percent[2])
  ) +
  theme_bw()

ggsave(
  "results/figures/cd14_monocyte_paired_pca.png",
  p, width = 8, height = 6, dpi = 180
)

q <- ggplot(meta, aes(n_cells, sizeFactor, color = condition)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = donor_id),
    vjust = -0.7, size = 3, check_overlap = TRUE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Cell coverage and normalization factors",
    x = "Eligible monocytes per sample (log scale)",
    y = "DESeq2 size factor (log scale)"
  ) +
  theme_bw()

ggsave(
  "results/figures/cd14_monocyte_coverage_size_factors.png",
  q, width = 8, height = 6, dpi = 180
)

res <- results(
  dds,
  contrast = c("condition", "stim", "ctrl"),
  alpha = 0.05
)

png(
  "results/figures/cd14_monocyte_primary_ma.png",
  width = 1400, height = 1000, res = 160
)
plotMA(
  res,
  alpha = 0.05,
  main = "CD14+ monocytes: stim versus ctrl"
)
dev.off()

cat("\n--- PCA variance ---\n")
cat("PC1:", round(variance_percent[1], 2), "%\n")
cat("PC2:", round(variance_percent[2], 2), "%\n")

cat("\nThree sample-level diagnostic figures saved.\n")
cat("No sample exclusions or model changes performed.\n")
