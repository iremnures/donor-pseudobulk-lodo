suppressPackageStartupMessages(library(DESeq2))

input <- readRDS(
  "data/processed/cd14_monocyte_pseudobulk_filtered.rds"
)

# Only the small genes x 16-sample pseudobulk matrix is made dense.
count_matrix <- as.matrix(input$counts)
stopifnot(
  all(is.finite(count_matrix)),
  all(count_matrix >= 0),
  all(count_matrix == round(count_matrix)),
  max(count_matrix) <= .Machine$integer.max
)
storage.mode(count_matrix) <- "integer"

samples <- input$samples
stopifnot(identical(colnames(count_matrix), rownames(samples)))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = samples,
  design = ~ donor_id + condition
)

dds <- DESeq(dds, parallel = FALSE)

res <- results(
  dds,
  contrast = c("condition", "stim", "ctrl"),
  alpha = 0.05,
  pAdjustMethod = "BH"
)

tab <- as.data.frame(res)
tab$ENSEMBL <- rownames(tab)
tab$SYMBOL <- input$genes$SYMBOL[
  match(tab$ENSEMBL, input$genes$ENSEMBL)
]

tab <- tab[, c(
  "ENSEMBL", "SYMBOL", "baseMean", "log2FoldChange",
  "lfcSE", "stat", "pvalue", "padj"
)]

write.csv(
  tab,
  "results/tables/cd14_monocyte_primary_deseq2.csv",
  row.names = FALSE
)

saveRDS(
  dds,
  "data/processed/cd14_monocyte_primary_dds.rds"
)

sample_summary <- as.data.frame(colData(dds))
write.csv(
  sample_summary,
  "results/tables/cd14_monocyte_primary_sample_metadata.csv",
  row.names = FALSE
)

significant <- is.finite(tab$padj) & tab$padj < 0.05

cat("\n--- Primary model summary ---\n")
cat("Genes in model:", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n")
cat("Donors:", nlevels(dds$donor_id), "\n")
cat("Finite p-values:", sum(is.finite(tab$pvalue)), "\n")
cat("Missing adjusted p-values:", sum(is.na(tab$padj)), "\n")
cat("FDR < 0.05:", sum(significant), "\n")
cat("Stim-direction:", sum(significant & tab$log2FoldChange > 0), "\n")
cat("Ctrl-direction:", sum(significant & tab$log2FoldChange < 0), "\n")

cat("\n--- Size factors ---\n")
print(data.frame(
  sample_id = colnames(dds),
  size_factor = sizeFactors(dds)
), row.names = FALSE)

cat("\n--- Top 10 by adjusted p-value ---\n")
ordered <- tab[order(tab$padj, na.last = TRUE), ]
print(head(ordered, 10), row.names = FALSE)

capture.output(
  sessionInfo(),
  file = "docs/session_info_primary_model.txt"
)

cat("\nPrimary model saved; sample-level QC and LODO remain pending.\n")
