suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- readRDS("data/processed/cd14_monocyte_primary_dds.rds")
results_tab <- read.csv(
  "results/tables/cd14_monocyte_primary_deseq2.csv",
  stringsAsFactors = FALSE
)

selected <- c("CCL8", "ISG20", "IRF7", "RPL6")

mapping <- results_tab[
  results_tab$SYMBOL %in% selected,
  c("ENSEMBL", "SYMBOL"),
  drop = FALSE
]

stopifnot(
  nrow(mapping) == length(selected),
  !anyDuplicated(mapping$SYMBOL),
  all(mapping$ENSEMBL %in% rownames(dds))
)

normalized <- counts(dds, normalized = TRUE)
raw <- counts(dds, normalized = FALSE)
meta <- as.data.frame(colData(dds))

tables <- lapply(seq_len(nrow(mapping)), function(i) {
  gene_id <- mapping$ENSEMBL[i]
  data.frame(
    gene_id = gene_id,
    gene_name = mapping$SYMBOL[i],
    sample_id = colnames(dds),
    donor_id = as.character(meta$donor_id),
    condition = as.character(meta$condition),
    raw_count = as.numeric(raw[gene_id, ]),
    normalized_count = as.numeric(normalized[gene_id, ])
  )
})

tab <- do.call(rbind, tables)
tab$condition <- factor(tab$condition, levels = c("ctrl", "stim"))
tab$gene_name <- factor(tab$gene_name, levels = selected)
tab$plot_value <- log2(tab$normalized_count + 1)

write.csv(
  tab,
  "results/tables/cd14_selected_genes_paired_counts.csv",
  row.names = FALSE
)

p <- ggplot(
  tab,
  aes(condition, plot_value, group = donor_id)
) +
  geom_line(color = "grey65") +
  geom_point(aes(color = condition), size = 2.5) +
  geom_text(
    data = tab[tab$condition == "stim", ],
    aes(label = donor_id),
    hjust = -0.2, size = 2.7, check_overlap = TRUE
  ) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 2) +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.7))) +
  labs(
    title = "Selected genes across paired donors",
    subtitle = "Each line connects the same donor",
    x = "Condition",
    y = "log2(normalized count + 1)"
  ) +
  theme_bw()

ggsave(
  "results/figures/cd14_selected_genes_paired.png",
  p, width = 9, height = 7, dpi = 180
)

summary <- do.call(rbind, lapply(selected, function(gene) {
  sub <- tab[as.character(tab$gene_name) == gene, ]
  control <- sub[sub$condition == "ctrl", ]
  stimulated <- sub[sub$condition == "stim", ]
  stimulated <- stimulated[
    match(control$donor_id, stimulated$donor_id), ]

  difference <- (
    stimulated$normalized_count - control$normalized_count
  )

  data.frame(
    gene_name = gene,
    stim_higher_pairs = sum(difference > 0),
    ctrl_higher_pairs = sum(difference < 0),
    equal_pairs = sum(difference == 0),
    zero_raw_count_samples = sum(sub$raw_count == 0)
  )
}))

write.csv(
  summary,
  "results/tables/cd14_selected_genes_direction_summary.csv",
  row.names = FALSE
)

cat("\n--- Descriptive paired directions ---\n")
print(summary, row.names = FALSE)
cat("\nPaired figure and count tables saved.\n")
