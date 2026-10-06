suppressPackageStartupMessages({
  library(muscData)
  library(SingleCellExperiment)
})

options(timeout = 1200)

path <- "data/raw/kang18_8vs8.rds"

if (file.exists(path)) {
  cat("Reading saved dataset...\n")
  sce <- readRDS(path)
} else {
  cat("Downloading/loading Kang18_8vs8...\n")
  sce <- muscData::Kang18_8vs8()
  saveRDS(sce, path)
  cat("Dataset saved.\n")
}

cat("\n--- Dataset dimensions: genes x cells ---\n")
print(dim(sce))

cat("\n--- Available assays ---\n")
print(assayNames(sce))

cat("\n--- Count matrix representation ---\n")
print(class(assay(sce, "counts")))
cat("Object size:", format(object.size(sce), units = "auto"), "\n")

meta <- as.data.frame(colData(sce))

cat("\n--- Metadata columns ---\n")
print(names(meta))

cat("\n--- First six cells ---\n")
print(head(meta))

required <- c("ind", "stim", "cell", "multiplets")
stopifnot(all(required %in% names(meta)))

cat("\n--- Donor assignments ---\n")
print(table(meta$ind, useNA = "ifany"))

cat("\n--- Donor x condition: cell counts ---\n")
print(table(meta$ind, meta$stim, useNA = "ifany"))

cat("\n--- Cell-type labels x condition ---\n")
print(table(meta$cell, meta$stim, useNA = "ifany"))

cat("\n--- Multiplet field ---\n")
print(table(meta$multiplets, useNA = "ifany"))

cat("\n--- Missing metadata values ---\n")
print(colSums(is.na(meta[required])))

meta$source_barcode <- colnames(sce)
meta$source_column <- seq_len(ncol(sce))

write.csv(
  meta,
  "results/tables/kang_cell_metadata_initial.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "docs/session_info_initial.txt"
)

cat("\nMetadata inspection saved. No cells removed.\n")
