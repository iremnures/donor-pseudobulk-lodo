suppressPackageStartupMessages(library(DESeq2))

run_dir <- trimws(readLines("repro_checks/latest_run.txt", n = 1))
stopifnot(
  dir.exists(run_dir),
  file.exists(file.path(run_dir, "execution_completed.txt"))
)

checks <- list()

record <- function(label, original, reproduced) {
  comparison <- all.equal(
    original, reproduced,
    tolerance = 1e-8,
    check.attributes = TRUE
  )
  passed <- isTRUE(comparison)

  checks[[length(checks) + 1]] <<- data.frame(
    check = label,
    passed = passed,
    detail = if (passed) "Match" else paste(comparison, collapse = "; ")
  )
}

raw_path <- "data/raw/kang18_8vs8.rds"
record(
  "Starting dataset fingerprint",
  unname(tools::md5sum(raw_path)),
  unname(tools::md5sum(file.path(run_dir, raw_path)))
)

csv_paths <- list.files(
  file.path(run_dir, "results/tables"),
  pattern = "\\.csv$",
  recursive = TRUE
)
stopifnot(length(csv_paths) > 0)

for (relative in csv_paths) {
  path <- file.path("results/tables", relative)
  if (!file.exists(path)) {
    record(paste("Table:", relative), "Original file exists", "Missing")
    next
  }

  original <- read.csv(path, check.names = FALSE)
  reproduced <- read.csv(
    file.path(run_dir, path), check.names = FALSE
  )
  record(paste("Table:", relative), original, reproduced)
}

model_paths <- c(
  "data/processed/cd14_monocyte_primary_dds.rds",
  file.path(
    "data/processed/lodo",
    paste0("no_", c("101", "107", "1015", "1016",
                   "1039", "1244", "1256", "1488"), "_dds.rds")
  )
)

for (path in model_paths) {
  stopifnot(
    file.exists(path),
    file.exists(file.path(run_dir, path))
  )

  original <- readRDS(path)
  reproduced <- readRDS(file.path(run_dir, path))
  label <- basename(path)

  record(paste(label, "raw counts"),
         counts(original), counts(reproduced))
  record(paste(label, "sample metadata"),
         as.data.frame(colData(original)),
         as.data.frame(colData(reproduced)))
  record(paste(label, "design"),
         design(original), design(reproduced))
  record(paste(label, "size factors"),
         sizeFactors(original), sizeFactors(reproduced))
  record(paste(label, "dispersions"),
         dispersions(original), dispersions(reproduced))

  contrast <- c("condition", "stim", "ctrl")
  record(
    paste(label, "DESeq2 results"),
    as.data.frame(results(original, contrast = contrast, alpha = 0.05)),
    as.data.frame(results(reproduced, contrast = contrast, alpha = 0.05))
  )
}

summary <- do.call(rbind, checks)
write.csv(
  summary,
  "docs/reproducibility_comparison.csv",
  row.names = FALSE
)

cat("\n--- Reproduction comparison ---\n")
print(summary[, c("check", "passed")], row.names = FALSE)

if (any(!summary$passed)) {
  cat("\n--- Differences requiring inspection ---\n")
  print(summary[!summary$passed, ], row.names = FALSE)
  stop("Some comparisons failed; inspect before claiming agreement.")
}

cat("\nALL COMPARISONS PASSED.\n")
cat("Numerical tolerance: 1e-8.\n")
cat("Scope: same saved input and current R/package environment.\n")
