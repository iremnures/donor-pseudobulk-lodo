# Donor-Level Pseudobulk RNA-seq and LODO Sensitivity

An educational reanalysis of IFN-beta-treatment-associated expression
differences in CD14+ monocytes from eight paired donors.

## Research question

Are the estimated expression differences consistent across donors,
and how do they change when one donor is excluded?

## Dataset and biological units

- Source: Kang et al., accessed through muscData::Kang18_8vs8().
- Original study: https://doi.org/10.1038/nbt.4042
- Eight lupus donors, each contributing control and stimulated cells.
- Primary analysis: 16 donor-condition CD14+ monocyte profiles.
- Biological replication is defined by donors, not individual cells.

Source singlet classifications and available cell-type labels were
used to select eligible cells. Each primary profile contained
126-779 CD14+ monocytes.

## Analysis

Raw counts were summed within donor, condition and cell type.
Control and stimulated profiles from the same donor remained separate.

Genes with counts >= 10 in at least eight of 16 profiles were retained.

DESeq2 used:

    ~ donor_id + condition

The contrast was stimulated versus control. Positive log2 fold
changes indicate stimulated direction.

Benjamini-Hochberg FDR correction was applied within the primary
CD14+ monocyte gene-testing family.

Eight leave-one-donor-out (LODO) models excluded both conditions
from one donor at a time. The primary gene universe was retained,
while normalization and dispersion parameters were re-estimated.

## Main results

- 3,707 genes tested.
- 2,334 genes with primary FDR < 0.05.
- 1,183 stimulated-direction and 1,151 control-direction genes.
- All 2,334 primary-significant genes retained direction in all
  eight LODO models.
- 1,977 also retained FDR < 0.05 in every LODO model.

JUP and AXL retained direction and significance while showing
sensitivity in effect magnitude. MRPS16 retained direction but
lost FDR significance in every donor-exclusion model.

## Effect magnitude and uncertainty

![Primary and LODO effect estimates](results/figures/cd14_sensitive_genes_lodo_estimates.png)

Intervals are approximate, unadjusted 95% Wald intervals.
They are not multiple-testing-adjusted confidence intervals.

## Interpretation limits

- LODO is internal sensitivity analysis, not independent validation.
- Control and stimulated cells were pooled separately in the source
  experiment; pool-specific technical effects cannot be separated
  from condition by the donor-blocked model.
- Pseudobulk does not establish that every cell responds identically.
- Normalized expression changes are relative, not absolute RNA
  changes per cell.
- Source cell-type annotations may affect the aggregated profiles.
- Selected gene examples are exploratory.
- This analysis does not test cell-type abundance or lupus versus
  healthy differences.

## Environment and execution

The verified environment used R 4.5.3 and DESeq2 1.50.2.
All recorded package versions are in docs/software_versions.csv.

Required R packages:
muscData, SingleCellExperiment, SummarizedExperiment, Matrix,
DESeq2 and ggplot2.

From the project root, create the output directories and run:

    mkdir -p data/raw data/processed results/tables results/figures docs logs
    Rscript scripts/check_environment.R

Then run numbered R scripts 01-10 in order.
Script 01 downloads/loads the source dataset if the saved input
is absent. Internet access is required for initial retrieval.

The Python toy example is separate from the real-data workflow
and requires pandas.

## Reproducibility verification

The numbered analysis workflow was rerun in a separate directory.
Generated CSV tables and primary/LODO model comparisons passed
with numerical tolerance 1e-8.

This check used the same saved input and current software environment.
Fresh installation and download were not tested.

Raw datasets, processed model objects, logs and reproduction
directories are excluded from Git.

## Documentation

- [Analysis plan](docs/analysis_plan.md)
- [Results and interpretation](docs/m09_results_report.md)
- [Reproducibility verification](docs/reproducibility.md)
- [Detailed comparison checks](docs/reproducibility_comparison.csv)
- [Software versions](docs/software_versions.csv)
