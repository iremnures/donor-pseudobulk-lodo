# M09 primary analysis plan

## Question
How does IFN-beta stimulation change expression within
CD14+ monocytes across paired donors?

## Dataset and analysis unit
Kang18_8vs8, accessed through muscData.
Use source singlet labels and nonmissing cell-type annotations.
Aggregate raw counts by biological sample and cell type.
Keep donor-condition samples separate.

Eight donors provide 16 paired samples.
Cells are not independent donor-level replicates.

## Primary scope
Cell type: CD14+ Monocytes.
Contrast: stim versus ctrl.
Minimum coverage: 50 eligible monocytes per sample.
This study-specific threshold was set after coverage inspection
but before differential-expression results.

Retain genes with raw pseudobulk counts >=10 in at least
8 of the 16 primary samples.

## Model and testing
DESeq2 design: ~ donor_id + condition.
Reference condition: ctrl.
Positive log2 fold change: higher expression in stim.

Primary test family: retained genes in CD14+ monocytes.
Apply Benjamini-Hochberg FDR correction, threshold 0.05.
Other cell types are outside the current primary test family.

## Leave-one-donor-out sensitivity
Exclude both condition samples of one donor in each run.
Keep the primary filtered gene universe.
Create fresh models and re-estimate normalization and dispersions.
Review effect direction, magnitude, standard error and FDR status.
FDR changes alone do not determine robustness.

LODO is within-cohort sensitivity analysis, not independent validation.

## Interpretation
The contrast concerns expression within an annotated cell type,
not changes in cell abundance or a lupus-versus-healthy comparison.
Source annotation uncertainty and pooled experimental processing
remain relevant limitations.

Document later deviations with reasons.
