# Batch-effect assessment — rat lung proteomics (reviewer comment 12)

Assesses whether sample-preparation day and LC-MS/MS injection order affect the
reported Aged vs Young differential abundance results.

## Folder structure

```
batch_effect_analysis/
├── README.md
├── code/
│   ├── 01_batch_effect_models.R     # per-protein models, LRT, PCA, variance partitioning, sensitivity analysis
│   ├── 02_supplementary_figure.R    # Supplementary Figure (panels a–d)
│   └── 03_supplementary_table.R     # Supplementary Table (sheets A, B, C1, C2)
├── data/                            # inputs (copied, not modified)
│   ├── UW.Madison.rat.lung.proteomics_data_tables.RData   # Script_3 output; data_imp = 3,600 proteins x 96 samples
│   ├── sample_metadata.xlsx                               # incl. Preparation Day and Injection Order
│   ├── Real_values_number_in_sample_group.xlsx            # observed-value counts (Method 4 filter)
│   └── 38 matrisome DEPs with unique features in 8 pairwise comparisons.txt
├── results/                         # output of 01
│   ├── batch_results.RData                       # cmp, tech, vpart, pc_assoc, robust, overall, workdata
│   ├── Batch_effect_sensitivity_analysis.xlsx    # all intermediate tables, one sheet each
│   ├── Batch_effect_assessment_figures.pdf       # diagnostic plots
│   └── 0X_*_workspace.RData                      # full workspace of each script (save.image)
├── figures/                         # output of 02
│   ├── Supplementary_Figure_batch_effect.pdf / .png
│   └── PCA_by_preparation_day.pdf                # Script_3-style PCA, one view per page
└── tables/                          # output of 03
    └── Supplementary_Table_batch_effect.xlsx
```

## How to run

Run from this folder (the one containing `code/` and `data/`):

```bash
Rscript code/01_batch_effect_models.R      # ~20–30 min (3,600 proteins x 2 models, Kenward–Roger df)
Rscript code/02_supplementary_figure.R
Rscript code/03_supplementary_table.R
```

R packages: lme4, lmerTest, emmeans, pbkrtest, dplyr, tidyr, readxl, openxlsx,
ggplot2, patchwork, ggbiplot (vqv/ggbiplot, as in Script_3).

## Models

- **Original (Model 0):** `Protein ~ Age + Sex + Region + Age:Sex + Age:Region + Sex:Region + (1 | PID)`
  — identical to Script_6; reproduces the reported 1,448 DEP calls / 670 proteins.
- **Adjusted (Model 1):** Model 0 + `PreparationDay` (factor, 4 levels) + `InjectionOrder` (continuous, scaled).

DEPs: Method 4 filter (≥ 5 observed values in either group), BH-FDR within each
comparison, FDR < 0.05 and |log2FC| > log2(1.5). Contrasts Aged − Young via emmeans
(Kenward–Roger df, as in Script_6).

## Analyses

1. **Design balance** — χ² / Fisher (biological factor × preparation day), Kruskal–Wallis (injection order by group).
2. **PCA** — `prcomp(t(data_imp))`, association of PC1–10 with each factor.
3. **Likelihood-ratio test** — Model 0 vs Model 1 (ML, 4 df) per protein, BH-adjusted.
4. **Variance partitioning** — per-protein fraction of variance by factor (variancePartition-style model).
5. **Sensitivity analysis** — Aged vs Young log2FC, FDR and DEP calls compared between models.

## Key results

| | |
|---|---|
| Design balance (preparation day) | Age p = 0.51, Sex p = 0.23, Region p = 0.53 (χ²) |
| Design balance (injection order) | Age p = 0.90, Sex p = 0.12, Region p = 0.74 (Kruskal–Wallis) |
| LRT, technical covariates | 2,626 / 3,600 proteins FDR < 0.05 |
| Variance explained (median) | preparation day 7.3%, injection order 1.9% |
| Age-effect concordance | Pearson r = 0.96–0.999 across 8 comparisons (0.984 pooled) |
| DEP calls retained | 1,336 / 1,448 (92.3%); 100% same direction |
| Matrisome (38 proteins) | 37 / 38 remain DEP in ≥ 1 comparison |
