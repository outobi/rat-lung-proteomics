# Rat lung proteomics: age, sex and regional variation

Analysis scripts, supporting data, pathway summaries and figures for the rat lung proteomics study. The repository covers protein/peptide detection, preprocessing, differential expression, pathway enrichment, age-by-sex and age-by-region interactions, cell-type interpretation and Shiny applications.

This is a manuscript analysis archive, not a single-command pipeline. Each numbered folder retains its own inputs and outputs. External annotation, enrichment and manual curation are separate steps.

## Sample terminology

- **PP:** peripubertal; **YA:** young adult.
- **F / M:** female / male.
- **B / C / P / T:** bottom / central / peripheral / top lung region.
- Legacy filenames and model labels use `Young` or `Y` for PP and `Aged` or `A` for YA. These are the same study groups, not additional age cohorts. Check sample metadata because character order varies among legacy codes.

## Workflow

![Analysis workflow from protein annotation and preprocessing to differential expression, enrichment, cell-type analysis and Shiny applications.](docs/workflow.svg)

Arrows describe data dependencies, not automatic file transfer. Script 2 is a coverage-summary branch, not a prerequisite for modeling. Script 9 uses separately prepared expression/metadata workbooks plus reference scRNA-seq data. The diagram uses multiline SVG text to avoid clipped labels in README viewers.

## Scripts and outputs

| Module | Purpose and principal files |
|---|---|
| [1 — Protein counts and annotation](Script_1_protein_number_count/) | `Script_1_protein_number_counts.R` reads the DIA-NN protein-group matrix, UniProt annotations and sample index. External matrisome annotation supplies `Matrisome_annotation_output.xlsx`. Outputs include `Matrisome_annotation_table.xlsx` and `Protein_counts_summary.xlsx`. |
| [2 — Peptide counts](Script_2_peptide_number_count/) | `Script_2_peptide_number_count.R` combines `Peptide.xlsx`, protein annotations and sample metadata to generate `Peptide_counts_summary.xlsx`. |
| [3 — Preprocessing and QC](Script_3_data_preprocessing/) | `Script_3_data_preprocessing.R` reads the `Raw_table_annotation` sheet of `Protein_counts_summary.xlsx` and sample metadata. Performs filtering, log2 transformation, normalization, imputation and QC/PCA. Exports processed intensity/data snapshots and `Real_values_number_in_sample_group.xlsx`. |
| [4 — Single-factor comparisons](Script_4_single_factor_DE/) | `Script_4_Single_factor_DE.R` estimates age, sex and regional contrasts and produces `Single_factor_DEP_summary.xlsx`, volcano plots, counts and heatmaps. It uses an additive age/sex/region mixed model with a subject random intercept; “single-factor” describes the reported contrasts, not an unadjusted one-predictor model. |
| [5 — Single-factor pathways](Script_5_single_factor_pathway/) | Subfolders contain YA-versus-PP ORA/GSEA, female-versus-male GSEA and regional GSEA. Scripts plot supplied GO/Reactome summaries; supporting workbooks include pathway annotations and leading-edge protein tables. |
| [6 — Multiple-factor comparisons](Script_6_multiple_factor_DE/) | `Script_6_multiple_factor_DE.R` fits the Method 4 mixed model with age, sex, region, pairwise interactions and an animal random intercept. Exports `age.vs.young4.xlsx`, `female.vs.male4.xlsx`, `region.vs.region4.xlsx`, DEP summaries and figures: 8 age, 8 sex and 24 regional contrasts. |
| [7 — Overlap and stratified ORA](Script_7_mulltiple_factor_pathway/) | Age-comparison protein sets stratified by sex/region, female-versus-male overlap analysis, and central-versus-bottom/peripheral ORA in PP females and males. Includes overlap plots, selected pathways and pathway-to-protein mappings. |
| [8 — Interaction analysis](Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/) | The `_DE.R` script exports age-by-sex and age-by-region interaction protein tables. The separate `GSEA/` script reads `gsea summary.xlsx` and plots curated interaction enrichment results. |
| [9 — Regional cell-type analysis](Script_9_region_cell_type_enrichment_with_scRNAseq/) | Two branches: regional protein markers mapped to reference scRNA-seq cell types, and scRNA-seq-derived cell-type marker scores across proteomic regions. Inputs include reference matrices, `expression matrix.xlsx` and `sample metadata.xlsx`. |
| [10 — Interactive exploration](Script_10_R_Shiny/) | `Script-10_R_Shiny.R` prepares helper data from multiple-factor results. `apps/protein_explorer/` and `apps/volcano_explorer/` display precomputed results without refitting models. |

### Pathway folder guide

**Script 5** separates `age_vs_young_ORA`, `age_vs_young_GSEA`, `female_vs_male_GSEA` and `region_vs_region_GSEA`. Regional summaries cover B–C, B–P, B–T, C–P, C–T and P–T. The age ORA folder contains two plotting-script/workbook pairs: the space-named script reads `ORA summary.xlsx`, while the underscore-named script reads `ORA_summary.xlsx`. Use the workbook referenced by the selected script.

**Script 7** uses the current organization:

```text
Script_7_mulltiple_factor_pathway/
├── age vs young in sex and region/
│   ├── S7_multiple_factor_pathway.R
│   ├── summary_table3.xlsx
│   ├── ORA summary aged vs young in sex.xlsx
│   └── ORA summary aged vs young in region.xlsx
├── female vs male/
│   ├── upset and venngraph.R
│   └── Female vs male DEP in 8 pairwise groups.xlsx
└── region vs region/
    ├── central vs bottom female young/
    ├── central vs bottom male young/
    ├── central vs peirpheral female young/
    └── central vs peripheral male young/
```

Original directory spellings are retained to preserve paths. Each regional ORA folder contains `ORA plot.R`, selected results and a pathway-summary/protein-mapping workbook.

**Script 8** separates interaction estimates from enrichment plotting. Its current GSEA script plots age-by-sex GO/Reactome, B–P GO, C–P GO and C–T GO/Reactome summaries. It does not plot a C–P Reactome result set. These interaction results are distinct from Script 5's regional main-effect comparisons.

## Running the analysis

1. Work in a copy to preserve archived outputs: scripts commonly write into their own folder.
2. Open the desired script in RStudio and inspect input filenames, worksheet names and package imports. Most scripts use `rstudioapi::getSourceEditorContext()`; adapt the working-directory setup for command-line execution.
3. Use archived processed snapshots for downstream inspection, or rerun preprocessing and explicitly refresh downstream copies.
4. Complete external annotation/enrichment and prepare the expected summary workbooks before running pathway plots.
5. Preserve sample order, protein-ID alignment, contrast direction and run provenance.

### Main data handoffs

| From | Input required by the next stage |
|---|---|
| Script 1 → Scripts 2 and 3 | `Protein_counts_summary.xlsx`, including `Raw_table_annotation`. Script 3 reads intensities directly from this worksheet; it does not require `expression_matrix.rds`. |
| Script 3 → Scripts 4, 6 and 8 | `UW.Madison.rat.lung.proteomics_data tables_<date>.RData` and `UW.Madison.rat.lung.proteomics_DIANN_protein intensity_<date>.xlsx`. |
| Script 3 → Scripts 6 and 8 | `Real_values_number_in_sample_group.xlsx` for eligibility based on observed, pre-imputation values. |
| Scripts 4 and 6 → Scripts 5 and 7 | DE results followed by external enrichment and/or preparation of comparison and pathway summaries. This handoff is not automated. |
| Script 8 DE → Script 8 GSEA plots | Interaction ranking results followed by external GSEA and curation of `GSEA/gsea summary.xlsx`. |
| Script 6 → Script 10 | `Script_6_multiple_factor_DE.RData`; after helper preparation, refresh matching copies inside each app folder. |

Supplied preprocessing snapshots are dated **2026-09-15**. Script 3 writes names using `Sys.Date()`, whereas downstream scripts read fixed snapshot names. Update references together when rerunning. Metadata filenames differ by module; `Sample_index.xlsx`, `Sample_metadata.xlsx` and `sample metadata.xlsx` are not interchangeable simply by name.

### Dependencies and Shiny applications

Use each script's package imports as its dependency list. Analyses use packages including `readxl`, `rio`, `writexl`, `openxlsx`, `dplyr`, `tidyr`, `ggplot2`, `lme4`, `lmerTest`, `emmeans`, `NormalyzerDE`, `missForest`, `imputeLCMD`, `Matrix` and `Seurat`, plus module-specific plotting packages. Some packages require Bioconductor installation.

The [Shiny environment record](Script_10_R_Shiny/sessionInfo.txt) documents that component, not a locked environment for all analyses. No project-wide `renv.lock` is supplied.

From the repository root in R:

```r
install.packages(c("shiny", "ggplot2", "DT", "dplyr", "tidyr", "ggtext"))
shiny::runApp("Script_10_R_Shiny/apps/protein_explorer")
# Stop the first app before launching the second.
shiny::runApp("Script_10_R_Shiny/apps/volcano_explorer")
```

Bundled helpers allow viewing archived results without rerunning model fitting. See the [Shiny README](Script_10_R_Shiny/README.md) for comparisons, precision and downloads. Author-provided hosted links: [protein explorer](https://feiwang-outobi.shinyapps.io/Protein-combined-analysis3/) and [volcano explorer](https://feiwang-outobi.shinyapps.io/DEP-volcano-plot3/).

## Interpreting enrichment tables

ORA uses selected protein lists; GSEA uses ranked inputs. Plotting modules consume prepared enrichment tables rather than rerun enrichment services. Preserve ranking direction, annotation version and run settings when regenerating results.

Pathway summaries contain pathway-level statistics. Protein-mapping tables identify associated input proteins; GSEA leading-edge tables identify the contributing subset. If leading-edge members come from a later rerun but statistics are retained from an earlier analysis, retain that provenance rather than imply an exact reconstruction of the original run. Interpret NES direction against the ranking contrast, particularly for interactions.

## Reproducibility notes

Checked on **2026-09-30**: all **23 R files** passed syntax-only parsing. This is not an end-to-end execution test; the analyses were not rerun for this documentation update.

- **Script 1:** the raw matrix, UniProt workspace and returned annotation are now included. Its `import()` calls still require `rio` to be loaded or namespace-qualified.
- **Script 3:** the PCA section retains a machine-specific `T:/.../ggbiplot.r` source path. Imputation has no explicit random seed; reruns may differ from saved values.
- **Script 4:** references `Sample_metadata.xlsx` although its directory supplies `Sample_index.xlsx`; later plotting sections also reference `DE summary1.xlsx`. Resolve the required filenames and worksheet layouts before rerunning.
- **Script 7:** rerunning Script 6 does not automatically refresh supplied comparison/intersection and enrichment workbooks.
- **Script 9:** both scripts now parse, but model fits, cell-type assignments and figure exports still require execution-level validation.
- **Portability:** existing names containing `:`, `<` or `>` need adjustment for Windows.

The older [syntax report](docs/SYNTAX_CHECK.md), [file-I/O index](docs/FILE_IO_INDEX.md) and [missing-input report](docs/MISSING_INPUTS.md) are historical snapshots and may contain resolved issues. Current scripts and this updated inventory take precedence.

## Archive and reuse

Supporting inputs, processed snapshots, pathway workbooks and figures remain alongside their scripts. The [added-file manifest](docs/added_files_manifest.csv) records a historical import, not every subsequent revision.

Two large Script 9 saved workspaces (~760 MB each) were never added to Git and have since been deleted, because the revised Script 9 step 1 / step 2 scripts regenerate all results from the GEO data and the proteomics matrix; the [large-file manifest](docs/large_files_not_in_git.csv) keeps their sizes and checksums for reference. Local histories, temporary Excel files and deployment records are also excluded.

This checkout does not include a license file or citation metadata. Formal release, archival deposition and reuse terms should be finalized with the manuscript; repository availability alone does not imply a license grant.
