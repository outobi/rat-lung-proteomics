# Shiny applications for the rat lung proteomics study

This folder contains two Shiny applications for exploring precomputed protein differential-expression results across age, sex and lung region, together with the script used to prepare their input data.

## Archive version and data source

- Proposed manuscript release: **v1.0.0** (not yet tagged or published).
- Local data snapshot: **2026-09-21**.
- Authoritative helper data: the four files in the root of `Script_10_R_Shiny`, regenerated using the current analysis workflow.
- The helper files in `apps/` were copied from this current working version. 
- The local applications have initialized successfully with these files. 

## Folder contents

```text
Script_10_R_Shiny/
├── README.md
├── sessionInfo.txt
├── Script-10_R_Shiny.R
├── Script_6_multiple_factor_DE.RData
├── helper age.RData
├── helper sex.RData
├── helper region.RData
├── helper volcano plot.RData
└── apps/
    ├── protein_explorer/
    │   ├── app.R
    │   ├── helper age.RData
    │   ├── helper sex.RData
    │   └── helper region.RData
    └── volcano_explorer/
        ├── app.R
        └── helper volcano plot.RData
```

The root helper files are the source copies; the copies inside each application make it independently runnable. If helper data are regenerated, refresh the corresponding SAME NAME files inside `apps/` before running or deploying the applications. 

## Run locally

Use R 4.4.3 or a compatible installation. In R or RStudio, set the working directory to this `Script_10_R_Shiny` folder. Install the application dependencies if needed:

```r
install.packages(c("shiny", "ggplot2", "DT", "dplyr", "tidyr", "ggtext"))
```

This installs currently available packages, not necessarily the recorded versions below. An `renv.lock` file has not yet been created.

Launch one application at a time:

```r
# Explore selected proteins across age, sex or regional comparisons.
shiny::runApp("apps/protein_explorer")

# Explore a selected comparison using an interactive volcano plot.
shiny::runApp("apps/volcano_explorer")
```

Stop the running application before launching the other. The applications load their bundled helper files automatically; running the preparation script is not required for viewing the archived results.

## Comparisons and direction

The archive contains 40 comparisons: 8 age comparisons, 8 sex comparisons and 24 regional comparisons.

| Family | Comparisons | Positive estimate | Negative estimate |
|---|---|---|---|
| Age | Aged versus young within each sex and region | Higher in aged | Higher in young |
| Sex | Female versus male within each age and region | Higher in female | Higher in male |
| Region | Six regional pairs within each age-by-sex group | Higher in the first region of the contrast | Higher in the second region |

The six regional pairs are Bottom–Central, Bottom–Peripheral, Bottom–Top, Central–Peripheral, Central–Top and Peripheral–Top. Estimates represent model contrasts on the log2 abundance scale. For example, a positive Bottom versus Central estimate indicates higher abundance in Bottom. Use the saved contrast and group metadata to identify a comparison; list positions differ between helper files.

## Precomputed analysis and helper files

Mixed-model fitting, contrast estimation, observed-value filtering and multiple-testing adjustment are performed upstream in the multiple-factor analysis. The Shiny applications do not refit these models or recalculate adjusted P values when users change selections or thresholds.

The current workflow retains proteins with at least five observed values in either group of a comparison, then applies Benjamini–Hochberg adjustment within that comparison. `Script-10_R_Shiny.R` loads `Script_6_multiple_factor_DE.RData` and packages the results for the apps.

| File | Contents |
|---|---|
| `helper age.RData` | `df1`, `df4` and `df6` for the eight age comparisons |
| `helper sex.RData` | `df1`, `df4` and `df6` for the eight sex comparisons |
| `helper region.RData` | `df1`, `df4` and `df6` for the 24 regional comparisons |
| `helper volcano plot.RData` | DEP result lists: `age.vs.young`, `female.vs.male`, `region.vs.region` |

For the protein explorer, `df1` contains annotations and comparison statistics, `df4` is a transposed display matrix of estimates and adjusted P values, and `df6` contains protein annotations. The preparation script rounds `df4` to four significant digits. Selected-protein downloads are derived from these display matrices and therefore should not be treated as full-precision model output.

In the protein helper tables, adjusted P values are populated from the DE results; missing adjusted P values indicate that a value was not supplied by that filtered result table. They should not be interpreted as zero. The volcano helper retains the eligible DE comparison tables.

## Using the applications

**Protein explorer:** select an analysis family and one or more gene symbols to view comparison estimates, adjusted P values and plots. “Download Selected Protein Data” exports the selected results as CSV.

**Volcano explorer:** select a comparison, adjust highlighting thresholds and inspect the results table. “Download Data” exports the selected comparison; “Download Plot” exports a PNG.

The volcano application retains its original sliders and 0.1 steps. The configured default values are `log2(1.5)` for the absolute estimate cutoff and `-log10(0.05)` for the −log10(FDR) cutoff, corresponding to the manuscript criteria. Slider display or movement may round values to the slider precision; use full-precision downloaded statistics for exact manuscript cutoffs. Changing a cutoff changes highlighting, not the precomputed statistical results. 

## Relationship to the manuscript

These applications support exploration of the multiple-factor differential-expression results underlying the study's age-, sex- and region-related analyses. They are companion interfaces, not substitutes for the archived analysis scripts, full result tables or figure-generation workflow. Their selected views and interactive thresholds may differ from the panels selected for the manuscript.

Public application links cited in the manuscript:

- Volcano explorer: https://feiwang-outobi.shinyapps.io/DEP-volcano-plot3/
- Protein explorer: https://feiwang-outobi.shinyapps.io/Protein-combined-analysis3/


## Software versions

See `sessionInfo.txt` for platform and namespace details.

| Application dependency | Version |
|---|---|
| shiny | 1.12.1 |
| ggplot2 | 4.0.1 |
| DT | 0.34.0 |
| dplyr | 1.1.4 |
| tidyr | 1.3.2 |
| ggtext | 0.1.2 |

Additional packages loaded by the preparation script: lme4 1.1.38; lmerTest 3.1.3; emmeans 2.0.1; rio 1.2.4; readxl 1.4.5; gridExtra 2.3; limma 3.62.2; reshape2 1.4.5; Rmisc 1.5.1; stringr 1.6.0; multcomp 1.4.29.

## License

A license has not yet been selected. Before public release, the authors should choose a license for the code and specify the reuse terms for the processed data. No license grant is implied by this README. Add the approved license text as `LICENSE` and update this section before publication.

## Release and citation

The intended manuscript release is `v1.0.0`. No GitHub release, Zenodo deposition or DOI has been created by this documentation step. Before publication, freeze the reviewed code and helper data, record dependencies, add the approved license, create the release and archive that exact release on Zenodo. Replace the placeholders below with the issued metadata:

> [Authors]. [Year]. [Repository title]. Version 1.0.0. Zenodo. [DOI].

Use the DOI of the archived manuscript version in the code-availability statement. Record any later changes as a new version.
