# Rat lung proteomics analysis

Code and processed data accompanying the rat lung proteomics manuscript.

This repository is being assembled for manuscript revision. The first module included is the Shiny application workflow; Scripts 1–9 will be added as their review is completed. This is not yet the complete computational workflow or a finalized manuscript release.

## Available module

[Script 10: Shiny applications](Script_10_R_Shiny/README.md) contains the data-preparation script, its saved input, and two applications with their processed helper data.

- **Protein explorer:** selected-protein comparisons across age, sex and lung region.
- **Volcano explorer:** interactive differential-expression results and downloads.

The helper data are the current working version chosen by the authors, not the older helper data from the original application folders.

## Run locally

From this repository's root in R:

```r
install.packages(c("shiny", "ggplot2", "DT", "dplyr", "tidyr", "ggtext"))
shiny::runApp("Script_10_R_Shiny/apps/protein_explorer")
# Stop the first app before launching the second.
shiny::runApp("Script_10_R_Shiny/apps/volcano_explorer")
```

See the module README and `sessionInfo.txt` for recorded local package versions, comparison directions, and data preparation. Package installation above uses currently available versions; a dependency lockfile has not yet been created.

## Hosted applications

- https://feiwang-outobi.shinyapps.io/DEP-volcano-plot3/
- https://feiwang-outobi.shinyapps.io/Protein-combined-analysis3/

The authors reported completing deployment of the current applications. Hosted behavior has not been independently revalidated as part of repository preparation.

## Release, citation and license

The planned manuscript release is `v1.0.0`, after Scripts 1–10 and their supporting data have been reviewed. No release or Zenodo DOI has been issued yet. Citation metadata and the code/data license will be added before public release; no license grant is implied at this stage.
