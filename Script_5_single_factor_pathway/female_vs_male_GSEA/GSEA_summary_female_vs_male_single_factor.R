# Load required packages.
library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)


# Visualize precomputed GSEA results for the sex comparison.
# Input: GSEA_summary_female_vs_male_single_factor.xlsx.
# Output: GSEA_female_vs_male_single_factor.pdf.
# This script plots enrichment results; it does not perform enrichment analysis.


# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Reactome pathways enriched in the female group -------------------------
# Read the results from the first worksheet.
down <- read_excel("GSEA_summary_female_vs_male_single_factor.xlsx")


# Sort by increasing absolute log-transformed FDR value.
# Larger absolute values appear higher on the y-axis.
down <- down %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(
    Description = factor(Description, levels = unique(Description))
  )


# Map enrichment score, reported log-transformed FDR, and leading-edge size.
# Preserve the original NES direction for this comparison.
p4 <- ggplot(
  down,
  aes(
    x = NES,
    y = Description,
    color = `Log10(FDR)`,
    size = `Leading Edge Number`
  )
) +
  geom_point() +
  scale_colour_gradient(
    low = "blue",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Female vs Male, enriched in the female group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# Display the Reactome enrichment plot (Figure S7b).
print(p4)


# Export the plot to PDF.
pdf("GSEA_female_vs_male_single_factor.pdf", width = 11, height = 8.5)
print(p4)
dev.off()
