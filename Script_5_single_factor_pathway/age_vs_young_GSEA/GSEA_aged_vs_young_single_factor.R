# Load required packages.
library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)


# Visualize precomputed GSEA results for the age comparison.
# Input: GSEA_summary_aged_vs_young_single_factor.xlsx.
# This section creates three plots; it does not perform enrichment analysis.


# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# GO pathways enriched in the aged group -------------------------------
up <- read_excel(
  "GSEA_summary_aged_vs_young_single_factor.xlsx",
  sheet = "GO aged"
)


# Sort within each GO category and preserve the order for the y-axis.
# Factor levels are displayed from bottom to top.
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Description = factor(Description, levels = unique(Description))
  )


# Map enrichment score, reported log-transformed FDR, and leading-edge size.
p1 <- ggplot(
  up,
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
  ggtitle("GSEA GO, Aged vs Young, enriched in the aged group") +
  # Separate GO categories with independently sized y-axis panels.
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )
#make Figure S4a
print(p1)


# GO pathways enriched in the young group -------------------------------
down <- read_excel(
  "GSEA_summary_aged_vs_young_single_factor.xlsx",
  sheet = "GO young"
)


# Sort within each GO category and preserve the order for the y-axis.
# Factor levels are displayed from bottom to top.
down <- down %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Description = factor(Description, levels = unique(Description))
  )


# Map enrichment score, reported log-transformed FDR, and leading-edge size.
# Reverse the NES sign for display, preserving the original values in the table.
p2 <- ggplot(
  down,
  aes(
    x = (-1) * NES,
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
  ggtitle("GSEA GO, Aged vs Young, enriched in the young group") +
  # Separate GO categories with independently sized y-axis panels.
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

#make Figure S4b
print(p2)


# Reactome pathways enriched in the young group -------------------------------
up <- read_excel(
  "GSEA_summary_aged_vs_young_single_factor.xlsx",
  sheet = "Reactome young"
)


# Sort by increasing absolute log-transformed FDR value.
# Larger absolute values appear higher on the y-axis.
up <- up %>%
  dplyr::arrange(abs(`Log10(FDR)`)) %>%
  dplyr::mutate(
    Description = factor(Description, levels = unique(Description))
  )


# Map enrichment score, reported log-transformed FDR, and leading-edge size.
# Reverse the NES sign for display, preserving the original values in the table.
p3 <- ggplot(
  up,
  aes(
    x = (-1) * NES,
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
  ggtitle("GSEA Reactome, Aged vs Young, enriched in the young group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# make Figure S4c
print(p3)

# summary PDF
pdf("GSEA_age_vs_young_single_factor.pdf", width = 11, height = 8.5)
print(p1)
print(p2)
print(p3)
dev.off()