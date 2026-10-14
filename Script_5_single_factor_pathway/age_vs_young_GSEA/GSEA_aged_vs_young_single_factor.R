# =============================================================================
# Script 5 (age_vs_young_GSEA): Dot plots of gene set enrichment analysis
#                               (GSEA) results, Aged vs Young (single-factor)
# =============================================================================
#
# Purpose:
#   Plot pre-computed GSEA results (GO and Reactome) for the Aged vs Young
#   comparison of the single-factor model in Script_4 (A = young adult,
#   Y = peripubertal). This script only draws figures; the enrichment
#   analysis itself was run outside R.
#
# NES direction:
#   Positive NES = gene set enriched in the Aged group (up in Aged - Young);
#   negative NES = enriched in the Young group. For the young-enriched plots
#   (Figures S4b, S4c) the x-axis shows -NES, so all points lie on the
#   positive side; the values in the table are unchanged.
#
# Input (in this folder):
#   - GSEA_summary_aged_vs_young_single_factor.xlsx
#       Selected significant gene sets, one sheet per plot:
#         "GO aged"        (4 GO terms, NES > 0)
#         "GO young"       (27 GO terms, NES < 0)
#         "Reactome young" (21 Reactome pathways, NES < 0)
#       Columns: GO (BP / CC / MF; not in the Reactome sheet), Gene Set,
#       Description, Size, Leading Edge Number, ES, NES, P Value, FDR, and
#       Log10(FDR).
#       Note: "Log10(FDR)" holds -log10(FDR) (positive; larger = more
#       significant). The column keeps this name because Excel reads a header
#       starting with "-" as a formula; the figure legends are labelled
#       "-Log10(FDR)" with labs() below.
#
# Output (written to this folder):
#   - GSEA_age_vs_young_single_factor.pdf
#       3 pages: Figure S4a, Figure S4b, Figure S4c (in this order).
#
# Plot layout (all pages): x = NES (or -NES), y = gene set, point color =
#   -log10(FDR), point size = Leading Edge Number. Gene sets are sorted by
#   -log10(FDR), so the most significant one is at the top.
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# =============================================================================

# Load required packages.
# readxl: read_excel(); dplyr: arrange() / group_by() / mutate();
# ggplot2: dot plots.
library(readxl)
library(dplyr)
library(ggplot2)


# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Step 1: GO gene sets enriched in the aged group (Figure S4a) -----------------
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


# Map enrichment score, -log10(FDR), and leading-edge size.
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
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  # Separate GO categories with independently sized y-axis panels.
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )
#make Figure S4a
print(p1)


# Step 2: GO gene sets enriched in the young group (Figure S4b) ----------------
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


# Map enrichment score, -log10(FDR), and leading-edge size.
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
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  # Separate GO categories with independently sized y-axis panels.
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

#make Figure S4b
print(p2)


# Step 3: Reactome pathways enriched in the young group (Figure S4c) -----------
up <- read_excel(
  "GSEA_summary_aged_vs_young_single_factor.xlsx",
  sheet = "Reactome young"
)


# Sort by increasing -log10(FDR).
# Larger values (more significant) appear higher on the y-axis.
up <- up %>%
  dplyr::arrange(abs(`Log10(FDR)`)) %>%
  dplyr::mutate(
    Description = factor(Description, levels = unique(Description))
  )


# Map enrichment score, -log10(FDR), and leading-edge size.
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
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# make Figure S4c
print(p3)


# Step 4: Save figures ---------------------------------------------------------

# Export the three plots to one PDF (pages: Figure S4a, S4b, S4c).
pdf("GSEA_age_vs_young_single_factor.pdf", width = 11, height = 8.5)
print(p1)
print(p2)
print(p3)
dev.off()


# Save the current workspace (plot objects p1-p3).
save.image("GSEA_aged_vs_young_single_factor.RData")
