# =============================================================================
# Script 5 (region_vs_region_GSEA): Dot plots of gene set enrichment analysis
#                                   (GSEA) results, region vs region
#                                   (single-factor model)
# =============================================================================
#
# Purpose:
#   Plot pre-computed GSEA results (GO and Reactome) for the regional
#   comparisons of the single-factor model in Script_4. This script only
#   draws figures; the enrichment analysis itself was run outside R.
#
# NES direction:
#   Comparisons follow the Script_4 contrasts (first region minus second):
#   Bottom - Central, Bottom - Peripheral, Bottom - Top, Central - Peripheral,
#   Central - Top, Peripheral - Top. Positive NES = enriched in the first
#   region; negative NES = enriched in the second region. Plots of
#   negative-NES gene sets show -NES on the x-axis, so the points lie on the
#   positive side; the values in the tables are unchanged.
#
# Inputs (in this folder):
#   - GSEA_summary_central_vs_peripheral.xlsx
#       Sheets "GO" (20 terms, NES > 0), "Reactome central" (3, NES > 0),
#       "Reactome peripheral" (1, NES < 0).
#   - GSEA_summary_bottom_vs_top.xlsx
#       Sheets "GO" (22 terms, NES < 0) and "Reactome" (6, NES < 0).
#       Sheets "GO BP", "GO CC" and "GO MF" are not used.
#   - Supplementary_table_region_vs_region_single_factor_GSEA.xlsx
#       9 sheets for the other comparisons; the sheet name gives the
#       direction (e.g. "GO bottom < central" = NES < 0, enriched in central).
#   Columns: GO (BP / CC / MF; not in Reactome sheets), Gene Set, Description,
#   Size, Leading Edge Number, ES, NES, P Value, FDR, and Log10(FDR).
#   Note: despite its name, the "Log10(FDR)" column contains -log10(FDR)
#   (Excel formula =-LOG10(FDR)): all values are positive, and a larger value
#   means a smaller FDR. The header has no leading "-" because Excel would
#   read it as a formula; the figure legends are relabelled "-Log10(FDR)"
#   with labs() below.
#
# Outputs (written to this folder):
#   - GSEA_region_vs_region_single_factor.pdf
#       14 pages, one per plot, in the order p1-p14:
#         Page  Plot  Comparison (enriched in)               Figure
#          1    p1    GO, Central vs Peripheral (central)       6d
#          2    p2    Reactome, Central vs Peripheral (central) S14e
#          3    p3    Reactome, Central vs Peripheral (peri)    S14f
#          4    p4    GO, Bottom vs Top (top)                   6e
#          5    p5    Reactome, Bottom vs Top (top)             S14d
#          6    p6    GO, Bottom vs Central (central)           S14a
#          7    p7    Reactome, Bottom vs Central (central)     S14b
#          8    p8    GO, Bottom vs Peripheral (peripheral)     S14c
#          9    p9    GO, Central vs Top (central)              S14g
#         10    p10   Reactome, Central vs Top (central)        S14h
#         11    p11   Reactome, Central vs Top (top)            S14i
#         12    p12   GO, Peripheral vs Top (peripheral)        S14j
#         13    p13   GO, Peripheral vs Top (top)               S14k
#         14    p14   Reactome, Peripheral vs Top (top)         S14l
#       Main figures: 6d, 6e; supplementary: S14a-S14l.
#   - GSEA_region_vs_region_single_factor.RData
#       Workspace with the plot objects p1-p14 and the last tables read.
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

# Main figures ----------------------------------------------------------------
# 2 main figures (6d-6e)
# Plot 1: GSEA GO, Central vs Peripheral, enriched in the central group.
up <- read_excel("GSEA_summary_central_vs_peripheral.xlsx",
  sheet = "GO"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction (positive-NES gene sets).
p1 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Central vs Peripheral, enriched in the central group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )
# Make Figure 6d
print(p1)


# Plot 2: GSEA Reactome, Central vs Peripheral, enriched in the central group.
up <- read_excel("GSEA_summary_central_vs_peripheral.xlsx",
  sheet = "Reactome central"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction (positive-NES gene sets).
p2 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Peripheral, enriched in the central group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14e
print(p2)


# Plot 3: GSEA Reactome, Central vs Peripheral, enriched in the peripheral group.
down <- read_excel("GSEA_summary_central_vs_peripheral.xlsx",
  sheet = "Reactome peripheral"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
down <- down %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p3 <- ggplot(down, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Peripheral, enriched in the peripheral group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14f
print(p3)


# Plot 4: GSEA GO, bottom vs top, enriched in the top group.
down <- read_excel("GSEA_summary_bottom_vs_top.xlsx",
  sheet = "GO"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
down <- down %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p4 <- ggplot(down, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Bottom vs Top, enriched in the top group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure 6e
print(p4)


# Plot 5: GSEA reactome, bottom vs top, enriched in the top group.
down <- read_excel("GSEA_summary_bottom_vs_top.xlsx",
  sheet = "Reactome"
)

# Shorten one long pathway name for display (matched by name, not row number).
# Original pathway name is too long: Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins
down$Description[down$Description == "Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins."] <- "ATP synthesis by chemiosmotic coupling"

# Sort terms and preserve the resulting y-axis order (bottom to top).
down <- down %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p5 <- ggplot(down, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA reactome, Bottom vs Top, enriched in the top group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14d
print(p5)


# Plot 6: GSEA GO, Bottom vs Central, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx", sheet = "GO bottom < central")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p6 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Bottom vs Central, enriched in the central group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14a
print(p6)


# Plot 7: GSEA Reactome, Bottom vs Central, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx", sheet = "Reactome bottom < central")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p7 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Bottom vs Central, enriched in the central group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14b
print(p7)


# Plot 8: GSEA GO, Bottom vs Peripheral, enriched in the peripheral group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx", sheet = "GO bottom < peri")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p8 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Bottom vs Peripheral, enriched in the peripheral group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14c
print(p8)


# Plot 9: GSEA GO, Central vs Top, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx", sheet = "GO central > top")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction (positive-NES gene sets).
p9 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Central vs Top, enriched in the central group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14g
print(p9)


# Plot 10: GSEA Reactome, Central vs Top, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx", sheet = "Reactome central > top")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction (positive-NES gene sets).
p10 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Top, enriched in the central group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14h
print(p10)


# Plot 11: GSEA Reactome, Central vs Top, enriched in the top group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx",
  sheet = "Reactome central < top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p11 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Top, enriched in the top group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14i
print(p11)


# Plot 12: GSEA GO, Peripheral vs Top, enriched in the peripheral group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx",
  sheet = "GO peri > top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction (positive-NES gene sets).
p12 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Peripheral vs Top, enriched in the peripheral group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14j
print(p12)


# Plot 13: GSEA GO, Peripheral vs Top, enriched in the top group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx",
  sheet = "GO peri < top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p13 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Peripheral vs Top, enriched in the top group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14k
print(p13)


# Plot 14: GSEA Reactome, Peripheral vs Top, enriched in the top group.
up <- read_excel("Supplementary_table_region_vs_region_single_factor_GSEA.xlsx",
  sheet = "Reactome peri < top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Reverse the NES sign for display (negative-NES gene sets), preserving the
# original values in the table.
p14 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Peripheral vs Top, enriched in the top group") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14l
print(p14)


# Export all 14 plots to a single PDF, with one plot per page.
# Wider pages provide space for pathway descriptions and plot legends.
region_plots <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14)

pdf("GSEA_region_vs_region_single_factor.pdf", width = 14, height = 8.5)
for (region_plot in region_plots) {
  print(region_plot)
}
dev.off()


# Save the current workspace (plot objects p1-p14).
save.image("GSEA_region_vs_region_single_factor.RData")
