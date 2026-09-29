# Load required packages.
library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)


# Visualize precomputed GSEA results for regional comparisons.
# Plot objects retain their identifiers; p13 is displayed before p12.
# This script plots supplied results; it does not run enrichment analysis.


# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# 4 MAIN FIGURES; 1 FIGURE, P3, MENTION IN THE TEXT
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

# Preserve the original NES direction and plot aesthetics.
p1 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Central vs Peripheral, enriched in the central group") +
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

# Preserve the original NES direction and plot aesthetics.
p2 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Peripheral, enriched in the central group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure 6e
print(p2)


# Plot 3: GSEA Reactome, Central vs Peripheral, enriched in the peripheral group.
down <- read_excel("GSEA_summary_central_vs_peripheral.xlsx",
  sheet = "Reactome peripheral"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
down <- down %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p3 <- ggplot(down, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Peripheral, enriched in the peripheral group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# mention in the main text, only one item
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

# Preserve the original NES direction and plot aesthetics.
p4 <- ggplot(down, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Bottom vs Top, enriched in the top group") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure 6f
print(p4)


# Plot 5: GSEA reactome, bottom vs top, enriched in the top group.
down <- read_excel("GSEA_summary_bottom_vs_top.xlsx",
  sheet = "Reactome"
)

# Retain the original manual display-label edit for the third input row.
# Original pathway name is too long: Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins
down$Description[3] <- "ATP synthesis by chemiosmotic coupling"

# Sort terms and preserve the resulting y-axis order (bottom to top).
down <- down %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p5 <- ggplot(down, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA reactome, Bottom vs Top, enriched in the top group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure 6g
print(p5)

# 9 SUPPLEMEANTRY FIGURES---------
# Plot 6: GSEA GO, Bottom vs Central, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO bottom < central")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p6 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Bottom vs Central, enriched in the central group") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14a
print(p6)


# Plot 7: GSEA Reactome, Bottom vs Central, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "Reactome bottom < central")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p7 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Bottom vs Central, enriched in the central group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14b
print(p7)


# Plot 8: GSEA GO, Bottom vs Peripheral, enriched in the peripheral group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO bottom < peri")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p8 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Bottom vs Peripheral, enriched in the peripheral group") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14c
print(p8)


# Plot 9: GSEA GO, Central vs Top, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO central > top")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p9 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Central vs Top, enriched in the central group") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14d
print(p9)


# Plot 10: GSEA Reactome, Central vs Top, enriched in the central group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "Reactome central > top")

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p10 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Top, enriched in the central group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14e
print(p10)


# Plot 11: GSEA Reactome, Central vs Top, enriched in the top group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx",
  sheet = "Reactome central < top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p11 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Central vs Top, enriched in the top group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14f
print(p11)


# Plot 13: GSEA GO, Peripheral vs Top, enriched in the peripheral group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx",
  sheet = "GO peri < top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p12 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Peripheral vs Top, enriched in the peripheral group") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14g
print(p12)


# Plot 12: GSEA GO, Peripheral vs Top, enriched in the top group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx",
  sheet = "GO peri > top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p13 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA GO, Peripheral vs Top, enriched in the top group") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

# Make Figure S14h
print(p13)


# Plot 14: GSEA Reactome, Peripheral vs Top, enriched in the top group.
up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx",
  sheet = "Reactome peri < top"
)

# Sort terms and preserve the resulting y-axis order (bottom to top).
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Preserve the original NES direction and plot aesthetics.
p14 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("GSEA Reactome, Peripheral vs Top, enriched in the top group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Make Figure S14i
print(p14)


# Export all 14 plots to a single PDF, with one plot per page.
# Wider pages provide space for pathway descriptions and plot legends.
region_plots <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p13, p12, p14)

pdf("GSEA_region_vs_region_single_factor.pdf", width = 14, height = 8.5)
for (region_plot in region_plots) {
  print(region_plot)
}
dev.off()
