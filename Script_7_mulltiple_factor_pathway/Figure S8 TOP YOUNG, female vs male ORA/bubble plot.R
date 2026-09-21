# Load required packages.
library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)


# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Plot representative ORA results for Figure S8 ---------------------------------
# Comparison: female versus male in the top region of young lungs.
# ORA input: 36 female-enriched proteins, submitted manually to PANTHER.
# Plotting tables were curated manually in Excel.

# Import GO enrichment results.

up <- read_excel("ORA list in female vs male top young.xlsx", sheet = "GO")

# Order terms within each GO category by increasing Log10(FDR).
up <- up %>%
  group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S8a
p1 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("ORA GO, Female vs Male in top young group, enriched in the female group") +
  # Display GO categories in separate panels.
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

print(p1)


# Plot Reactome enrichment results ---------------------------------------------

up <- read_excel("ORA list in female vs male top young.xlsx", sheet = "Reactome")

# Order terms by increasing absolute Log10(FDR); larger values appear at the top.
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S8b
p2 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    space = "Lab", na.value = "grey50", guide = "colourbar",
    aesthetics = "colour"
  ) +
  ggtitle("ORA Reactome, Female vs Male in top young group, enriched in the female group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p2)


# Save Figures S8a and S8b in a single PDF, one plot per page.
pdf("top young female > male, represatitive ORA and Reactome.pdf", width = 11, height = 8.5)

print(p1)
print(p2)

dev.off()
