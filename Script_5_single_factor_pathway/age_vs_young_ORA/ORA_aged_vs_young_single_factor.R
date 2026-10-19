# =============================================================================
# Script 5 (age_vs_young_ORA): Dot plots of over-representation analysis (ORA)
#                              results, Aged vs Young (single-factor model)
# =============================================================================
#
# Purpose:
#   Plot pre-computed ORA results (GO and Reactome) for proteins that differ
#   between Aged (A = 7-11m) and Young (Y = 7w) rats in the
#   single-factor model of Script_4. This script only draws figures; the
#   enrichment analysis itself was run outside R.
#
# Input (in this folder):
#   - ORA_summary.xlsx
#       One sheet per selection of enriched terms. Columns: GO (BP / CC / MF;
#       not in the Reactome sheet), Description, Ref size (term size in the
#       background), Size (input proteins in the term), Expect, Ratio
#       (enrichment ratio = Size / Expect), P Value, FDR, and Log10(FDR).
#       Note: "Log10(FDR)" holds -log10(FDR) (positive; larger = more
#       significant; Excel formula =-LOG10(FDR)). The column keeps this name
#       because Excel reads a header starting with "-" as a formula; the
#       figure legends are labelled "-Log10(FDR)" with labs() below.
#       Sheets used: "aged GO CC", "aged GO select", "young Reactome select",
#       "young GO select". Sheets "aged  GO BP" and "aged GO MF" are not used.
#
# Outputs (written to this folder):
#   - ORA_age_vs_young_single_factor.pdf
#       4 pages: Figure 2d, Figure S3a, Figure 3f, Figure S3b (in this order).
#   - ORA_aged_vs_young_single_factor.RData
#       Plot objects p1-p4 and the last table read (up).
#
# Plot layout (all pages): x = Ratio, y = term, point color = -log10(FDR),
#   point size = Size. Terms are sorted by -log10(FDR), so the most
#   significant term is at the top.
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# =============================================================================

# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages.
# readxl: read_excel(); dplyr: arrange() / group_by() / mutate();
# ggplot2: dot plots.
library(readxl)
library(dplyr)
library(ggplot2)


# Step 1: Enriched in the aged group -------------------------------------------

# GO cellular component terms enriched in the aged group (Figure 2d).
up <- read_excel("ORA_summary.xlsx", sheet = "aged GO CC")

# Sort terms by increasing -log10(FDR); ggplot draws the last factor level
# at the top, so the most significant term appears at the top of the plot.
up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot (Figure 2d).
# Point-size scale fixed to 0-90 (largest Size in this sheet: 86).
p1 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + scale_size_continuous(breaks = c(1, 30, 60, 90),
                                       limits = c(0, 90)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO CC, Aged vs Young, enriched in the aged group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p1)


# Selected GO terms (BP, CC, MF) enriched in the aged group (Figure S3a).
up <- read_excel("ORA_summary.xlsx", sheet = "aged GO select")

# Sort terms by increasing -log10(FDR) within each GO category.
up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot (Figure S3a), one facet per GO category.
# Point-size scale fixed to 0-90 (largest Size in this sheet: 86).
p2 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + scale_size_continuous(breaks = c(1, 30, 60, 90),
                                        limits = c(0, 90)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO, Aged vs Young, enriched in the aged group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  # Facet by GO category (BP / CC / MF) to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p2)


# Step 2: Enriched in the young group ------------------------------------------

# Selected Reactome pathways enriched in the young group (Figure 3f).
up <- read_excel("ORA_summary.xlsx", sheet = "young Reactome select")

# Sort terms by increasing -log10(FDR) (most significant at the top).
up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot (Figure 3f).
# Point-size scale fixed to 0-20 (largest Size in this sheet: 17).
p3 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + scale_size_continuous(breaks = c(4,8,12,16,20),
                                       limits = c(0, 20)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA Reactome, Aged vs Young, enriched in the young group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p3)


# Selected GO terms (BP, CC, MF) enriched in the young group (Figure S3b).
up <- read_excel("ORA_summary.xlsx", sheet = "young GO select")

# Sort terms by increasing -log10(FDR) within each GO category.
up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot (Figure S3b), one facet per GO category.
# Point-size scale fixed to 0-80 (largest Size in this sheet: 66).
p4 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +  scale_size_continuous(breaks = c(1, 20, 40, 60, 80),
                                         limits = c(0, 80)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO, Aged vs Young, enriched in the young group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  # Facet by GO category (BP / CC / MF) to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p4)


# Step 3: Save figures and workspace -------------------------------------------

# Export the four plots to one PDF (pages: Figure 2d, S3a, 3f, S3b).
pdf("ORA_age_vs_young_single_factor.pdf", width = 11, height = 8.5)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()


# Save the current workspace (plot objects p1-p4).
save.image("ORA_aged_vs_young_single_factor.RData")
