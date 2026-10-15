# =============================================================================
# Script 7 (ORA region vs region): ORA dot plots, Central vs Peripheral,
#   female young (peripubertal) rats (multiple-factor model)
# =============================================================================
#
# Purpose:
#   Plot selected over-representation analysis (ORA) results for the proteins
#   higher in the central region than in the peripheral region within female
#   young rats. This script only draws figures; the ORA was run in PANTHER
#   (Overrepresentation Test, Fisher's exact test, FDR correction; reference
#   list: all Rattus norvegicus genes).
#
# ORA input (not stored in this folder):
#   77 UniProt IDs = Script_6 region DEPs (Central - Peripheral, female young:
#   FDR < 0.05, estimate > 0 and |log2 FC| > log2(1.5)), i.e. higher in
#   Central. They were submitted to PANTHER; the complete PANTHER output is
#   not kept here, only the selected terms below.
#
# Input (in this folder):
#   - Central_vs_Peripheral_Female_Young_ORA_selected.xlsx
#       Selected significant terms; sheets used: "GO" (BP / CC / MF terms with
#       a GO column for faceting) and "Reactome". Columns: Description,
#       REFLIST, Size, Expect, over/under, Ratio (fold enrichment), P value,
#       FDR, Log10(FDR).
#       Note: the "Log10(FDR)" column contains -log10(FDR) (positive; larger =
#       more significant); the legends are labelled "-Log10(FDR)".
#
# Outputs (written to this folder):
#   - Central vs Peripheral in female young group summary ORA.pdf
#       2 pages: GO dot plot (faceted by BP / CC / MF), Reactome dot plot.
#   - ORA plot.RData
#       Workspace with the plot objects p1, p2.
#
# Plot layout: x = Ratio (fold enrichment), y = term, point color =
#   -log10(FDR), point size = Size (input proteins in the term). Terms are
#   sorted by -log10(FDR), so the most significant one is at the top.
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

## ORA example ----------------------------------------------------------------
# GO terms enriched in proteins higher in central than peripheral, female young (n = 77 input proteins).

up <- read_excel("Central_vs_Peripheral_Female_Young_ORA_selected.xlsx", sheet = "GO")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p1 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + 
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO, Central vs Peripheral in female young group, enriched in the central group") +
  # Facet by GO category (BP / CC / MF) to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p1)


## reactome -------
## up in central -------

up <- read_excel("Central_vs_Peripheral_Female_Young_ORA_selected.xlsx", sheet = "Reactome")

# Sort terms by increasing -log10(FDR); ggplot draws the last factor level at
# the top, so the most significant term appears at the top of the plot.
up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p2 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA Reactome, Central vs Peripheral in female young group, enriched in the central group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p2)


pdf("Central vs Peripheral in female young group summary ORA.pdf", width = 11, height = 8.5)
print(p1)
print(p2)
dev.off()


# Save the current workspace (plot objects p1, p2).
save.image("ORA plot.RData")

