# =============================================================================
# Script 5 (female_vs_male_GSEA): Dot plot of gene set enrichment analysis
#                                 (GSEA) results, Female vs Male (single-factor)
# =============================================================================
#
# Purpose:
#   Plot pre-computed GSEA results (Reactome) for the Female vs Male
#   comparison of the single-factor model in Script_4. This script only draws
#   the figure; the enrichment analysis itself was run outside R.
#
# NES direction:
#   Positive NES = gene set enriched in the Female group (up in Female - Male);
#   negative NES = enriched in the Male group. Both significant pathways have
#   positive NES (female-enriched), so the x-axis shows NES unchanged.
#
# Input (in this folder):
#   - GSEA_summary_female_vs_male_single_factor.xlsx
#       First (only) sheet "Reactome": 2 significant Reactome pathways
#       (Platelet degranulation; Response to elevated platelet cytosolic
#       Ca2+), NES 2.14-2.16, FDR 0.028-0.036.
#       Columns: Gene Set, Description, Size, Leading Edge Number, ES, NES,
#       P Value, FDR, and Log10(FDR).
#       Note: despite its name, the "Log10(FDR)" column contains -log10(FDR)
#       (Excel formula =-LOG10(FDR)): all values are positive, and a larger
#       value means a smaller FDR. The header has no leading "-" because
#       Excel would read it as a formula; the figure legend is relabelled
#       "-Log10(FDR)" with labs() below.
#
# Outputs (written to this folder):
#   - GSEA_female_vs_male_single_factor.pdf
#       1 page: Figure S7b.
#   - GSEA_female_vs_male_single_factor.RData
#       Workspace with the plot object p4 and the table read (down).
#
# Plot layout: x = NES, y = pathway, point color = -log10(FDR),
#   point size = Leading Edge Number. Pathways are sorted by -log10(FDR),
#   so the most significant one is at the top.
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# =============================================================================

# Load required packages.
# readxl: read_excel(); dplyr: arrange() / mutate(); ggplot2: dot plot.
library(readxl)
library(dplyr)
library(ggplot2)


# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Step 1: Reactome pathways enriched in the female group (Figure S7b) ----------
# Read the results from the first worksheet ("Reactome").
down <- read_excel("GSEA_summary_female_vs_male_single_factor.xlsx")


# Sort by increasing -log10(FDR).
# Larger values (more significant) appear higher on the y-axis.
down <- down %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(
    Description = factor(Description, levels = unique(Description))
  )


# Map enrichment score, -log10(FDR), and leading-edge size.
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
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# Display the Reactome enrichment plot (Figure S7b).
print(p4)


# Step 2: Save the figure ------------------------------------------------------

# Export the plot to PDF (Figure S7b).
pdf("GSEA_female_vs_male_single_factor.pdf", width = 11, height = 8.5)
print(p4)
dev.off()


# Save the current workspace (plot object p4).
save.image("GSEA_female_vs_male_single_factor.RData")

