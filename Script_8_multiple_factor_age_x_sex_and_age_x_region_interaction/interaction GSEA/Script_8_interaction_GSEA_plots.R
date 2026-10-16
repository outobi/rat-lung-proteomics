# =============================================================================
# Script 8 (interaction GSEA): Dot plots of GSEA results for the Age x Sex and
#   Age x Region interaction contrasts (multiple-factor model)
# =============================================================================
#
# Purpose:
#   Plot selected gene set enrichment analysis (GSEA) results for the
#   interaction contrasts of "interaction contrast of contrast"
#   (age_sex_and_age_region_interaction.R). This script only draws figures;
#   the GSEA was run in WebGestalt (https://www.webgestalt.org/) with all
#   3,600 proteins ranked by the interaction estimate:
#     - Age x Sex:    (Aged - Young)Female - (Aged - Young)Male
#     - Age x Region: (Aged - Young)Region1 - (Aged - Young)Region2
#   (Aged = young adult; Young = peripubertal.)
#
# NES direction:
#   Positive NES = gene set enriched where the age effect is larger in the
#   first group (Female, or Region1); negative NES = larger in the second
#   group (Male, or Region2). P4 (negative NES) shows -NES on the x-axis, so
#   all plots have points on the positive side; table values are unchanged.
#
# Input (in this folder):
#   - gsea summary.xlsx
#       Selected significant gene sets, one sheet per plot:
#         "female vs male go"        (2 GO BP / CC terms, NES > 0)
#         "female vs male reactome"  (3 Reactome pathways, NES > 0)
#         "bottom vs peri GO"        (5 GO MF terms, NES > 0)
#         "central vs peri GO"       (1 GO CC term, NES < 0)
#         "central vs top GO"        (3 GO CC terms, NES > 0)
#         "central vs top reactome"  (9 Reactome pathways, NES > 0)
#       Sheet "ITEM SUMMARY" is not used.
#       Columns: GO (BP / CC / MF; not in Reactome sheets), Gene Set,
#       Description, Size, Leading Edge Number, ES, NES, P Value, FDR,
#       Log10(FDR).
#       Note: the "Log10(FDR)" column contains -log10(FDR) (positive; larger =
#       more significant); the legends are labelled "-Log10(FDR)".
#
# Outputs (written to this folder):
#   - P1 ... P6 GSEA <...>.pdf  (one PDF per plot, file names as below)
#   - Script_8_interaction_GSEA_plots.RData
#       Workspace with the plot objects p1-p6.
#
# Plot layout (all plots): x = NES (P4: -NES), y = gene set, point color = -log10(FDR),
#   point size = Leading Edge Number. Gene sets are sorted by -log10(FDR),
#   so the most significant one is at the top.
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

## GSEA example ----------------------------------------------------------------
# aged vs young 
# up in female GO BP CC, no MF, N = 2 -----------------------------

up <- read_excel("gsea summary.xlsx", sheet = "female vs male go")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p1 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +  scale_size_continuous(breaks = c(6, 7),range = c(12, 14))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group") +
  # Facet by GO category to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p1)

ggsave(
  filename = "P1 GSEA GO, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group.pdf",
  plot = p1,
  device = "pdf",
  width = 11,
  height = 8.5
)

# up in female REACTOME , N = 3-----------------------------

up <- read_excel("gsea summary.xlsx", sheet = "female vs male reactome")


up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p2 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() + scale_size_continuous(breaks = c(7, 15),range = c(7, 15))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA Reactome, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p2)

ggsave(
  filename = "P2 GSEA Reactome, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group.pdf",
  plot = p2,
  device = "pdf",
  width = 11,
  height = 8.5
)


# age * region interaction gsea plot------------

# bottom vs peripheral GO MF, N = 5

up <- read_excel("gsea summary.xlsx", sheet = "bottom vs peri GO")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p3 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() + scale_size_continuous(range = c(3, 11.5), breaks =c(5, 15, 25), limits = c(1,30))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO MF, (Aged vs Young)Bottom - (Aged vs Young)Peripheral, enriched in the bottom group") +
  # Facet by GO category to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p3)

ggsave(
  filename = "P3 GSEA GO MF, (Aged vs Young)Bottom - (Aged vs Young)Peripheral, enriched in the bottom group.pdf",
  plot = p3,
  device = "pdf",
  width = 11,
  height = 8.5
)



# central vs peripheral GO CC, N = 1

up <- read_excel("gsea summary.xlsx", sheet = "central vs peri GO")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
# Negative NES (larger age effect in peripheral): reverse the sign for display,
# preserving the original values in the table.
p4 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Peripheral, enriched in the peripheral group") +
  # Facet by GO category to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  labs(x = "-NES", colour = "-Log10(FDR)") + # Axis and legend titles
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p4)

ggsave(
  filename = "P4 GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Peripheral, enriched in the peripheral group.pdf",
  plot = p4,
  device = "pdf",
  width = 11,
  height = 8.5
)






# central vs top GO CC, N = 3

up <- read_excel("gsea summary.xlsx", sheet = "central vs top GO")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p5 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +  scale_size_continuous(range = c(2, 10), breaks =c(2, 22, 42), limits = c(2,42))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group") +
  # Facet by GO category to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p5)

ggsave(
  filename = "P5 GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group.pdf",
  plot = p5,
  device = "pdf",
  width = 11,
  height = 8.5
)


# central vs top GO Reactome, N = 9

up <- read_excel("gsea summary.xlsx", sheet = "central vs top reactome")


up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p6 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() + scale_size_continuous(range = c(6, 8), breaks =c(25, 30, 35), limits = c(25,35))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA Reactome, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group") +
  labs(colour = "-Log10(FDR)") + # Legend title (column holds -log10(FDR))
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p6)
ggsave(
  filename = "P6 GSEA Reactome, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group.pdf",
  plot = p6,
  device = "pdf",
  width = 11,
  height = 8.5
)


# Save the current workspace (plot objects p1-p6).
save.image("Script_8_interaction_GSEA_plots.RData")
