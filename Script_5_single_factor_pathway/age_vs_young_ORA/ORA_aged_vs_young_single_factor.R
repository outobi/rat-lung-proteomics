setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)

## ORA example ----------------------------------------------------------------
# aged vs young 



# up in aged GO CC----------------------

up <- read_excel("ORA_summary.xlsx", sheet = "aged GO CC")

# Order Description based on the absolute value of Log10(FDR) in descending order
# This places the most significant (largest absolute Log10(FDR)) at the top
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p1 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO CC, Aged vs Young, enriched in the aged group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# make Figure 2d
plot(p1)

up <- read_excel("ORA_summary.xlsx", sheet = "aged GO summary")

up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p2 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO summary, Aged vs Young, enriched in the aged group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

# make Figure S3a
plot(p2)


## reactome -------
## up in aged -------


up <- read_excel("ORA_summary.xlsx", sheet = "young Reactome")

# Order Description based on the absolute value of Log10(FDR) in descending order
# This places the most significant (largest absolute Log10(FDR)) at the top
up <- up %>%
  dplyr::arrange((abs(`Log10(FDR)`))) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p3 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA Reactome, Aged vs Young, enriched in the young group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# make Figure 2f
plot(p3)

# up in young GO summary--------

up <- read_excel("ORA_summary.xlsx", sheet = "young GO summary")


up <- up %>%
  dplyr::group_by(GO) %>%
  dplyr::arrange(`Log10(FDR)`, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p4 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO summary, Aged vs Young, enriched in the young group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

# make Figure S3b
plot(p4)


# summary PDF
pdf("ORA_age_vs_young_single_factor.pdf", width = 11, height = 8.5)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()