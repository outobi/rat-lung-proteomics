setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)

## ORA example ----------------------------------------------------------------
# central > peripheral female young, n = 48

up <- read_excel("Central_vs_Bottom_Female_Young_ORA_selected.xlsx", sheet = "GO")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_grodown= TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p1 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + 
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO, Central vs Bottom in female young group, enriched in the central group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p1)


## reactome -------
## up in central -------

up <- read_excel("Central_vs_Bottom_Female_Young_ORA_selected.xlsx", sheet = "Reactome")

# Order Description based on the absolute value of Log10(FDR) in descending order
# This places the most significant (largest absolute Log10(FDR)) at the top
up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p2 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA Reactome, Central vs Bottom in male young group, enriched in the central group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p2)


pdf("Central vs Bottom in female young group summary barplot ORA.pdf", width = 11, height = 8.5)
print(p1)
print(p2)
dev.off()

