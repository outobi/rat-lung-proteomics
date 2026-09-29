setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)

## ORA example ----------------------------------------------------------------
# aged vs young 



# up in aged GO CC----------------------

up <- read_excel("ORA summary.xlsx", sheet = "aged GO CC")

# Order Description based on the absolute value of Log10(FDR) in descending order
# This places the most significant (largest absolute Log10(FDR)) at the top
up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot Figure 2d
p1 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + scale_size_continuous(breaks = c(1, 30, 60, 90),
                                       limits = c(0, 90)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO CC, Aged vs Young, enriched in the aged group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p1)

up <- read_excel("ORA summary.xlsx", sheet = "aged GO select")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_grodown= TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot, Figure S3a
p2 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + scale_size_continuous(breaks = c(1, 30, 60, 90),
                                        limits = c(0, 90)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO, Aged vs Young, enriched in the aged group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p2)


## reactome -------
## up in aged -------


up <- read_excel("ORA summary.xlsx", sheet = "young Reactome select")

# Order Description based on the absolute value of Log10(FDR) in descending order
# This places the most significant (largest absolute Log10(FDR)) at the top
up <- up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot Figure 3f
p3 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() + scale_size_continuous(breaks = c(4,8,12,16,20),
                                       limits = c(0, 20)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA Reactome, Aged vs Young, enriched in the young group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p3)

# up in young GO SELECT--------

up <- read_excel("ORA summary.xlsx", sheet = "young GO select")


up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_grodown= TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot,Figure S3b
p4 <- ggplot(up, aes(x = Ratio, y = Description, color = `Log10(FDR)`, size = Size)) +
  geom_point() +  scale_size_continuous(breaks = c(1, 20, 40, 60, 80),
                                         limits = c(0, 80)) +
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("ORA GO, Aged vs Young, enriched in the young group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p4)


pdf("Aged vs young summary ORA single factor.pdf", width = 11, height = 8.5)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()


save.image("ORA aged vs young single factor.RData")




