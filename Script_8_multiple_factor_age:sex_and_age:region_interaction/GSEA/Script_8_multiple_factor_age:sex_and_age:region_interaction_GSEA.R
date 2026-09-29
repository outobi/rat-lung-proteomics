setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readxl)
library(rio)
library(dplyr)
library(ggplot2)
library(ggrepel)

## GSEA example ----------------------------------------------------------------
# aged vs young 
# up in female GO BP CC, no MF, N = 2 -----------------------------

up <- read_excel("gsea summary.xlsx", sheet = "female vs male go")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_grodown= TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p1 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +  scale_size_continuous(breaks = c(6, 7),range = c(12, 14))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
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
p2 <- ggplot(up, aes(x = (-1) * NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() + scale_size_continuous(breaks = c(7, 15),range = c(7, 15))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA Reactome, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group") +
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
  arrange(`Log10(FDR)`, .by_grodown= TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p3 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() + scale_size_continuous(range = c(3, 11.5), breaks =c(5, 15, 25), limits = c(1,30))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO MF, (Aged vs Young)Bottom - (Aged vs Young)Peripheral, enriched in the bottom group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
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



# central vs peripheral GO CC, N = 3

up <- read_excel("gsea summary.xlsx", sheet = "central vs peri GO")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_grodown= TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p4 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +  scale_size_continuous(range = c(2, 10), breaks =c(2, 22, 42), limits = c(2,42))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Peripheral, enriched in the central group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), # Adjust text size for readability
        strip.text.y = element_text(angle = 0, hjust = 0)) # Keep facet labels horizontal

print(p4)

ggsave(
  filename = "P4 GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Peripheral, enriched in the central group.pdf",
  plot = p4,
  device = "pdf",
  width = 11,
  height = 8.5
)






# central vs top GO CC, N = 3

up <- read_excel("gsea summary.xlsx", sheet = "central vs top GO")

up <- up %>%
  group_by(GO) %>%
  arrange(`Log10(FDR)`, .by_grodown= TRUE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Create the plot
p5 <- ggplot(up, aes(x = NES, y = Description, color = `Log10(FDR)`, size = `Leading Edge Number`)) +
  geom_point() +  scale_size_continuous(range = c(2, 10), breaks =c(2, 22, 42), limits = c(2,42))+
  scale_colour_gradient(low = "blue", high = "red",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  ggtitle("GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group") +
  # Facet by FactorType to create separate clusters
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") +
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
  ggtitle("GSEA GO Reactome, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Adjust text size for readability

# Print the plot
print(p6)
ggsave(
  filename = "P6 GSEA GO Reactome, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group.pdf",
  plot = p6,
  device = "pdf",
  width = 11,
  height = 8.5
)





save.image("GSEA plotting summary.RData")
