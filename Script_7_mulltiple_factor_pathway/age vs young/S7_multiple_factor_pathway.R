# Load required packages.
library(ggplot2)
library(ggrepel)
library(readxl)
library(UpSetR)
library(dplyr)


# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Import protein-level comparison results ---------------------------------------
data0 <- read_excel("summary_table3.xlsx")

data0 <- as.data.frame(data0)

rownames(data0) <- data0$`rownames(summary_table2)`
data0 <- data0[c(2:9)]

# Select proteins upregulated in aged tissues.
Bottom_Female <- rownames(data0[data0$`Bottom Female` == "1", ])
Central_Female <- rownames(data0[data0$`Central Female` == "1", ])
Peripheral_Female <- rownames(data0[data0$`Peripheral Female` == "1", ])
Top_Female <- rownames(data0[data0$`Top Female` == "1", ])

Bottom_Male <- rownames(data0[data0$`Bottom Male` == "1", ])
Central_Male <- rownames(data0[data0$`Central Male` == "1", ])
Peripheral_Male <- rownames(data0[data0$`Peripheral Male` == "1", ])
Top_Male <- rownames(data0[data0$`Top Male` == "1", ])

# Collect protein sets for the eight sex-by-region comparisons.
list_of_proteins <- list(
  Bottom_Female = Bottom_Female,
  Central_Female = Central_Female,
  Peripheral_Female = Peripheral_Female,
  Top_Female = Top_Female,

  Bottom_Male = Bottom_Male,
  Central_Male = Central_Male,
  Peripheral_Male = Peripheral_Male,
  Top_Male = Top_Male
)

# Convert the aged-enriched protein lists to a binary membership data frame.
upset_object <- fromList(list_of_proteins)

# Make Figure S11a.
upset_plot_aged <- upset(
  upset_object,
  nsets = 8,
  order.by = "freq"
)

# Save the aged-enriched UpSet plot without an extra blank page.
pdf("Figure_S11a_UpSet_up_in_age.pdf", width = 10, height = 7)
print(upset_plot_aged, newpage = FALSE)

dev.off()

# Select proteins upregulated in young tissues.
Bottom_Female <- rownames(data0[data0$`Bottom Female` == "-1", ])
Central_Female <- rownames(data0[data0$`Central Female` == "-1", ])
Peripheral_Female <- rownames(data0[data0$`Peripheral Female` == "-1", ])
Top_Female <- rownames(data0[data0$`Top Female` == "-1", ])

Bottom_Male <- rownames(data0[data0$`Bottom Male` == "-1", ])
Central_Male <- rownames(data0[data0$`Central Male` == "-1", ])
Peripheral_Male <- rownames(data0[data0$`Peripheral Male` == "-1", ])
Top_Male <- rownames(data0[data0$`Top Male` == "-1", ])

# Collect protein sets for the eight sex-by-region comparisons.
list_of_proteins <- list(
  Bottom_Female = Bottom_Female,
  Central_Female = Central_Female,
  Peripheral_Female = Peripheral_Female,
  Top_Female = Top_Female,

  Bottom_Male = Bottom_Male,
  Central_Male = Central_Male,
  Peripheral_Male = Peripheral_Male,
  Top_Male = Top_Male
)

# Convert the young-enriched protein lists to a binary membership data frame.
upset_object <- fromList(list_of_proteins)

# Make Figure S11b.
upset_plot_young <- upset(
  upset_object,
  nsets = 8,
  order.by = "freq"
)

# Save the young-enriched UpSet plot without an extra blank page.
pdf("Figure_S11b_UpSet_up_in_young.pdf", width = 10, height = 7)
print(upset_plot_young, newpage = FALSE)

dev.off()


# Export aged-enriched protein sets by sex and region ----------------------------
Bottom_Female <- rownames(data0[data0$`Bottom Female` == "1", ])
Central_Female <- rownames(data0[data0$`Central Female` == "1", ])
Peripheral_Female <- rownames(data0[data0$`Peripheral Female` == "1", ])
Top_Female <- rownames(data0[data0$`Top Female` == "1", ])

Bottom_Male <- rownames(data0[data0$`Bottom Male` == "1", ])
Central_Male <- rownames(data0[data0$`Central Male` == "1", ])
Peripheral_Male <- rownames(data0[data0$`Peripheral Male` == "1", ])
Top_Male <- rownames(data0[data0$`Top Male` == "1", ])

Female_up <- union(Bottom_Female, Central_Female)
Female_up <- union(Female_up, Peripheral_Female)
Female_up <- union(Female_up, Top_Female) # Previously recorded count: 211.
test1 <- as.data.frame(Female_up)

Male_up <- union(Bottom_Male, Central_Male)
Male_up <- union(Male_up, Peripheral_Male)
Male_up <- union(Male_up, Top_Male) # Previously recorded count: 157.
test2 <- as.data.frame(Male_up)

Overlap_sex_up <- intersect(Female_up, Male_up) # 90 shared proteins.
export(test1, "female_up.xlsx") # 1. Aged > young in females, pooled across regions.
export(test2, "male_up.xlsx") # 2. Aged > young in males, pooled across regions.
# make venn graph in S11a manually

Bottom_up <- union(Bottom_Female, Bottom_Male) # Previously recorded count: 103.
Central_up <- union(Central_Female, Central_Male) # Previously recorded count: 187.
Peripheral_up <- union(Peripheral_Female, Peripheral_Male) # Previously recorded count: 98.
Top_up <- union(Top_Female, Top_Male) # Previously recorded count: 94.
test1 <- as.data.frame(Bottom_up)
test2 <- as.data.frame(Central_up)
test3 <- as.data.frame(Peripheral_up)
test4 <- as.data.frame(Top_up)

export(test1, "bottom_up.xlsx") # 3. Aged > young in bottom, pooled across sexes.
export(test2, "central_up.xlsx") # 4. Aged > young in central, pooled across sexes.
export(test3, "peripheral_up.xlsx") # 5. Aged > young in peripheral, pooled across sexes.
export(test4, "top_up.xlsx") # 6. Aged > young in top, pooled across sexes.
# make venn graph in S11b manually

# Export young-enriched protein sets by sex and region ---------------------------
Bottom_Female <- rownames(data0[data0$`Bottom Female` == "-1", ])
Central_Female <- rownames(data0[data0$`Central Female` == "-1", ])
Peripheral_Female <- rownames(data0[data0$`Peripheral Female` == "-1", ])
Top_Female <- rownames(data0[data0$`Top Female` == "-1", ])

Bottom_Male <- rownames(data0[data0$`Bottom Male` == "-1", ])
Central_Male <- rownames(data0[data0$`Central Male` == "-1", ])
Peripheral_Male <- rownames(data0[data0$`Peripheral Male` == "-1", ])
Top_Male <- rownames(data0[data0$`Top Male` == "-1", ])

Female_down <- union(Bottom_Female, Central_Female)
Female_down <- union(Female_down, Peripheral_Female)
Female_down <- union(Female_down, Top_Female) # Previously recorded count: 334.
test3 <- as.data.frame(Female_down)

Male_down <- union(Bottom_Male, Central_Male)
Male_down <- union(Male_down, Peripheral_Male)
Male_down <- union(Male_down, Top_Male) # Previously recorded count: 178.
test4 <- as.data.frame(Male_down)

Overlap_sex_down <- intersect(Female_down, Male_down) # 119 shared proteins.
export(test3, "female_down.xlsx") # 7. Young > aged in females, pooled across regions.
export(test4, "male_down.xlsx") # 8. Young > aged in males, pooled across regions.
# make venn graph in S11a manually

Bottom_down <- union(Bottom_Female, Bottom_Male) # Previously recorded count: 84.
Central_down <- union(Central_Female, Central_Male) # Previously recorded count: 249.
Peripheral_down <- union(Peripheral_Female, Peripheral_Male) # Previously recorded count: 86.
Top_down <- union(Top_Female, Top_Male) # Previously recorded count: 206.
test1 <- as.data.frame(Bottom_down)
test2 <- as.data.frame(Central_down)
test3 <- as.data.frame(Peripheral_down)
test4 <- as.data.frame(Top_down)

export(test1, "bottom_down.xlsx") # 9. Young > aged in bottom, pooled across sexes.
export(test2, "central_down.xlsx") # 10. Young > aged in central, pooled across sexes.
export(test3, "peripheral_down.xlsx") # 11. Young > aged in peripheral, pooled across sexes.
export(test4, "top_down.xlsx") # 12. Young > aged in top, pooled across sexes.
# make venn graph in S11b manually

# Plot representative ORA results by sex ----------------------------------------
# Read ORA summary tables curated manually in Excel.
# Import aged-enriched ORA results.
BP_female_up <- read_excel(
  "ORA summary on sex.xlsx",
  sheet = "female aged"
)
BP_female_up$Description <- sub("\\(.*", "", BP_female_up$Description)

BP_male_up <- read_excel(
  "ORA summary on sex.xlsx",
  sheet = "male aged"
)
BP_male_up$Description <- sub("\\(.*", "", BP_male_up$Description)

# Order terms by increasing absolute Log10(FDR); larger values appear at the top.
BP_female_up <- BP_female_up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure 5d
p1 <- ggplot(BP_female_up, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in aged groups only in female") +
  coord_flip() +
  ylim(0, 8) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p1)

BP_male_up <- BP_male_up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure 5f
p2 <- ggplot(BP_male_up, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#485AA7") +
  ggtitle("ORA, up in aged groups only in male") +
  coord_flip() +
  ylim(0, 14) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p2)

# Import young-enriched ORA results.
BP_female_down <- read_excel(
  "ORA summary on sex.xlsx",
  sheet = "female young"
)
BP_female_down$Description <- sub("\\(.*", "", BP_female_down$Description)

BP_male_down <- read_excel(
  "ORA summary on sex.xlsx",
  sheet = "male young"
)
BP_male_down$Description <- sub("\\(.*", "", BP_male_down$Description)

# Order terms by increasing absolute Log10(FDR); larger values appear at the top.
BP_female_down <- BP_female_down %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure 5e
p3 <- ggplot(BP_female_down, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in young groups only in female") +
  coord_flip() +
  ylim(0, 5) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p3)

BP_male_down <- BP_male_down %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure 5g
p4 <- ggplot(BP_male_down, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#485AA7") +
  ggtitle("ORA, up in young groups only in male") +
  coord_flip() +
  ylim(0, 8) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p4)

# Save all four ORA plots in a single PDF.
pdf("Age_vs_young_in_female_or_male_ORA_summary.pdf", width = 10, height = 7)

print(p1)
print(p2)
print(p3)
print(p4)

dev.off()


# Plot representative ORA results by region -------------------------------------
# Read ORA summary tables curated manually in Excel.
# Import aged-enriched ORA results for the four lung regions.

BP_bottom_up <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "bottom up in aged"
)
BP_bottom_up$Description <- sub("\\(.*", "", BP_bottom_up$Description)

BP_central_up <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "central up in aged"
)
BP_central_up$Description <- sub("\\(.*", "", BP_central_up$Description)

BP_peripheral_up <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "peri up in aged"
)
BP_peripheral_up$Description <- sub("\\(.*", "", BP_peripheral_up$Description)

BP_top_up <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "top up in aged"
)
BP_top_up$Description <- sub("\\(.*", "", BP_top_up$Description)

# Order terms by increasing absolute Log10(FDR); larger values appear at the top.
BP_bottom_up <- BP_bottom_up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S18a
p11 <- ggplot(BP_bottom_up, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in aged groups only in bottom") +
  coord_flip() +
  ylim(0, 8) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p11)

BP_central_up <- BP_central_up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S18b
p22 <- ggplot(BP_central_up, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in aged groups only in central") +
  coord_flip() +
  ylim(0, 8) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p22)

# Order terms by increasing absolute Log10(FDR); larger values appear at the top.
BP_peripheral_up <- BP_peripheral_up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S18c
p33 <- ggplot(BP_peripheral_up, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in aged groups only in peripheral") +
  coord_flip() +
  ylim(0, 8) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p33)

BP_top_up <- BP_top_up %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S18d
p44 <- ggplot(BP_top_up, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in aged groups only in top") +
  coord_flip() +
  ylim(0, 8) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p44)

# Save Figures S18a-d in a single PDF, one plot per page.
pdf("Age_>_young_in_regions_ORA_summary.pdf", width = 10, height = 7)

print(p11)
print(p22)
print(p33)
print(p44)

dev.off()

# Import young-enriched ORA results for the four lung regions.
# Read ORA summary tables curated manually in Excel.
BP_bottom_down <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "bottom up in young"
)
BP_bottom_down$Description <- sub("\\(.*", "", BP_bottom_down$Description)

BP_central_down <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "central up in young"
)
BP_central_down$Description <- sub("\\(.*", "", BP_central_down$Description)

BP_peripheral_down <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "peri up in young"
)
BP_peripheral_down$Description <- sub("\\(.*", "", BP_peripheral_down$Description)

BP_top_down <- read_excel(
  "ORA summary on region.xlsx",
  sheet = "top up in young"
)
BP_top_down$Description <- sub("\\(.*", "", BP_top_down$Description)

# Order terms by increasing absolute Log10(FDR); larger values appear at the top.
BP_bottom_down <- BP_bottom_down %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S19a
p55 <- ggplot(BP_bottom_down, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in young groups only in bottom") +
  coord_flip() +
  ylim(0, 10) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p55)

BP_central_down <- BP_central_down %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S19b
p66 <- ggplot(BP_central_down, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in young groups only in central") +
  coord_flip() +
  ylim(0, 10) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p66)

# Order terms by increasing absolute Log10(FDR); larger values appear at the top.
BP_peripheral_down <- BP_peripheral_down %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S19c
p77 <- ggplot(BP_peripheral_down, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in young groups only in peripheral") +
  coord_flip() +
  ylim(0, 10) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p77)

BP_top_down <- BP_top_down %>%
  arrange((abs(`Log10(FDR)`))) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Make Figure S19d
p88 <- ggplot(BP_top_down, aes(y = `Log10(FDR)`, x = Description)) +
  geom_bar(stat = "identity", fill = "#CD2034") +
  ggtitle("ORA, up in young groups only in top") +
  coord_flip() +
  ylim(0, 10) + # Set axis limits.
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p88)

# Save Figures S19a-d in a single PDF, one plot per page.
pdf("Age_<_young_in_regions_ORA_summary.pdf", width = 10, height = 7)

print(p55)
print(p66)
print(p77)
print(p88)

dev.off()
