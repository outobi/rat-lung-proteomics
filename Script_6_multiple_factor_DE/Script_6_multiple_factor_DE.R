# =============================================================================
# Script 6: Multiple-factor differential abundance (within-subgroup contrasts)
# =============================================================================
#
# Purpose:
#   For each of the 3,600 filtered and imputed protein groups from Script_3,
#   fit a linear mixed model with all two-way interactions and a random
#   intercept per animal:
#       Protein ~ Age + Sex + Region + Age:Region + Sex:Region + Age:Sex
#                 + (1 | PID)
#   and test each factor within subgroups of the other two (emmeans):
#     - Age:    Aged - Young     within each Region x Sex   ( 8 comparisons)
#     - Sex:    Female - Male    within each Region x Age   ( 8 comparisons)
#     - Region: 6 region pairs   within each Age x Sex      (24 comparisons)
#   (A = Aged = 7-11m; Y = Young = 7w.)
#   Then count DEPs and draw bar plots, volcano plots, a heatmap, boxplots
#   and a correlation map for the age comparisons.
#
# Method 4 filter and multiple testing:
#   For each comparison, keep proteins with >= 5 observed (non-imputed)
#   values in at least one of the two groups compared (counts from
#   Real_values_number_in_sample_group.xlsx), then BH-adjust the p-values
#   within that comparison. Unfiltered results are kept in
#   *_comparisons_unfiltered (not adjusted).
#   DEP: FDR (p.value.adjusted) < 0.05 and |estimate| > log2(1.5).
#   Positive estimate = higher in the first group of the contrast.
#
# Inputs (all in this folder):
#   - UW.Madison.rat.lung.proteomics_data_tables.RData  (output of Script_3)
#       data_table$data_imp: imputed, normalized log2 intensities
#       (3,600 proteins x S1-S96); info_table; matrisome_anno_table.
#   - UW.Madison.rat.lung.proteomics_DIANN_protein_intensity.xlsx
#       (output of Script_3) Sheet "data_imp": the 7 annotation columns are
#       attached to the results; columns 8-103 (S1-S96) are used for the
#       boxplots. Rows are in the same order as data_table$data_imp.
#   - Sample_metadata.xlsx
#       96 rows in S1-S96 order (Sample.ID, Sample.group, PID, Age, Region,
#       Sex, ...); used as the model table (PID = animal ID).
#   - Real_values_number_in_sample_group.xlsx  (output of Script_3)
#       Number of observed values (0-6) per protein in each of the 16 groups.
#   - Sample_group_order.txt
#       Sample order (Sample.group, Sample.ID) for the boxplots.
#   - 38 matrisome DEPs with unique features in 8 pairwise comparisons.txt
#       Uniprot.ID and Gene.symbol of the 38 matrisome proteins shown in the
#       heatmap and boxplots.
#
# Outputs (written to this folder):
#   - age.vs.young4.xlsx, female.vs.male4.xlsx, region.vs.region4.xlsx
#       Method 4 results, one sheet per comparison (8 / 8 / 24): annotations,
#       group labels, observed-value counts, contrast estimate, SE, df,
#       t.ratio, p.value and p.value.adjusted.
#   - Method 4 DE_protein_summary_all_comparisons_FDR0.05_FC1.5.xlsx
#       DEP counts (total and matrisome, up / down) per comparison.
#   - DE_protein_summary_plots.pdf            (Figures 3a, 5b, S15a-d)
#   - DE_matrisome_protein_summary_plots.pdf  (Figures 3b, 5c, S15e-h)
#   - age_vs_young_DEP_multiple_factor_volcano_plot.pdf
#       8 pages, one per age comparison (Figures 3d-g, S5a-d).
#   - summary_table3.xlsx
#       Direction (1 = up in Aged, -1 = up in Young, 0 = not a DEP or
#       filtered out) of the union of age DEPs (670 proteins in the reported
#       run) across the 8 age comparisons.
#   - Age_vs_young_matrisome_DEP_heatmap.pdf  (Figure 3h)
#       log2 fold changes of the 38 matrisome proteins (0 = not a DEP or
#       filtered out).
#   - All_Proteins_Boxplots.pdf               (Figures 5h, 5i, S6a-d)
#       One page per protein (38), all 6 samples per group.
#   - Age_vs_young_DEP_correlation_map.pdf    (Figure 3c)
#       Pearson correlations of the DEP direction matrix between the 8 age
#       comparisons.
#   - Script_6_multiple_factor_DE.RData
#       Full R workspace, including intermediate objects.
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# Model fitting takes several minutes (3,600 lmer fits, 3 emmeans each).
# =============================================================================

# Load required packages.
# lme4 + lmerTest: lmer(); emmeans: emmeans() / contrast();
# dplyr + tidyr: table handling (bind_rows, gather, ...); rio: export();
# readxl: read_excel(); readr: read_delim(); ggplot2 + ggrepel: plots;
# gridExtra: grid.arrange(); pheatmap + RColorBrewer: heatmap;
# reshape2: melt().
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(tidyr)
library(rio)
library(readxl)
library(readr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(reshape2)


# Set the working directory to this script's folder (RStudio only).
# data7 contains normalized, imputed log2 abundances; observed-value counts
# are from before imputation. workdata must contain one row per data7
# column, in the same sample order (S1-S96).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load the saved imputed matrices and sample information.
load("UW.Madison.rat.lung.proteomics_data_tables.RData")


# Read the annotated imputed-intensity worksheet for output annotations.
# The first seven columns are assumed to contain protein annotations.
data_imp <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein_intensity.xlsx",
  sheet = "data_imp"
)

data7 <- data_table$data_imp

workdata <- read_excel("Sample_metadata.xlsx")

sample_group_number <- read_excel("Real_values_number_in_sample_group.xlsx")
sample_group_number2 <- sample_group_number[sample_group_number$Uniprot %in% rownames(data7), ]

sample_order <- c(
  "YTM", "YCM", "YPM", "YBM", "YTF", "YCF", "YPF", "YBF",
  "ATM", "ACM", "APM", "ABM", "ATF", "ACF", "APF", "ABF"
)


# Preserve sample-factor objects used by later analysis sections.
f <- info_table$Sample.Group
Region <- info_table$Sample.Region
Age <- info_table$Sample.Age
Sex <- info_table$Sample.Sex
ID <- as.factor(info_table$Sample.PID)


# Use the prepared protein IDs and annotation columns.
protein_ids <- rownames(data7)
df0 <- as.data.frame(data_imp[, 1:7])


# Define the 8 age, 8 sex, and 24 region comparisons in their original order.
age_definitions <- data.frame(
  Name = c(
    "bottom.female",
    "central.female",
    "peri.female",
    "top.female",
    "bottom.male",
    "central.male",
    "peri.male",
    "top.male"
  ),
  Group1 = c("ABF", "ACF", "APF", "ATF", "ABM", "ACM", "APM", "ATM"),
  Group2 = c("YBF", "YCF", "YPF", "YTF", "YBM", "YCM", "YPM", "YTM")
)

sex_definitions <- data.frame(
  Name = c(
    "bottom.age",
    "central.age",
    "peri.age",
    "top.age",
    "bottom.young",
    "central.young",
    "peri.young",
    "top.young"
  ),
  Group1 = c("ABF", "ACF", "APF", "ATF", "YBF", "YCF", "YPF", "YTF"),
  Group2 = c("ABM", "ACM", "APM", "ATM", "YBM", "YCM", "YPM", "YTM")
)

region_definitions <- data.frame(
  Name = c(
    "bottom.vs.central.female.age",
    "bottom.vs.peri.female.age",
    "bottom.vs.top.female.age",
    "central.vs.peri.female.age",
    "central.vs.top.female.age",
    "peri.vs.top.female.age",
    "bottom.vs.central.male.age",
    "bottom.vs.peri.male.age",
    "bottom.vs.top.male.age",
    "central.vs.peri.male.age",
    "central.vs.top.male.age",
    "peri.vs.top.male.age",
    "bottom.vs.central.female.young",
    "bottom.vs.peri.female.young",
    "bottom.vs.top.female.young",
    "central.vs.peri.female.young",
    "central.vs.top.female.young",
    "peri.vs.top.female.young",
    "bottom.vs.central.male.young",
    "bottom.vs.peri.male.young",
    "bottom.vs.top.male.young",
    "central.vs.peri.male.young",
    "central.vs.top.male.young",
    "peri.vs.top.male.young"
  ),
  Group1 = c(
    "ABF", "ABF", "ABF", "ACF", "ACF", "APF",
    "ABM", "ABM", "ABM", "ACM", "ACM", "APM",
    "YBF", "YBF", "YBF", "YCF", "YCF", "YPF",
    "YBM", "YBM", "YBM", "YCM", "YCM", "YPM"
  ),
  Group2 = c(
    "ACF", "APF", "ATF", "APF", "ATF", "ATF",
    "ACM", "APM", "ATM", "APM", "ATM", "ATM",
    "YCF", "YPF", "YTF", "YPF", "YTF", "YTF",
    "YCM", "YPM", "YTM", "YPM", "YTM", "YTM"
  )
)


# Fit one mixed model per protein and obtain all three contrast families.
# Retain the original two-way interactions and subject-specific random intercept.
# Preserve the original emmeans calls and unadjusted contrast P values.
all_comparisons <- vector("list", nrow(data7))
all_comparisons2 <- vector("list", nrow(data7))
all_comparisons3 <- vector("list", nrow(data7))

for (i in seq_len(nrow(data7))) {
  workdata$Protein <- as.numeric(data7[i, ])
  res <- lmerTest::lmer(
    Protein ~ Age + Sex + Region + Age:Region + Sex:Region + Age:Sex + (1 | PID),
    data = workdata
  )

  age_emmeans <- emmeans(res, trt.vs.ctrl ~ Age | Region + Sex, adjust = "none")
  sex_emmeans <- emmeans(res, trt.vs.ctrl ~ Sex | Region + Age, adjust = "none")
  region_emmeans <- emmeans(res, trt.vs.ctrl ~ Region | Sex + Age, adjust = "none")

  all_comparisons[[i]] <- as.data.frame(summary(contrast(age_emmeans, "pairwise")))
  all_comparisons2[[i]] <- as.data.frame(summary(contrast(sex_emmeans, "pairwise")))
  all_comparisons3[[i]] <- as.data.frame(summary(contrast(region_emmeans, "pairwise")))
}

all_comparisons_df <- bind_rows(all_comparisons)
all_comparisons_df2 <- bind_rows(all_comparisons2)
all_comparisons_df3 <- bind_rows(all_comparisons3)


# Build annotated comparison tables and verify contrast direction and group labels.
# Keep unfiltered raw results for the later Shiny preparation; do not adjust them.
build_comparisons <- function(results, definitions, family) {
  regions <- c(B = "Bottom", C = "Central", P = "Peripheral", T = "Top")
  ages <- c(A = "Aged", Y = "Young")
  sexes <- c(F = "Female", M = "Male")
  n_comparisons <- nrow(definitions)
  stopifnot(nrow(results) == length(protein_ids) * n_comparisons)

  comparisons <- lapply(seq_len(n_comparisons), function(j) {
    group1 <- definitions$Group1[j]
    group2 <- definitions$Group2[j]
    comparison <- results[seq(j, nrow(results), by = n_comparisons), , drop = FALSE]
    rownames(comparison) <- protein_ids

    expected_contrast <- switch(family,
      age = "Aged - Young",
      sex = "Female - Male",
      region = paste(regions[substr(group1, 2, 2)], regions[substr(group2, 2, 2)], sep = " - ")
    )
    stopifnot(all(as.character(comparison$contrast) == expected_contrast))
    if (family != "region") {
      stopifnot(all(as.character(comparison$Region) == regions[substr(group1, 2, 2)]))
    }
    if (family != "age") {
      stopifnot(all(as.character(comparison$Age) == ages[substr(group1, 1, 1)]))
    }
    if (family != "sex") {
      stopifnot(all(as.character(comparison$Sex) == sexes[substr(group1, 3, 3)]))
    }

    metadata <- data.frame(
      Index = as.character(j),
      Comparison = paste(group1, "VS", group2),
      "First group before VS" = sample_group_number2[[group1]],
      "Second group after VS" = sample_group_number2[[group2]],
      check.names = FALSE
    )
    cbind(df0, metadata, comparison)
  })
  names(comparisons) <- definitions$Name
  comparisons
}

apply_method4 <- function(comparison) {
  count1 <- comparison[["First group before VS"]]
  count2 <- comparison[["Second group after VS"]]
  stopifnot(!anyNA(count1), !anyNA(count2))
  comparison <- comparison[count1 >= 5 | count2 >= 5, , drop = FALSE]
  comparison$p.value.adjusted <- p.adjust(comparison$p.value, method = "BH")
  comparison
}

age_comparisons_unfiltered <- build_comparisons(all_comparisons_df, age_definitions, "age")
sex_comparisons_unfiltered <- build_comparisons(all_comparisons_df2, sex_definitions, "sex")
region_comparisons_unfiltered <- build_comparisons(all_comparisons_df3, region_definitions, "region")

age.vs.young4 <- lapply(age_comparisons_unfiltered, apply_method4)
female.vs.male4 <- lapply(sex_comparisons_unfiltered, apply_method4)
region.vs.region4 <- lapply(region_comparisons_unfiltered, apply_method4)


# Export Method 4 results as combined tables and comparison-specific worksheets.
age_summary4 <- bind_rows(age.vs.young4)
sex_summary4 <- bind_rows(female.vs.male4)
region_summary4 <- bind_rows(region.vs.region4)

export(age.vs.young4, "age.vs.young4.xlsx")
export(female.vs.male4, "female.vs.male4.xlsx")
export(region.vs.region4, "region.vs.region4.xlsx")


# Count DEPs using FDR < 0.05 and absolute model estimate > log2(1.5).
# Positive estimates indicate higher abundance in the first group of a contrast.
summarize_DE <- function(comparison_list, comparison_names) {
  stopifnot(length(comparison_list) == length(comparison_names))
  counts <- lapply(comparison_list, function(comparison) {
    significant <- !is.na(comparison$p.value.adjusted) & comparison$p.value.adjusted < 0.05
    up <- significant & !is.na(comparison$estimate) & comparison$estimate > log2(1.5)
    down <- significant & !is.na(comparison$estimate) & comparison$estimate < -log2(1.5)
    matrisome <- !is.na(comparison$Matrisome.Category)
    data.frame(
      Up_Total = sum(up),
      Down_Total = sum(down),
      Up_Matrisome = sum(up & matrisome),
      Down_Matrisome = sum(down & matrisome)
    )
  })
  results <- cbind(Comparison = comparison_names, bind_rows(counts))
  results$Total_DE <- results$Up_Total + results$Down_Total
  results$Total_DE_Matrisome <- results$Up_Matrisome + results$Down_Matrisome
  results
}

age_comparison_names <- c(
  "Bottom Female", "Central Female", "Peripheral Female", "Top Female",
  "Bottom Male", "Central Male", "Peripheral Male", "Top Male"
)

age_results <- summarize_DE(age.vs.young4, age_comparison_names)

sex_comparison_names <- c(
  "Bottom Aged", "Central Aged", "Peripheral Aged", "Top Aged",
  "Bottom Young", "Central Young", "Peripheral Young", "Top Young"
)

sex_results <- summarize_DE(female.vs.male4, sex_comparison_names)

region_comparison_names <- c(
  "Bottom vs Central Aged Female", "Bottom vs Peripheral Aged Female",
  "Bottom vs Top Aged Female", "Central vs Peripheral Aged Female",
  "Central vs Top Aged Female", "Peripheral vs Top Aged Female",
  "Bottom vs Central Aged Male", "Bottom vs Peripheral Aged Male",
  "Bottom vs Top Aged Male", "Central vs Peripheral Aged Male",
  "Central vs Top Aged Male", "Peripheral vs Top Aged Male",
  "Bottom vs Central Young Female", "Bottom vs Peripheral Young Female",
  "Bottom vs Top Young Female", "Central vs Peripheral Young Female",
  "Central vs Top Young Female", "Peripheral vs Top Young Female",
  "Bottom vs Central Young Male", "Bottom vs Peripheral Young Male",
  "Bottom vs Top Young Male", "Central vs Peripheral Young Male",
  "Central vs Top Young Male", "Peripheral vs Top Young Male"
)

region_results <- summarize_DE(region.vs.region4, region_comparison_names)

results_list4 <- list(
  "Age_Comparisons" = age_results,
  "Sex_Comparisons" = sex_results,
  "Region_Comparisons" = region_results
)

export(results_list4, file = "Method 4 DE_protein_summary_all_comparisons_FDR0.05_FC1.5.xlsx")


# Generate all summary plots from Method 4 results.
# Prepare age-comparison plots.
age_comparison_names <- c(
  "Bottom Male", "Bottom Female",
  "Top Male", "Top Female",
  "Central Male", "Central Female",
  "Peripheral Male", "Peripheral Female"
)

age_plot_data <- age_results %>%
  dplyr::select(Comparison, Up_Total, Down_Total) %>%
  gather(key = "Category", value = "Number", -Comparison) %>%
  mutate(
    Direction = case_when(
      grepl("Down", Category) ~ "Up in Young",
      grepl("Up", Category) ~ "Up in Aged"
    )
  ) %>%
  mutate(Direction = factor(Direction, levels = c("Up in Young", "Up in Aged"))) %>%
  mutate(Comparison = factor(Comparison, levels = age_comparison_names))


# Make Figure 3a
p1 <- ggplot(age_plot_data, aes(x = Comparison, y = Number, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  ylim(0, 250) +
  labs(
    title = "Age Comparisons",
    x = "",
    y = "Number of Proteins",
    fill = ""
  ) +
  geom_text(aes(label = ifelse(Number == 0, "0", Number)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 3
  )


# Prepare sex-comparison plots.
sex_comparison_names <- c(
  "Bottom Young", "Bottom Aged",
  "Top Young", "Top Aged",
  "Central Young", "Central Aged",
  "Peripheral Young", "Peripheral Aged"
)

sex_plot_data <- sex_results %>%
  dplyr::select(Comparison, Up_Total, Down_Total) %>%
  gather(key = "Category", value = "Number", -Comparison) %>%
  mutate(
    Direction = case_when(
      grepl("Down", Category) ~ "Up in Male",
      grepl("Up", Category) ~ "Up in Female"
    )
  ) %>%
  mutate(Direction = factor(Direction, levels = c("Up in Male", "Up in Female"))) %>%
  mutate(Comparison = factor(Comparison, levels = sex_comparison_names))


# Make Figure 5b
p2 <- ggplot(sex_plot_data, aes(x = Comparison, y = Number, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  ylim(0, 250) +
  labs(
    title = "Sex Comparisons",
    x = "",
    y = "Number of Proteins",
    fill = ""
  ) +
  geom_text(aes(label = ifelse(Number == 0, "0", Number)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 3
  )


# Prepare region-comparison plots for each age and sex.
young_male_data <- region_results[19:24, ]
young_female_data <- region_results[13:18, ]
aged_male_data <- region_results[7:12, ]
aged_female_data <- region_results[1:6, ]


# Plot total DEPs for six region contrasts.
plot_region_comparison <- function(data, title) {
  plot_data <- data %>%
    dplyr::select(Comparison, Up_Total, Down_Total) %>%
    gather(key = "Category", value = "Number", -Comparison) %>%
    mutate(
      Direction = case_when(
        grepl("Down", Category) ~ "Up in Second Region",
        grepl("Up", Category) ~ "Up in First Region"
      )
    ) %>%
    mutate(Direction = factor(Direction, levels = c("Up in Second Region", "Up in First Region"))) %>%
    mutate(Comparison = gsub(
      " vs ", "\nvs\n",
      sub(" (Aged|Young) (Female|Male)$", "", Comparison)
    ))

  ggplot(plot_data, aes(x = Comparison, y = Number, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top",
      axis.text = element_text(size = 8),
      plot.margin = margin(b = 10, t = 10)
    ) +
    ylim(0, 250) +
    labs(
      title = title,
      x = "",
      y = "Number of Proteins",
      fill = ""
    ) +
    geom_text(aes(label = ifelse(Number == 0, "0", Number)),
      position = position_dodge(width = 0.7),
      vjust = -0.5,
      size = 3
    )
}


# Make Figure S15a,b,c,d
p3 <- plot_region_comparison(young_male_data, "Region Comparisons in Young Male")
p4 <- plot_region_comparison(young_female_data, "Region Comparisons in Young Female")
p5 <- plot_region_comparison(aged_male_data, "Region Comparisons in Aged Male")
p6 <- plot_region_comparison(aged_female_data, "Region Comparisons in Aged Female")


# Save the six total-protein plots together. Make Figure 3a, 5b, S15a,b,c,d
pdf("DE_protein_summary_plots.pdf", width = 12, height = 18)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
dev.off()


# Prepare age-comparison plots.
age_comparison_names <- c(
  "Bottom Male", "Bottom Female",
  "Top Male", "Top Female",
  "Central Male", "Central Female",
  "Peripheral Male", "Peripheral Female"
)


# Prepare the matrisome age-comparison plot.
age_plot_data_matrisome <- age_results %>%
  dplyr::select(Comparison, Up_Matrisome, Down_Matrisome) %>%
  gather(key = "Category", value = "Number", -Comparison) %>%
  mutate(
    Direction = case_when(
      grepl("Down", Category) ~ "Up in Young",
      grepl("Up", Category) ~ "Up in Aged"
    )
  ) %>%
  mutate(Direction = factor(Direction, levels = c("Up in Young", "Up in Aged"))) %>%
  mutate(Comparison = factor(Comparison, levels = age_comparison_names))


# Make Figure 3b
p7 <- ggplot(age_plot_data_matrisome, aes(x = Comparison, y = Number, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  ylim(0, 50) +
  labs(
    title = "Age Comparisons - Matrisome Proteins",
    x = "",
    y = "Number of Proteins",
    fill = ""
  ) +
  geom_text(aes(label = ifelse(Number == 0, "0", Number)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 3
  )


# Prepare sex-comparison plots.
sex_comparison_names <- c(
  "Bottom Young", "Bottom Aged",
  "Top Young", "Top Aged",
  "Central Young", "Central Aged",
  "Peripheral Young", "Peripheral Aged"
)


# Prepare the matrisome sex-comparison plot.
sex_plot_data_matrisome <- sex_results %>%
  dplyr::select(Comparison, Up_Matrisome, Down_Matrisome) %>%
  gather(key = "Category", value = "Number", -Comparison) %>%
  mutate(
    Direction = case_when(
      grepl("Down", Category) ~ "Up in Male",
      grepl("Up", Category) ~ "Up in Female"
    )
  ) %>%
  mutate(Direction = factor(Direction, levels = c("Up in Male", "Up in Female"))) %>%
  mutate(Comparison = factor(Comparison, levels = sex_comparison_names))


# Make Figure 5c
p8 <- ggplot(sex_plot_data_matrisome, aes(x = Comparison, y = Number, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  ylim(0, 50) +
  labs(
    title = "Sex Comparisons - Matrisome Proteins",
    x = "",
    y = "Number of Proteins",
    fill = ""
  ) +
  geom_text(aes(label = ifelse(Number == 0, "0", Number)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 3
  )


# Plot matrisome DEPs for six region contrasts.
plot_region_comparison_matrisome <- function(data, title) {
  plot_data <- data %>%
    dplyr::select(Comparison, Up_Matrisome, Down_Matrisome) %>%
    gather(key = "Category", value = "Number", -Comparison) %>%
    mutate(
      Direction = case_when(
        grepl("Down", Category) ~ "Up in Second Region",
        grepl("Up", Category) ~ "Up in First Region"
      )
    ) %>%
    mutate(Direction = factor(Direction, levels = c("Up in Second Region", "Up in First Region"))) %>%
    mutate(Comparison = gsub(
      " vs ", "\nvs\n",
      sub(" (Aged|Young) (Female|Male)$", "", Comparison)
    ))

  ggplot(plot_data, aes(x = Comparison, y = Number, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top",
      axis.text = element_text(size = 8),
      plot.margin = margin(b = 10, t = 10)
    ) +
    ylim(0, 50) +
    labs(
      title = paste0(title, " - Matrisome Proteins"),
      x = "",
      y = "Number of Proteins",
      fill = ""
    ) +
    geom_text(aes(label = ifelse(Number == 0, "0", Number)),
      position = position_dodge(width = 0.7),
      vjust = -0.5,
      size = 3
    )
}


# Make Figure S15e,f,g,h
p9 <- plot_region_comparison_matrisome(young_male_data, "Region Comparisons in Young Male")
p10 <- plot_region_comparison_matrisome(young_female_data, "Region Comparisons in Young Female")
p11 <- plot_region_comparison_matrisome(aged_male_data, "Region Comparisons in Aged Male")
p12 <- plot_region_comparison_matrisome(aged_female_data, "Region Comparisons in Aged Female")


# Save the six matrisome plots together. Make Figure 3b, 5c, S15e,f,g,h
pdf("DE_matrisome_protein_summary_plots.pdf", width = 12, height = 18)
grid.arrange(p7, p8, p9, p10, p11, p12, ncol = 2)
dev.off()


# Plot volcano plots for the eight multifactor aged-versus-young comparisons.

# 10 proteins upregulated in aged groups across all eight comparisons.
example_up <- c(
  "Ftl1", "Clu", "Fth1", "Slc4a1", "Fmo1",
  "Flt1", "Spp2", "Bst1", "Xpnpep2", "Hamp"
)

# 11 proteins upregulated in young groups across all eight comparisons.
example_down <- c(
  "Loxl2", "Pcolce", "Ltbp2", "Col3a1", "Serpinh1",
  "Ikbip", "Rbms1", "Fam151a", "Plod2", "Eln", "Fbln5"
)

example_combined <- c(example_down, example_up)


# Labels follow the order of the eight aged-versus-young comparisons.
age_vs_young_index_all <- c(
  "Bottom Female", "Central Female", "Peripheral Female", "Top Female",
  "Bottom Male", "Central Male", "Peripheral Male", "Top Male"
)


# Make Figure 3d,e,f,g and Figure S5a,b,c,d
pdf("age_vs_young_DEP_multiple_factor_volcano_plot.pdf", width = 12, height = 8)

for (i in 1:8) {
  age_DE <- age.vs.young4[[i]]

  res1 <- age_DE
  res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)
  res2 <- res1[res1$p.value.adjusted < 0.05, ]
  res2 <- res2[abs(res2$estimate) > 0.5849625, ]
  res3 <- res2[res2$Gene.symbol %in% example_combined, ]

  fc_cut <- 1.5 # Set the fold-change cutoff.

  pp1 <- ggplot(res1, aes(estimate, LogPValue))
  pp1 <- pp1 + geom_point(shape = 1, color = "gray") + coord_cartesian(xlim = c(-4, 4), ylim = c(-0.5, 20))
  pp1 <- pp1 + geom_point(
    data = res2,
    aes(x = estimate, y = LogPValue, colour = estimate),
    shape = 19,
    size = 2
  ) + scale_colour_gradient2(
    low = "green",
    mid = "lightyellow",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) #+ scale_colour_gradient(high="red",low = "green")
  pp1 <- pp1 + geom_point(data = res3, aes(x = estimate, y = LogPValue), colour = "black", shape = 1, size = 2)
  pp1 <- pp1 + geom_text_repel(
    data = res3,
    aes(x = estimate, y = LogPValue, label = Gene.symbol),
    size = 4.5,
    arrow = arrow(length = unit(0.01, 'npc')),
    point.padding = unit(0.6, 'lines'),
    max.overlaps = 1600
  )
  pp1 <- pp1 + geom_hline(yintercept = -log10(0.05), colour = "red", size = 0.5, linetype = 2)
  pp1 <- pp1 + geom_vline(xintercept = 0.5849625, colour = "red", size = 0.5, linetype = 2)
  pp1 <- pp1 + geom_vline(xintercept = -0.5849625, colour = "red", size = 0.5, linetype = 2)

  pp1 <- pp1 + theme_bw(base_size = 16) + labs(x = "estimate", y = "-Log10(FDR)") + ggtitle(age_vs_young_index_all[i])
  plot(pp1)
}
dev.off()


# Summarize DEPs across the eight pairwise aged-versus-young comparisons.
# Collect the union of 670 DEPs.
i = 1
gene_set <- c()
for (i in 1:8) {
  age_DE <- age.vs.young4[[i]]

  res1 <- age_DE
  res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)
  # res2 <- res1[res1$estimate>1|res1$estimate<(-1),]
  res2 <- res1[res1$p.value.adjusted < 0.05, ]
  res2 <- res2[abs(res2$estimate) > 0.5849625, ]
  gene_set <- union(gene_set, res2$Uniprot.ID)
}
# The union contains 670 proteins.

i = 1
summary_table <- list()
for (i in 1:8) {
  age_DE <- age.vs.young4[[i]]

  res1 <- age_DE
  res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)
  # res2 <- res1[res1$estimate>1|res1$estimate<(-1),]
  res2 <- res1[res1$p.value.adjusted < 0.05, ]
  res2 <- res2[abs(res2$estimate) > 0.5849625, ]
  res3 <- res2[c(1, 15)]

  summary_table[[i]] <- data.frame("index" = seq_along(gene_set), "Uniprot" = gene_set, "estimate" = 0, stringsAsFactors = FALSE)

  # A zero in the summary table indicates either exclusion by filtering
  # or no significant differential expression between aged and young groups.

  for (j in seq_along(gene_set)) {
    if (summary_table[[i]]$Uniprot[j] %in% res3$Uniprot.ID) {
      index_number <- which(res3$Uniprot == summary_table[[i]]$Uniprot[j])
      summary_table[[i]]$estimate[j] <- res3$estimate[index_number]
    }
  }
}

summary_table1 <- cbind(
  summary_table[[1]], summary_table[[2]]$estimate,
  summary_table[[3]]$estimate,
  summary_table[[4]]$estimate,
  summary_table[[5]]$estimate,
  summary_table[[6]]$estimate,
  summary_table[[7]]$estimate,
  summary_table[[8]]$estimate
)

colnames(summary_table1)[c(3:10)] <- age_vs_young_index_all
rownames(summary_table1) <- summary_table1$Uniprot
summary_table1 <- summary_table1[-c(1, 2)]
summary_table2 <- ifelse(summary_table1 > 0, 1, ifelse(summary_table1 < 0, -1, 0))
summary_table2 <- as.data.frame(summary_table2)
summary_table3 <- cbind(rownames(summary_table2), summary_table2)

matrisome_anno_table1 <- matrisome_anno_table[matrisome_anno_table$Uniprot.ID %in% summary_table3$`rownames(summary_table2)`, ]
matrisome_anno_table1 <- matrisome_anno_table1[match(summary_table3$`rownames(summary_table2)`, matrisome_anno_table1$Uniprot.ID), ]

export(summary_table3, "summary_table3.xlsx")

# Plot the 38 selected matrisome DEPs across eight comparisons.

matrisome_dep_list <- read_delim("38 matrisome DEPs with unique features in 8 pairwise comparisons.txt",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)
# summary_table1 contains log2 fold changes for 670 proteins across eight comparisons.
# summary_table4 contains the corresponding values for 38 selected matrisome proteins.
# Zero indicates exclusion by filtering or no significant differential expression.
summary_table4 <- summary_table1[match(matrisome_dep_list$Uniprot.ID, rownames(summary_table1)), ]
summary_table4 <- summary_table4[c(1, 5, 2, 6, 3, 7, 4, 8)] # Column order: BF, BM, CF, CM, PF, PM, TF, TM.
rownames(summary_table4) <- matrisome_dep_list$Gene.symbol


# Plot the selected matrisome log2 fold changes.
# Zero indicates exclusion by the second filter or no significant difference
# at the fold-change and FDR thresholds of 1.5 and 0.05.
# Discrete color scale with a cutoff at |log2 FC| = log2(1.5) = 0.585.
# Bin edges are the 11 equal-width bins over the data range used by the
# original pheatmap default (range -2.13 to 1.77 in the reported run);
# the two middle bins are replaced by one white bin from -0.585 to 0.585
# (the zeros: not a DEP or filtered out). Result: 5 blue bins (higher in
# Young = 7w), white, 4 red bins (higher in Aged = 7-11m).
heatmap_fc_cut <- log2(1.5)
heatmap_matrix <- as.matrix(summary_table4)
heatmap_bins <- seq(min(heatmap_matrix), max(heatmap_matrix), length.out = 12)
heatmap_breaks <- c(heatmap_bins[1:5], -heatmap_fc_cut, heatmap_fc_cut, heatmap_bins[9:12])
rdbu <- rev(brewer.pal(n = 11, name = "RdBu")) # dark blue ... white ... dark red
heatmap_colors <- c(rdbu[1:5], "white", rdbu[8:11])

# Column annotation bars: lung region (B / C / P / T) and sex (F / M).
# pheatmap draws the last column on top, so Sex is listed first to place
# the Region bar above the Sex bar.
heatmap_col_anno <- data.frame(
  Sex = factor(c("F", "M", "F", "M", "F", "M", "F", "M"), levels = c("F", "M")),
  Region = factor(c("B", "B", "C", "C", "P", "P", "T", "T"), levels = c("B", "C", "P", "T")),
  row.names = colnames(summary_table4)
)
heatmap_anno_colors <- list(
  Region = c(B = "#F47F26", C = "#A5562D", P = "#A9C0D9", T = "#874F9F"),
  Sex = c(F = "#CA2435", M = "#485AA7")
)

# Gene symbols in upper case, as in the figure.
rownames(heatmap_matrix) <- toupper(rownames(heatmap_matrix))

pppp <- pheatmap(
  heatmap_matrix,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  legend_breaks = c(-2, -1, -heatmap_fc_cut, heatmap_fc_cut, 1),
  legend_labels = c("-2", "-1", "-0.585", "0.585", "1"),
  annotation_col = heatmap_col_anno,
  annotation_colors = heatmap_anno_colors,
  annotation_names_col = FALSE, # Hide the "Region" / "Sex" labels next to the bars.
  fontsize_row = 8, # Row font size
  show_rownames = TRUE, # Show row names.
  show_colnames = FALSE, # Columns are identified by the annotation bars.
  cluster_cols = FALSE, # Keep the original column order (BF, BM, CF, ... TM).
  cluster_rows = FALSE, # Keep the original row order.
  cellwidth = 14, cellheight = 9, # Fixed cell size.
  border_color = "grey50", # Draw grey cell borders.
  main = "Matrisome DEPs, 7-11M vs 7W" # Main title of the heatmap
)


# Make Figure 3h
print(pppp)

# Save the heatmap (Figure 3h).
pdf("Age_vs_young_matrisome_DEP_heatmap.pdf", width = 11, height = 8.5)
print(pppp)
dev.off()


# Plot the 38 selected matrisome proteins using imputed log2 abundances.
# Reuse the protein list and annotated intensity table loaded above.
boxplot_intensities <- as.data.frame(
  data_imp[match(matrisome_dep_list$Uniprot.ID, data_imp$Uniprot.ID), 8:103]
)
rownames(boxplot_intensities) <- matrisome_dep_list$Gene.symbol


# Use the original sample order, with young and aged samples adjacent.
boxplot_sample_order <- read.table("Sample_group_order.txt", header = TRUE)
boxplot_intensities <- boxplot_intensities[, boxplot_sample_order$Sample.ID]
boxplot_groups <- c(
  "YBF", "ABF", "YBM", "ABM", "YCF", "ACF", "YCM", "ACM",
  "YPF", "APF", "YPM", "APM", "YTF", "ATF", "YTM", "ATM"
)


# Save one protein per page, showing all six biological replicates per group.
# Make Figure 5h, 5i, S6a,b,c,d
pdf("All_Proteins_Boxplots.pdf", width = 11.5, height = 8)

for (i in seq_len(nrow(boxplot_intensities))) {
  boxplot_data <- data.frame(
    Group = factor(boxplot_sample_order$Sample.group, levels = boxplot_groups),
    Value = as.numeric(boxplot_intensities[i, ])
  )
  boxplot_data$Age <- ifelse(substr(boxplot_data$Group, 1, 1) == "Y", "Young", "Aged")

  pp2 <- ggplot(boxplot_data, aes(x = Group, y = Value)) +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_point(
      aes(colour = Age),
      position = position_jitter(width = 0.2, height = 0, seed = 1),
      alpha = 0.8
    ) +
    scale_colour_manual(values = c(Young = "#3AAA69", Aged = "#BC5DA3")) +
    coord_cartesian(ylim = range(boxplot_data$Value, na.rm = TRUE) + c(-0.2, 0.2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = rownames(boxplot_intensities)[i],
      x = "Group",
      y = "Log2 protein abundance"
    )

  print(pp2)
}

dev.off()


# Calculate correlations across the 670-protein DEP direction matrix.
# Preserve the original input: -1 = up in young, 1 = up in aged, and
# 0 = not significant or excluded by the comparison-specific filter.
# These correlations describe DEP directions, not continuous log2 fold changes.
correlation_input <- as.matrix(summary_table2)
cor_summary <- cor(correlation_input, method = "pearson")
cor_summary1 <- cor(correlation_input, method = "spearman")
cor_summary_long <- reshape2::melt(cor_summary) # using Pearson in bubble plot

pp3 <- ggplot(cor_summary_long, aes(x = Var1, y = Var2, size = value, color = value)) +
  geom_point(alpha = 1) +
  scale_size(
    range = c(1, 10),
    limits = c(0.2, 1.0),
    breaks = seq(0.2, 1.0, by = 0.2),
    labels = c("0.2", "0.4", "0.6", "0.8", "1.0")
  ) +
  scale_color_viridis_c(
    option = "inferno",
    limits = c(0.2, 1.0),
    breaks = seq(0.2, 1.0, by = 0.2),
    labels = c("0.2", "0.4", "0.6", "0.8", "1.0")
  ) +
  labs(x = "Row", y = "Column", size = "Correlation", color = "Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_fixed()


# Make Figure 3c.
print(pp3)
ggsave("Age_vs_young_DEP_correlation_map.pdf", pp3, width = 8, height = 7)

# save the current workspace. 
save.image("Script_6_multiple_factor_DE.RData")
