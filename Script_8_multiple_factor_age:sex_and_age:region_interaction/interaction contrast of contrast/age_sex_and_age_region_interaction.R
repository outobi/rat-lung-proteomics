# =============================================================================
# Script 8 (interaction contrast of contrast): Age x Sex and Age x Region
#   interaction contrasts (multiple-factor model)
# =============================================================================
#
# Purpose:
#   Test whether the age effect (Aged - Young; Aged = young adult,
#   Young = peripubertal) differs between sexes or between lung regions.
#   For each of the 3,600 filtered and imputed protein groups from Script_3,
#   fit the same linear mixed model as Script_6:
#       Protein ~ Age + Sex + Region + Age:Region + Sex:Region + Age:Sex
#                 + (1 | Subject)
#   and compute two "contrast of contrasts" with emmeans:
#     - Age x Sex:    (Aged - Young)Female - (Aged - Young)Male,
#                     averaged over the 4 regions (1 contrast per protein).
#     - Age x Region: (Aged - Young)Region1 - (Aged - Young)Region2,
#                     averaged over sex, for the 6 region pairs
#                     (Bottom - Central, Bottom - Peripheral, Bottom - Top,
#                      Central - Peripheral, Central - Top, Peripheral - Top).
#   p-values are unadjusted per protein and BH-adjusted: across the 3,600
#   proteins for Age x Sex, and across all 21,600 Age x Region contrasts
#   together (3,600 proteins x 6 region pairs) for Age x Region.
#   No observed-value filter is applied (all 3,600 proteins).
#   Significant: FDR < 0.05 and |estimate| > log2(1.5) = 0.585.
#   Positive estimate = larger age effect in Female (Age x Sex) or in
#   Region1 (Age x Region).
#
# Inputs (all in this folder):
#   - UW.Madison.rat.lung.proteomics_data_tables.RData  (output of Script_3)
#       data_table$data_imp: imputed, normalized log2 intensities
#       (3,600 proteins x S1-S96); info_table: sample metadata in S1-S96
#       order.
#   - UW.Madison.rat.lung.proteomics_DIANN_protein_intensity.xlsx
#       (output of Script_3) Sheet "data_imp": the 7 annotation columns are
#       attached to the results (rows in the same order as data_imp).
#
# Outputs (written to this folder):
#   - age sex interaction protein list.xlsx
#       Age x Sex contrast for all 3,600 proteins (annotations, estimate, SE,
#       df, t.ratio, p.value, FDR). Reported run: 23 larger in Female,
#       32 larger in Male.
#   - age region interaction protein list.xlsx
#       Age x Region contrasts, one sheet per region pair (3,600 proteins
#       each). Reported run (Region1 / Region2 larger): Bottom-Central 30/30,
#       Bottom-Peripheral 28/35, Bottom-Top 40/30, Central-Peripheral 31/39,
#       Central-Top 40/28, Peripheral-Top 49/27.
#   - age sex interaction volcano plot.pdf
#   - age region interaction volcano plot <region1> vs <region2>.pdf  (6 files)
#   - Age_region_interaction_volcano_summary_2x3.pdf  (the 6 region plots)
#   - data_summary.RData
#       Result tables, summaries and plot objects (not the whole workspace;
#       the combined 2 x 3 plot object is not saved because of its size).
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# Model fitting takes several minutes (3,600 lmer fits).
# =============================================================================

# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages.
# readxl: read_excel(); rio: export(); openxlsx: write.xlsx() (multi-sheet);
# lme4 + lmerTest: lmer(); emmeans: emmeans() / contrast();
# dplyr: bind_rows() / group_by() / summarise();
# ggplot2 + ggrepel: volcano plots; patchwork: combined 2 x 3 figure.
library(readxl)
library(rio)
library(openxlsx)
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)


# Step 1: Load data and build the model table ----------------------------------

# Load the saved matrices (data_table) and sample information (info_table).
load("UW.Madison.rat.lung.proteomics_data_tables.RData")

# data7: 3,600 proteins (rows, Uniprot.ID) x 96 samples (S1-S96).
data7 <- data_table$data_imp

# file02: the 7 annotation columns, attached to the results by row position.
file02 <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein_intensity.xlsx",
  sheet = "data_imp"
)[, 1:7]
stopifnot(identical(file02$Uniprot.ID, rownames(data7)))

# Model table from info_table (96 rows in S1-S96 order, the same order as the
# columns of data7). Subject = animal ID (each animal contributes 4 regions);
# Protein is filled with one protein per loop iteration.
workdata <- data.frame(
  Subject = info_table$Sample.PID,
  Age = as.character(info_table$Sample.Age),
  Region = as.character(info_table$Sample.Region),
  Sex = as.character(info_table$Sample.Sex),
  Protein = NA_real_
)


# Step 2: Fit the model and compute both interaction contrasts -----------------

# The same model is used for both contrasts, so it is fitted once per protein.
fc_cut <- log2(1.5) # 0.585; |estimate| cutoff
all_interactions_age_sex <- vector("list", nrow(data7))
all_interactions_age_region <- vector("list", nrow(data7))

for (i in seq_len(nrow(data7))) {
  workdata$Protein <- as.numeric(data7[i, ])

  res <- lmer(
    Protein ~ Age + Sex + Region +
      Age:Region + Sex:Region + Age:Sex +
      (1 | Subject),
    data = workdata
  )

  # Age x Sex: estimated marginal means for Age x Sex, averaged over the
  # four regions; contrast = (Aged - Young in Female) - (Aged - Young in Male).
  emm_age_sex <- emmeans(res, ~ Age * Sex)
  interaction_age_sex <- contrast(emm_age_sex, interaction = "pairwise", adjust = "none")
  all_interactions_age_sex[[i]] <- as.data.frame(summary(interaction_age_sex))

  # Age x Region: estimated marginal means for Age x Region, averaged over
  # sex; contrasts compare the age effect between each pair of regions.
  emm_age_region <- emmeans(res, ~ Age * Region)
  interaction_age_region <- contrast(emm_age_region, interaction = "pairwise", adjust = "none")
  all_interactions_age_region[[i]] <- as.data.frame(summary(interaction_age_region))
}


# Step 3: Age x Sex results ----------------------------------------------------

# Combine contrasts and apply BH correction across the 3,600 proteins.
age_sex_interaction_df <- bind_rows(all_interactions_age_sex)
age_sex_interaction_df$FDR <- p.adjust(age_sex_interaction_df$p.value, method = "BH")
rownames(age_sex_interaction_df) <- rownames(data7)

# Count significant proteins (FDR < 0.05 and |estimate| > 0.585).
age_sex_sig_summary <- age_sex_interaction_df %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::summarise(
    Female_greater = sum(estimate > fc_cut, na.rm = TRUE),
    Male_greater   = sum(estimate < -fc_cut, na.rm = TRUE),
    Total          = dplyr::n() # all FDR < 0.05, without the estimate cutoff
  )
age_sex_sig_summary
# Reported run: Female_greater 23, Male_greater 32, Total 55

# Attach annotations and export.
age_sex_interaction_df <- cbind(file02, age_sex_interaction_df)
export(age_sex_interaction_df, "age sex interaction protein list.xlsx")


# Step 4: Age x Region results -------------------------------------------------

# Combine contrasts and apply BH correction across all 21,600 Age x Region
# contrasts together (3,600 proteins x 6 region pairs), as in the reported
# analysis. (The original code grouped by region pair, but with plyr attached
# its mutate() ignored the grouping, so the reported FDRs are across all
# contrasts; this is now written explicitly.)
age_region_interaction_df <- bind_rows(all_interactions_age_region)
age_region_interaction_df$FDR <- p.adjust(age_region_interaction_df$p.value, method = "BH")

# Count significant proteins per region pair.
age_region_sig_summary <- age_region_interaction_df %>%
  dplyr::group_by(Region_pairwise) %>%
  dplyr::summarise(
    n_tested = dplyr::n(),
    n_FDR_005 = sum(FDR < 0.05 & abs(estimate) > fc_cut, na.rm = TRUE),
    Region1_greater = sum(FDR < 0.05 & estimate > fc_cut, na.rm = TRUE),
    Region2_greater = sum(FDR < 0.05 & estimate < -fc_cut, na.rm = TRUE)
  )
age_region_sig_summary
# Reported run:
#   Region_pairwise      n_tested n_FDR_005 Region1_greater Region2_greater
#   Bottom - Central         3600        60              30              30
#   Bottom - Peripheral      3600        63              28              35
#   Bottom - Top             3600        70              40              30
#   Central - Peripheral     3600        70              31              39
#   Central - Top            3600        68              40              28
#   Peripheral - Top         3600        76              49              27

# Split into one table per region pair (sheet names "Bottom_vs_Central", ...),
# attach annotations (rows within each pair are in data7 order) and export.
age_region_list <- split(age_region_interaction_df, age_region_interaction_df$Region_pairwise)
names(age_region_list) <- gsub(" - ", "_vs_", names(age_region_list))
age_region_list <- lapply(age_region_list, function(x) {
  cbind(file02, x)
})

openxlsx::write.xlsx(
  age_region_list,
  file = "age region interaction protein list.xlsx",
  rowNames = FALSE
)


# Step 5: Volcano plots --------------------------------------------------------
# Each plot: all proteins (gray open circles); significant proteins
# (FDR < 0.05 and |estimate| > 0.585) colored by estimate (red = positive,
# green = negative); significant matrisome proteins outlined in black and
# labelled. Dashed red lines: FDR = 0.05 and estimate = +/- 0.585.
# The title gives the number of significant proteins in each direction.

plot_interaction_volcano <- function(res1, title_prefix, group1, group2) {
  res1$LogPValue <- (-1) * log10(res1$FDR)
  res2 <- res1[res1$FDR < 0.05, ]
  res2 <- res2[abs(res2$estimate) > fc_cut, ]
  res3 <- res2[!is.na(res2$Matrisome.Sub.category), ] # Matrisome proteins.
  n_up <- sum(res2$estimate > fc_cut)
  n_down <- sum(res2$estimate < -fc_cut)

  p <- ggplot(res1, aes(estimate, LogPValue))
  p <- p + geom_point(shape = 1, color = "gray")
  p <- p + geom_point(data = res2, aes(x = estimate, y = LogPValue, colour = estimate), shape = 19, size = 2) +
    scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")
  p <- p + geom_point(data = res3, aes(x = estimate, y = LogPValue), colour = "black", shape = 1, size = 2)
  p <- p + geom_text_repel(data = res3, aes(x = estimate, y = LogPValue, label = Gene.symbol), size = 4.5, arrow = arrow(length = unit(0.01, "npc")), point.padding = unit(0.6, "lines"), max.overlaps = 1600)
  p <- p + geom_hline(yintercept = -log10(0.05), colour = "red", size = 0.5, linetype = 2)
  p <- p + geom_vline(xintercept = fc_cut, colour = "red", size = 0.5, linetype = 2)
  p <- p + geom_vline(xintercept = -fc_cut, colour = "red", size = 0.5, linetype = 2)
  p + theme_bw(base_size = 16) + labs(x = "estimate", y = "-Log10(FDR)") +
    ggtitle(paste0(title_prefix, " interaction volcano plot, up in red ", group1, ", n = ", n_up,
                   "; up in green ", group2, ", n = ", n_down))
}

# Age x Sex volcano plot.
p_age_sex <- plot_interaction_volcano(age_sex_interaction_df, "Age×Sex", "female", "male")
plot(p_age_sex)
ggsave(filename = "age sex interaction volcano plot.pdf", plot = p_age_sex, device = "pdf", width = 11, height = 8.5)

# Age x Region volcano plots, one per region pair (p1-p6).
p1 <- plot_interaction_volcano(age_region_list$Bottom_vs_Central, "Age×Region", "bottom", "central")
p2 <- plot_interaction_volcano(age_region_list$Bottom_vs_Peripheral, "Age×Region", "bottom", "peripheral")
p3 <- plot_interaction_volcano(age_region_list$Bottom_vs_Top, "Age×Region", "bottom", "top")
p4 <- plot_interaction_volcano(age_region_list$Central_vs_Peripheral, "Age×Region", "central", "peripheral")
p5 <- plot_interaction_volcano(age_region_list$Central_vs_Top, "Age×Region", "central", "top")
p6 <- plot_interaction_volcano(age_region_list$Peripheral_vs_Top, "Age×Region", "peripheral", "top")

region_volcano_files <- c(
  "age region interaction volcano plot bottom vs central.pdf",
  "age region interaction volcano plot bottom vs peripheral.pdf",
  "age region interaction volcano plot bottom vs top.pdf",
  "age region interaction volcano plot central vs peripheral.pdf",
  "age region interaction volcano plot central vs top.pdf",
  "age region interaction volcano plot peripheral vs top.pdf"
)
region_volcano_plots <- list(p1, p2, p3, p4, p5, p6)
for (k in seq_along(region_volcano_plots)) {
  plot(region_volcano_plots[[k]])
  ggsave(filename = region_volcano_files[k], plot = region_volcano_plots[[k]], device = "pdf", width = 11, height = 8.5)
}

# Combine the 6 region plots in a 2 x 3 layout.
p_summary <-
  (p1 | p2 | p3) /
  (p4 | p5 | p6)
p_summary

ggsave(
  filename = "Age_region_interaction_volcano_summary_2x3.pdf",
  plot = p_summary,
  device = pdf,
  width = 39,
  height = 26
)


# Step 6: Save results ---------------------------------------------------------

# Save the result tables, summaries and individual plot objects.
# (p_summary is not saved: it holds 6 copies of the data, ~110 MB.)
save(
  age_sex_interaction_df, age_sex_sig_summary,
  age_region_interaction_df, age_region_sig_summary, age_region_list,
  p_age_sex, p1, p2, p3, p4, p5, p6,
  file = "data_summary.RData"
)

