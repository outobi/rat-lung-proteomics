# Load required packages.
library(limma)
library(rio)
library(tidyr)
library(readxl)
library(lme4)
library(lmerTest)
library(reshape2)
library(emmeans)
library(Rmisc)
library(stringr)
library(dplyr)
library(multcomp)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)


# Differential abundance analysis of age, sex, and regional main effects.
# Model: additive fixed effects with a subject-specific random intercept.
# No interaction terms or subgroup-specific models are fitted.
# Inputs: saved imputed data, annotated intensity workbook, and sample metadata.
# Outputs: comparison tables, volcano plots, count bar plots, a heatmap,
# an effect-size bar plot, and an R workspace.
# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load the saved imputed matrices and sample information.
load("UW.Madison.rat.lung.proteomics_data tables_2026-09-15.RData")


# Read the annotated imputed-intensity worksheet for output annotations.
# The first seven columns are assumed to contain protein annotations.
data_imp <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein intensity_2026-09-15.xlsx",
  sheet = "data_imp"
)


# Use the saved imputed matrix for model fitting.
# Protein rows must match the annotation worksheet in the same order.
data7 <- data_table$data_imp


# Extract experimental variables from the saved sample information.
# The models below use the separately imported workdata table.
f <- info_table$Sample.Group
Region <- info_table$Sample.Region
Age <- info_table$Sample.Age
Sex <- info_table$Sample.Sex
ID <- as.factor(info_table$Sample.PID)


# Read model metadata in the same order as the sample columns in data7.
# Column 5 must be the Protein response column; Subject identifies each animal.
workdata <- read_excel("Sample_metadata.xlsx")


# Age main-effect analysis -----------------------------------------------
# Fit the additive model for each of 3,600 proteins, adjusting for sex and region.
i = 1
age_contrast_list <- list()
for (i in 1:3600) {
  # Assign the current protein intensities to the response column.
  workdata[5] <- t(data7[i, ])

  # Account for repeated regional measurements within each subject.
  res <- lmer(Protein ~ Age + Sex + Region + (1 | Subject), data = workdata)

  # Request marginal means and treatment-versus-control contrasts.
  age_emmeans_result <- emmeans(res, spec = trt.vs.ctrl ~ Age, adjust = "none")

  # Apply the pairwise contrast operation to the returned emmeans object.
  age_pairwise_contrasts <- contrast(age_emmeans_result, method = "pairwise")

  # Store the contrast summary for this protein.
  age_contrast_list[[i]] <- as.data.frame(summary(age_pairwise_contrasts))
}


# Combine age contrasts and apply Benjamini-Hochberg correction across proteins.
age_contrast_results <- bind_rows(age_contrast_list)
age_contrast_results$p.value.adjusted <- p.adjust(age_contrast_results$p.value, method = "BH")
row.names(age_contrast_results) <- rownames(data7)


# Append protein annotations by row position to the age results.
# The output contrast labels determine the direction of comparison.
age.vs.young <- cbind(data_imp[c(1:7)], age_contrast_results)


# Sex main-effect analysis -----------------------------------------------
# Fit the additive model for each of 3,600 proteins, adjusting for age and region.
i = 1
sex_contrast_list <- list()
for (i in 1:3600) {
  # Assign the current protein intensities to the response column.
  workdata[5] <- t(data7[i, ])

  # Account for repeated regional measurements within each subject.
  res <- lmer(Protein ~ Age + Sex + Region + (1 | Subject), data = workdata)

  # Request marginal means and treatment-versus-control contrasts.
  sex_emmeans_result <- emmeans(res, spec = trt.vs.ctrl ~ Sex, adjust = "none")

  # Apply the pairwise contrast operation to the returned emmeans object.
  sex_pairwise_contrasts <- contrast(sex_emmeans_result, method = "pairwise")

  # Store the contrast summary for this protein.
  sex_contrast_list[[i]] <- as.data.frame(summary(sex_pairwise_contrasts))
}


# Combine sex contrasts and apply Benjamini-Hochberg correction across proteins.
sex_contrast_results <- bind_rows(sex_contrast_list)
sex_contrast_results$p.value.adjusted <- p.adjust(sex_contrast_results$p.value, method = "BH")
row.names(sex_contrast_results) <- rownames(data7)


# Append protein annotations by row position to the sex results.
female.vs.male <- cbind(data_imp[c(1:7)], sex_contrast_results)


# Regional main-effect analysis ------------------------------------------
# Fit the additive model for each of 3,600 proteins, adjusting for age and sex.
i = 1
region_contrast_list <- list()
for (i in 1:3600) {
  # Assign the current protein intensities to the response column.
  workdata[5] <- t(data7[i, ])

  # Account for repeated regional measurements within each subject.
  res <- lmer(Protein ~ Age + Sex + Region + (1 | Subject), data = workdata)

  # Request marginal means and treatment-versus-control contrasts.
  region_emmeans_result <- emmeans(res, trt.vs.ctrl ~ Region, adjust = "none")

  # Apply the pairwise contrast operation to the returned emmeans object.
  region_pairwise_contrasts <- contrast(region_emmeans_result, method = "pairwise")

  # Store the contrast summary for this protein.
  region_contrast_list[[i]] <- as.data.frame(summary(region_pairwise_contrasts))
}


# Combine regional contrasts from all proteins.
region_contrast_results <- bind_rows(region_contrast_list)


# Extract bottom versus central results and adjust across proteins.
# Row selection assumes six contrasts per protein in the order used below.
# Verify this order against the returned contrast labels and factor levels.
sequence1 <- seq(from = 1, to = 3600 * 6, by = 6)
comparison.bottom.vs.central <- region_contrast_results[sequence1, ]
comparison.bottom.vs.central$p.value.adjusted <- p.adjust(comparison.bottom.vs.central$p.value, method = "BH")
row.names(comparison.bottom.vs.central) <- rownames(data7)


# Extract bottom versus peri results and adjust across proteins.
sequence1 <- seq(from = 2, to = 3600 * 6, by = 6)
comparison.bottom.vs.peri <- region_contrast_results[sequence1, ]
comparison.bottom.vs.peri$p.value.adjusted <- p.adjust(comparison.bottom.vs.peri$p.value, method = "BH")
row.names(comparison.bottom.vs.peri) <- rownames(data7)


# Extract bottom versus top results and adjust across proteins.
sequence1 <- seq(from = 3, to = 3600 * 6, by = 6)
comparison.bottom.vs.top <- region_contrast_results[sequence1, ]
comparison.bottom.vs.top$p.value.adjusted <- p.adjust(comparison.bottom.vs.top$p.value, method = "BH")
row.names(comparison.bottom.vs.top) <- rownames(data7)


# Extract central versus peri results and adjust across proteins.
sequence1 <- seq(from = 4, to = 3600 * 6, by = 6)
comparison.central.vs.peri <- region_contrast_results[sequence1, ]
comparison.central.vs.peri$p.value.adjusted <- p.adjust(comparison.central.vs.peri$p.value, method = "BH")
row.names(comparison.central.vs.peri) <- rownames(data7)


# Extract central versus top results and adjust across proteins.
sequence1 <- seq(from = 5, to = 3600 * 6, by = 6)
comparison.central.vs.top <- region_contrast_results[sequence1, ]
comparison.central.vs.top$p.value.adjusted <- p.adjust(comparison.central.vs.top$p.value, method = "BH")
row.names(comparison.central.vs.top) <- rownames(data7)


# Extract peri versus top results and adjust across proteins.
sequence1 <- seq(from = 6, to = 3600 * 6, by = 6)
comparison.peri.vs.top <- region_contrast_results[sequence1, ]
comparison.peri.vs.top$p.value.adjusted <- p.adjust(comparison.peri.vs.top$p.value, method = "BH")
row.names(comparison.peri.vs.top) <- rownames(data7)


# Append protein annotations and collect the six regional comparison tables.
region.vs.region <- list(
  "bottom.vs.central" = cbind(data_imp[c(1:7)], comparison.bottom.vs.central),
  "bottom.vs.peri" = cbind(data_imp[c(1:7)], comparison.bottom.vs.peri),
  "bottom.vs.top" = cbind(data_imp[c(1:7)], comparison.bottom.vs.top),
  "central.vs.peri" = cbind(data_imp[c(1:7)], comparison.central.vs.peri),
  "central.vs.top" = cbind(data_imp[c(1:7)], comparison.central.vs.top),
  "peri.vs.top" = cbind(data_imp[c(1:7)], comparison.peri.vs.top)
)


# Collect age, sex, and regional results into eight worksheets.
# These tables retain all tested proteins; no significance cutoff is applied here.
single_factor_DEP_summary <- list(
  "age" = age.vs.young,
  "sex" = female.vs.male,
  "bot vs central" = region.vs.region$bottom.vs.central,
  "bot vs peri" = region.vs.region$bottom.vs.peri,
  "bot vs top" = region.vs.region$bottom.vs.top,
  "central vs peri" = region.vs.region$central.vs.peri,
  "central vs top" = region.vs.region$central.vs.top,
  "peri vs top" = region.vs.region$peri.vs.top
)


# Export the comparison results and adjusted p-values.
export(single_factor_DEP_summary, file = "Single_factor_DEP_summary.xlsx")

age_DE <- read_excel("Single_factor_DEP_summary.xlsx", sheet = "age")


# Age volcano plot -----------------------------------------------------

res1 <- age_DE
res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)

res2 <- res1[res1$p.value.adjusted < 0.05, ]
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category), ] # Retain non-missing matrisome subcategory annotations.

fc_cut <- 1.5 # Record the fold-change cutoff.

p1 <- ggplot(res1, aes(estimate, LogPValue))
p1 <- p1 + geom_point(shape = 1, color = "gray")
p1 <- p1 + geom_point(data = res2, aes(x = estimate, y = LogPValue, colour = estimate), shape = 19, size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")
p1 <- p1 + geom_point(data = res3, aes(x = estimate, y = LogPValue), colour = "black", shape = 1, size = 2)
p1 <- p1 + geom_text_repel(data = res3, aes(x = estimate, y = LogPValue, label = Gene.symbol), size = 4.5, arrow = arrow(length = unit(0.01, 'npc')), point.padding = unit(0.6, 'lines'), max.overlaps = 1600)

p1 <- p1 + geom_hline(yintercept = -log10(0.05), colour = "red", size = 0.5, linetype = 2)
p1 <- p1 + geom_vline(xintercept = 0.5849625, colour = "red", size = 0.5, linetype = 2)
p1 <- p1 + geom_vline(xintercept = -0.5849625, colour = "red", size = 0.5, linetype = 2)

p1 <- p1 + theme_bw(base_size = 16) + labs(x = "estimate", y = "-Log10(FDR)") + ggtitle("Age vs young vocalno")
plot(p1) # Display the age volcano plot (Figure 2b).


# Sex volcano plot -----------------------------------------------------

sex_DE <- read_excel("Single_factor_DEP_summary.xlsx", sheet = "sex")
res1 <- sex_DE
res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)

res2 <- res1[res1$p.value.adjusted < 0.05, ]
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category), ] # Retain non-missing matrisome subcategory annotations.

fc_cut <- 1.5 # Record the fold-change cutoff.

p2 <- ggplot(res1, aes(estimate, LogPValue))
p2 <- p2 + geom_point(shape = 1, color = "gray")
p2 <- p2 + geom_point(data = res2, aes(x = estimate, y = LogPValue, colour = estimate), shape = 19, size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")
p2 <- p2 + geom_point(data = res3, aes(x = estimate, y = LogPValue), colour = "black", shape = 1, size = 2)
p2 <- p2 + geom_text_repel(data = res3, aes(x = estimate, y = LogPValue, label = Gene.symbol), size = 4.5, arrow = arrow(length = unit(0.01, 'npc')), point.padding = unit(0.6, 'lines'), max.overlaps = 1600)

p2 <- p2 + geom_hline(yintercept = -log10(0.05), colour = "red", size = 0.5, linetype = 2)
p2 <- p2 + geom_vline(xintercept = 0.5849625, colour = "red", size = 0.5, linetype = 2)
p2 <- p2 + geom_vline(xintercept = -0.5849625, colour = "red", size = 0.5, linetype = 2)

p2 <- p2 + theme_bw(base_size = 16) + labs(x = "estimate", y = "-Log10(FDR)") + ggtitle("Female vs male vocalno")
plot(p2) # Display the sex volcano plot (Figure S7a).


# Export age, sex, and regional volcano plots to a multipage PDF.
# Figure references: 2b, 6b, 6c, S7a, and S13.
pdf("20260915_volcano_plots_all_comparisons.pdf", width = 12, height = 8)

plot(p1)

plot(p2)


# Map regional comparison names to result worksheets for volcano plotting.
region_comparisons <- list(
  "bottom.vs.central" = "bot vs central",
  "bottom.vs.peri" = "bot vs peri",
  "bottom.vs.top" = "bot vs top",
  "central.vs.peri" = "central vs peri",
  "central.vs.top" = "central vs top",
  "peri.vs.top" = "peri vs top"
)

for (comparison_name in names(region_comparisons)) {
  res1 <- read_excel("Single_factor_DEP_summary.xlsx", sheet = region_comparisons[[comparison_name]])
  res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)
  res2 <- res1[res1$p.value.adjusted < 0.05, ]
  res2 <- res2[abs(res2$estimate) > 0.5849625, ]
  res3 <- res2[!is.na(res2$Matrisome.Sub.category), ]

  p <- ggplot(res1, aes(estimate, LogPValue)) +
    geom_point(shape = 1, color = "gray") +
    geom_point(data = res2, aes(x = estimate, y = LogPValue, colour = estimate), shape = 19, size = 2) +
    scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0) +
    geom_point(data = res3, aes(x = estimate, y = LogPValue), colour = "black", shape = 1, size = 2) +
    geom_text_repel(
      data = res3, aes(x = estimate, y = LogPValue, label = Gene.symbol), size = 4.5,
      arrow = arrow(length = unit(0.01, 'npc')), point.padding = unit(0.6, 'lines'), max.overlaps = 1600
    ) +
    geom_hline(yintercept = -log10(0.05), colour = "red", size = 0.5, linetype = 2) +
    geom_vline(xintercept = c(-0.5849625, 0.5849625), colour = "red", size = 0.5, linetype = 2) +
    theme_bw(base_size = 16) +
    labs(x = "estimate", y = "-Log10(FDR)") +
    ggtitle(paste(comparison_name, "volcano"))
  plot(p)
}


# Close the PDF graphics device.
dev.off()


# Differentially abundant protein counts --------------------------------
# Count proteins with FDR < 0.05, without a fold-change threshold.
# Positive and negative estimates define the two directions of each contrast.
count_DEPs <- function(data) {
  total_sig <- sum(data$p.value.adjusted < 0.05)
  up_reg <- sum(data$p.value.adjusted < 0.05 & data$estimate > 0)
  down_reg <- sum(data$p.value.adjusted < 0.05 & data$estimate < 0)

  # Count the subset with non-missing matrisome subcategory annotations.
  matrisome_sig <- sum(data$p.value.adjusted < 0.05 & !is.na(data$Matrisome.Sub.category))
  matrisome_up <- sum(data$p.value.adjusted < 0.05 & data$estimate > 0 & !is.na(data$Matrisome.Sub.category))
  matrisome_down <- sum(data$p.value.adjusted < 0.05 & data$estimate < 0 & !is.na(data$Matrisome.Sub.category))

  return(c(
    total_sig = total_sig,
    up_regulated = up_reg,
    down_regulated = down_reg,
    matrisome_total = matrisome_sig,
    matrisome_up = matrisome_up,
    matrisome_down = matrisome_down
  ))
}


# Initialize the table of counts across comparisons.
dep_counts <- data.frame()


# Add counts for the age comparison.
age_counts <- count_DEPs(age_DE)
dep_counts <- rbind(
  dep_counts,
  data.frame(
    comparison = "Age vs Young",
    t(age_counts)
  )
)


# Add counts for the sex comparison.
sex_counts <- count_DEPs(sex_DE)
dep_counts <- rbind(
  dep_counts,
  data.frame(
    comparison = "Female vs Male",
    t(sex_counts)
  )
)


# Add counts for regional comparisons from the specified workbook.
for (comparison_name in names(region_comparisons)) {
  data <- read_excel("DE summary1.xlsx", sheet = region_comparisons[[comparison_name]])
  counts <- count_DEPs(data)
  dep_counts <- rbind(
    dep_counts,
    data.frame(
      comparison = comparison_name,
      t(counts)
    )
  )
}


# Display the count summary.
print(dep_counts)


# Count proteins meeting both FDR and absolute log2 fold-change thresholds.
count_DEPs_FC <- function(data, fc_threshold = 0.5849625) {
  # Require FDR < 0.05 and an absolute estimate above the log2 fold-change threshold.
  total_sig <- sum(data$p.value.adjusted < 0.05 & abs(data$estimate) > fc_threshold)
  up_reg <- sum(data$p.value.adjusted < 0.05 & data$estimate > fc_threshold)
  down_reg <- sum(data$p.value.adjusted < 0.05 & data$estimate < -fc_threshold)

  # Count the subset with non-missing matrisome subcategory annotations.
  matrisome_sig <- sum(data$p.value.adjusted < 0.05 &
    abs(data$estimate) > fc_threshold &
    !is.na(data$Matrisome.Sub.category))
  matrisome_up <- sum(data$p.value.adjusted < 0.05 &
    data$estimate > fc_threshold &
    !is.na(data$Matrisome.Sub.category))
  matrisome_down <- sum(data$p.value.adjusted < 0.05 &
    data$estimate < -fc_threshold &
    !is.na(data$Matrisome.Sub.category))

  return(c(
    total_sig = total_sig,
    up_regulated = up_reg,
    down_regulated = down_reg,
    matrisome_total = matrisome_sig,
    matrisome_up = matrisome_up,
    matrisome_down = matrisome_down
  ))
}


# Initialize the table of counts using both significance and effect-size criteria.
dep_counts_FC <- data.frame()


# Set the absolute log2 fold-change threshold to approximately log2(1.5).
fc_cutoff <- 0.5849625 # log2(1.5)


# Add counts for the age comparison.
age_counts <- count_DEPs_FC(age_DE, fc_cutoff)
dep_counts_FC <- rbind(
  dep_counts_FC,
  data.frame(
    comparison = "Age vs Young",
    t(age_counts)
  )
)


# Add counts for the sex comparison.
sex_counts <- count_DEPs_FC(sex_DE, fc_cutoff)
dep_counts_FC <- rbind(
  dep_counts_FC,
  data.frame(
    comparison = "Female vs Male",
    t(sex_counts)
  )
)


# Add counts for regional comparisons from the specified workbook.
for (comparison_name in names(region_comparisons)) {
  data <- read_excel("DE summary1.xlsx", sheet = region_comparisons[[comparison_name]])
  counts <- count_DEPs_FC(data, fc_cutoff)
  dep_counts_FC <- rbind(
    dep_counts_FC,
    data.frame(
      comparison = comparison_name,
      t(counts)
    )
  )
}


# Display the count summary.
print(dep_counts_FC)


# Plot differential-abundance counts ------------------------------------
# Reshape counts into positive and negative contrast directions for plotting.
plot_data <- data.frame(
  comparison = rep(dep_counts_FC$comparison, 2),
  regulation = c(
    rep("Up-regulated", nrow(dep_counts_FC)),
    rep("Down-regulated", nrow(dep_counts_FC))
  ),
  count = c(dep_counts_FC$up_regulated, dep_counts_FC$down_regulated)
)


# Set the display order of the eight comparisons.
plot_data$comparison <- factor(plot_data$comparison,
  levels = c(
    "Age vs Young", "Female vs Male",
    "bottom.vs.central", "bottom.vs.peri",
    "bottom.vs.top", "central.vs.peri",
    "central.vs.top", "peri.vs.top"
  )
)


# Plot protein counts by comparison and direction.
p_dep <- ggplot(plot_data, aes(x = comparison, y = count, fill = regulation)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 120, by = 20)) + # Use a fixed y-axis range.
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Comparison",
    y = "Number of differentially expressed proteins",
    title = "DEPs in different comparisons (FDR < 0.05, |FC| > 1.5)"
  ) +
  geom_text(aes(label = count),
    position = position_dodge(width = 0.7),
    vjust = -0.5
  )


# Export the protein-count bar plot (Figure S12a).
pdf("20260915_DEP_numbers_barplot.pdf", width = 10, height = 6)
print(p_dep)
dev.off()


# Prepare the corresponding counts for matrisome proteins.
plot_data_matrisome <- data.frame(
  comparison = rep(dep_counts_FC$comparison, 2),
  regulation = c(
    rep("Up-regulated", nrow(dep_counts_FC)),
    rep("Down-regulated", nrow(dep_counts_FC))
  ),
  count = c(dep_counts_FC$matrisome_up, dep_counts_FC$matrisome_down)
)

plot_data_matrisome$comparison <- factor(plot_data_matrisome$comparison,
  levels = c(
    "Age vs Young", "Female vs Male",
    "bottom.vs.central", "bottom.vs.peri",
    "bottom.vs.top", "central.vs.peri",
    "central.vs.top", "peri.vs.top"
  )
)


# Plot matrisome protein counts by comparison and direction.
p_dep_matrisome <- ggplot(plot_data_matrisome, aes(x = comparison, y = count, fill = regulation)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 120, by = 20)) + # Use a fixed y-axis range.
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Comparison",
    y = "Number of differentially expressed matrisome proteins",
    title = "Matrisome DEPs in different comparisons (FDR < 0.05, |FC| > 1.5)"
  ) +
  geom_text(aes(label = count),
    position = position_dodge(width = 0.7),
    vjust = -0.5
  )


# Export the matrisome-count bar plot (Figure S12b).
pdf("20260915_matrisome_DEP_numbers_barplot.pdf", width = 10, height = 6)
print(p_dep_matrisome)
dev.off()


# Matrisome heatmap and effect-size bar plot -----------------------------
# Select age-associated matrisome proteins with FDR < 0.05
# and absolute log2 fold change > 0.5849625.
# Define the experimental-group order.
# The actual heatmap column order is read from Sample_index.xlsx below.
group_order <- c("ABF", "ACF", "APF", "ATF", "ABM", "ACM", "APM", "ATM", "YBF", "YCF", "YPF", "YTF", "YBM", "YCM", "YPM", "YTM")

matrisome_list <- c()

res1 <- age_DE
res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)
res2 <- res1[res1$p.value.adjusted < 0.05, ]

res2 <- res2[(res2$estimate > 0.5849625) | (res2$estimate < (-1) * 0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category), ]
matrisome_list <- c(res3$Uniprot.ID, matrisome_list) # Collect selected protein identifiers.

sample_group <- read_excel("Sample_index.xlsx")
sample_group_order <- unlist(sample_group$Index)

hm_df_new <- data7[match(matrisome_list, rownames(data7)), ]
hm_df_new <- hm_df_new[, sample_group_order]


# Scale each protein row and cluster rows by correlation distance.
# Keep sample columns in the supplied order.
p8 <- pheatmap(
  hm_df_new,
  rev(brewer.pal(n = 11, name = "RdBu")),
  scale = c('row'),
  gaps_col = NULL, fontsize_row = 8,
  show_rownames = F, show_colnames = T, cluster_cols = F, cluster_rows = T,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation", border_color = "grey95",
  cutree_rows = 2, # Cut the row dendrogram into two clusters.
  main = ""
)


# Extract the clustered protein order and match it to the age results.
gene_order <- p8$tree_row$order
ordered_genes <- rownames(hm_df_new)[gene_order]
ordered_genes_table <- age_DE[match(ordered_genes, age_DE$Uniprot.ID), ]
mean_table <- ordered_genes_table[c(3, 9)]
mean_table$type <- ifelse(mean_table$estimate > 0, "Positive", "Negative")


# Reverse the clustered order for the horizontal effect-size bar plot.
# Preserve this order in the gene-symbol factor levels.

mean_table <- mean_table[rev(seq_len(nrow(mean_table))), ]
mean_table$Gene.symbol <- factor(mean_table$Gene.symbol, levels = mean_table$Gene.symbol)


# Plot model estimates on the log2 intensity scale, colored by sign.
p9 <- ggplot(mean_table, aes(x = Gene.symbol, y = estimate, fill = type)) +
  geom_col() +
  scale_fill_manual(
    name = "Bar plot",
    values = c("Positive" = "#BC5DA3", "Negative" = "#3AAA69")
  ) +
  # Rotate axis labels for readability.
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p9 <- p9 + ylim(-2, 2) + coord_flip()


# Display the effect-size bar plot.
plot(p9)


# Export the heatmap and bar plot as separate PDF pages (Figure 2c).
pdf("20260915_matrisome_DEP_heatmap_and_barplot.pdf", width = 10, height = 6)
print(p8)
print(p9)
dev.off()


# Save the current workspace, including intermediate results.
save.image("Single_factor_DEP_summary_data_summary.RData")
