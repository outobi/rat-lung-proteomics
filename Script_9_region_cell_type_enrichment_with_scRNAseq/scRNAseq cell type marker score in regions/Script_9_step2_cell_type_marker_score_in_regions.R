# =============================================================================
# Script 9 (scRNAseq cell type marker score in regions), step 2:
#   Cell-type marker scores in the proteomics data, compared between regions
# =============================================================================
#
# Purpose:
#   Estimate how strongly each lung cell type is represented in each lung
#   region of the proteomics data, using scRNA-seq marker genes (step 1) as
#   cell-type signatures, and test for regional differences.
#
# Workflow:
#   1. Marker selection (from step 1): keep cell-type markers with
#      p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 >= 0.25 and
#      pct.1 - pct.2 >= 0.15 whose gene is among the 3,600 quantified
#      proteins, then the top 30 per cell type by avg_log2FC
#      (21 cell types x 30 markers).
#   2. Score: z-score each protein across the 96 samples; a cell type's score
#      in a sample is the mean z-score of its markers.
#   3. For each cell type, fit Score ~ Age + Sex + Region + (1 | PID) and
#      obtain (a) the Region omnibus F-test (BH across the 21 cell types),
#      (b) estimated regional means (emmeans) for the heatmap, and
#      (c) pairwise region contrasts (unadjusted; BH across cell types within
#      each contrast).
#   4. Airway epithelial (the only cell type with Region FDR < 0.05 in the
#      reported run, FDR = 3.2e-4): Tukey-adjusted pairwise comparisons and a
#      boxplot of the per-sample scores.
#
# Inputs (in this folder):
#   - Raredon_rat_celltype_FindAllMarkers.csv  (output of step 1)
#   - expression matrix.xlsx
#       Imputed, normalized log2 protein intensities (identical to sheet
#       "data_imp" of Script_3's protein intensity workbook): 7 annotation
#       columns + S1-S96; 3,600 proteins.
#   - sample metadata.xlsx
#       Sample.ID, Sample.group, PID, Age, Region, Sex (identical to
#       Script_3's Sample_metadata.xlsx).
#
# Outputs (written to this folder):
#   - Raredon_rat_celltype_top30_proteomics_markers.csv
#       The 630 selected markers (top 30 per cell type).
#   - celltype_region_LMM_FDR.csv
#       Region omnibus test per cell type (F, df, p, BH FDR, Significant).
#   - celltype_region_emmeans.csv
#       Estimated regional mean score per cell type (emmean, SE, df, CI).
#   - heatmap summary.pdf
#       Heatmap of the regional emmeans (cell types with Region FDR < 0.05
#       marked with "*").
#   - airway epithelial cells score.pdf
#       Airway epithelial score per sample by region (lines join the 4
#       regions of each animal), with Tukey-adjusted significance marks.
#   - Script_9_step2_cell_type_marker_score_in_regions.RData
#       Full R workspace.
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# =============================================================================

# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages.
# readxl: read_excel(); dplyr + tidyr: table handling;
# lme4 + lmerTest: lmer() / anova(); emmeans: emmeans() / contrast() / pairs();
# ggplot2 + ggpubr: airway boxplot with significance marks;
# pheatmap: heatmap.
library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(pheatmap)


# Step 1: Read inputs ----------------------------------------------------------

meta <- read_excel("sample metadata.xlsx")
expr <- as.data.frame(read_excel("expression matrix.xlsx"))
rownames(expr) <- expr$Uniprot.ID

celltype_markers <- read.csv("Raredon_rat_celltype_FindAllMarkers.csv", check.names = FALSE)

# Keep the cell-type order of the Seurat output (factor levels in order of
# first appearance), as when step 1 and step 2 ran in one session.
celltype_markers$cluster <- factor(celltype_markers$cluster, levels = unique(celltype_markers$cluster))


# Step 2: Select cell-type markers ---------------------------------------------

# Specific scRNA-seq markers.
celltype_markers_specific <- celltype_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.5,
    pct.1 >= 0.25,
    (pct.1 - pct.2) >= 0.15
  )
table(celltype_markers_specific$cluster)

# Keep markers that are among the quantified proteins.
celltype_markers_protein <- celltype_markers_specific %>%
  filter(gene %in% expr$Gene.symbol)
table(celltype_markers_protein$cluster)
# Reported run:
#   VCAM1+ endothelial 136, Alveolar macrophage 471, ATI 328,
#   Vascular endothelial 87, T cell 66, NK cell 78, Monocyte 146, B cell 34,
#   ATII 167, Mast cell 45, Airway epithelial 362, Dendritic cell 231,
#   Sox9+ epithelial 293, Neutrophil 47, Mrc1+ macrophage 45,
#   Interstitial macrophage 131, Capillary endothelial 153, Mesothelial 483,
#   Fibroblast 200, Smooth muscle 99, pDC/DC 279

# Top 30 per cell type by avg_log2FC.
celltype_marker_top30 <- celltype_markers_protein %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30, with_ties = FALSE) %>%
  ungroup()

write.csv(celltype_marker_top30, "Raredon_rat_celltype_top30_proteomics_markers.csv", row.names = FALSE)


# Step 3: Cell-type scores per sample ------------------------------------------

# Proteins with a gene symbol ("-" = no symbol); rows = gene symbols,
# columns = S1-S96 (the first 7 columns are annotations).
expr1 <- expr[expr$Gene.symbol != "-", ]
expr_mat <- as.data.frame(expr1[, -c(1:7)])
rownames(expr_mat) <- expr1$Gene.symbol

# Gene-wise z-score across the 96 samples.
expr_z <- t(scale(t(expr_mat)))

# Score = mean z-score of the markers of each cell type.
marker_list <- split(celltype_marker_top30$gene, celltype_marker_top30$cluster)

celltype_score <- sapply(marker_list, function(markers) {
  genes <- intersect(markers, rownames(expr_z))
  colMeans(expr_z[genes, , drop = FALSE], na.rm = TRUE)
})

celltype_score <- as.data.frame(celltype_score)
celltype_score$Sample.ID <- rownames(celltype_score)
celltype_score <- merge(meta, celltype_score, by = "Sample.ID")

celltype_cols <- names(marker_list) # 21 cell types


# Step 4: Mixed model per cell type --------------------------------------------

# One fit per cell type; the omnibus test, regional means and pairwise
# contrasts all come from this fit.
celltype_region_omnibus <- list()
emm_list <- list()
celltype_lmm_results <- list()

for (ct in celltype_cols) {
  dat <- celltype_score %>%
    dplyr::select(PID, Age, Sex, Region, all_of(ct)) %>%
    dplyr::rename(Score = all_of(ct))

  fit <- lmer(Score ~ Age + Sex + Region + (1 | PID), data = dat)

  # (a) Region omnibus F-test.
  aov_tab <- anova(fit)
  celltype_region_omnibus[[ct]] <- data.frame(
    Cell.type = ct,
    F.value = aov_tab["Region", "F value"],
    df1 = aov_tab["Region", "NumDF"],
    df2 = aov_tab["Region", "DenDF"],
    p.value = aov_tab["Region", "Pr(>F)"]
  )

  # (b) Estimated regional means, adjusted for age and sex.
  emm <- emmeans(fit, ~ Region)
  emm_df <- as.data.frame(emm)
  emm_df$Cell.type <- ct
  emm_list[[ct]] <- emm_df

  # (c) Pairwise region contrasts (unadjusted).
  tmp <- as.data.frame(contrast(emm, method = "pairwise", adjust = "none"))
  tmp$Cell.type <- ct
  celltype_lmm_results[[ct]] <- tmp
}

# Omnibus test: BH across the 21 cell types.
celltype_region_omnibus <- bind_rows(celltype_region_omnibus) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  arrange(FDR) %>%
  mutate(Significant = FDR < 0.05)

celltype_region_omnibus %>% filter(FDR < 0.05) # Airway epithelial only

emm_region <- bind_rows(emm_list)

# Pairwise contrasts: BH across cell types within each contrast (written
# with ave() rather than group_by() + mutate(), which plyr can silently
# override; see Script_8).
celltype_region_results <- bind_rows(celltype_lmm_results)
celltype_region_results$FDR <- ave(
  celltype_region_results$p.value,
  celltype_region_results$contrast,
  FUN = function(p) p.adjust(p, method = "BH")
)

write.csv(emm_region, "celltype_region_emmeans.csv", row.names = FALSE)
write.csv(celltype_region_omnibus, "celltype_region_LMM_FDR.csv", row.names = FALSE)


# Step 5: Airway epithelial: post-hoc comparisons and boxplot ------------------

fit_airway <- lmer(`Airway epithelial` ~ Age + Sex + Region + (1 | PID), data = celltype_score)
emm_airway <- emmeans(fit_airway, ~ Region)

# Tukey-adjusted pairwise region comparisons.
pair_airway_df <- as.data.frame(pairs(emm_airway, adjust = "tukey")) %>%
  mutate(
    p.signif = case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
pair_airway_df

celltype_score$Region <- factor(celltype_score$Region, levels = c("Bottom", "Central", "Peripheral", "Top"))

p_airway <- ggplot(celltype_score, aes(x = Region, y = `Airway epithelial`)) +
  geom_line(aes(group = PID), alpha = 0.20, linewidth = 0.5) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.35) +
  geom_point(size = 2, alpha = 0.75) +
  theme_bw() +
  labs(x = "Lung region", y = "Airway epithelial marker score") +
  theme(
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

# Significance brackets for the Tukey comparisons with p < 0.05.
stat_airway_sig <- pair_airway_df %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".* - ", "", contrast)
  ) %>%
  filter(p.value < 0.05)
stat_airway_sig$y.position <- seq(1.65, by = 0.18, length.out = nrow(stat_airway_sig))

p_airway <- p_airway +
  stat_pvalue_manual(
    stat_airway_sig,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  ) +
  coord_cartesian(ylim = c(-1.2, 2.5), clip = "off")

p_airway
ggsave("airway epithelial cells score.pdf", p_airway, width = 11, height = 8.5)


# Step 6: Heatmap of regional means --------------------------------------------

heatmap_mat <- emm_region %>%
  dplyr::select(Cell.type, Region, emmean) %>%
  pivot_wider(names_from = Region, values_from = emmean) %>%
  as.data.frame()
rownames(heatmap_mat) <- heatmap_mat$Cell.type
heatmap_mat$Cell.type <- NULL

# Cell types grouped by compartment.
cell_order <- c(
  "ATI", "ATII", "Airway epithelial", "Sox9+ epithelial",
  "Fibroblast", "Smooth muscle", "Mesothelial",
  "Capillary endothelial", "Vascular endothelial", "VCAM1+ endothelial",
  "Alveolar macrophage", "Interstitial macrophage", "Mrc1+ macrophage",
  "Monocyte", "Dendritic cell", "pDC/DC", "Neutrophil", "Mast cell",
  "T cell", "B cell", "NK cell"
)
heatmap_mat <- heatmap_mat[cell_order, ]

# Mark cell types with Region FDR < 0.05.
fdr_vec <- celltype_region_omnibus$FDR[match(rownames(heatmap_mat), celltype_region_omnibus$Cell.type)]
rownames(heatmap_mat) <- ifelse(fdr_vec < 0.05, paste0(rownames(heatmap_mat), " *"), rownames(heatmap_mat))

p_heatmap <- pheatmap(
  as.matrix(heatmap_mat),
  scale = "none", # emmeans are already on the z-score scale
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA,
  angle_col = 45
)

pdf("heatmap summary.pdf", width = 11, height = 8.5)
print(p_heatmap)
dev.off()


# Save the current workspace.
save.image("Script_9_step2_cell_type_marker_score_in_regions.RData")
