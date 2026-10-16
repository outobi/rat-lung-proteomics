# =============================================================================
# Script 9 (region protein marker in scRNAseq cell type), step 2:
#   Region-specific proteins and their expression across scRNA-seq cell types
# =============================================================================
#
# Purpose:
#   Find proteins that are higher in one lung region than in the other three
#   (within each age group) and ask in which lung cell types the
#   corresponding genes are expressed, using the annotated rat lung scRNA-seq
#   reference from step 1.
#
# Workflow:
#   1. Region-specific proteins: within Young (peripubertal) and within Aged
#      (young adult) samples separately, compare each region with the other
#      three regions for each of the 3,600 proteins (two-sided Wilcoxon
#      rank-sum test, unpaired; BH across proteins). Log2FC = mean log2
#      intensity in the region minus mean in the other regions.
#      Region-specific = FDR < 0.05 and Log2FC > 0.25.
#   2. Keep proteins with a gene symbol that is present in the scRNA-seq
#      gene list. Reported run (FDR/Log2FC -> with symbol -> in scRNA-seq):
#      Young Central 113 -> 112 -> 105; Young Peripheral 103 -> 103 -> 92;
#      Aged Bottom 1 -> 1 -> 1; Aged Central 19 -> 19 -> 17;
#      Aged Peripheral 3 -> 3 -> 3; Aged Top 1 -> 1 -> 1;
#      Young Bottom and Young Top: none.
#   3. Signature score per cell type: log-normalize the scRNA-seq counts
#      (all 10,953 cells), z-score each signature gene across all cells, and
#      average the z-scores over the signature genes and the cells of each
#      cell type (21 cell types; and 5 compartments: Epithelial,
#      Endothelial, Stromal, Myeloid, Lymphoid). Positive = signature genes
#      enriched in that cell type, negative = depleted.
#
# Correction (2026-10): in the original script the per-cell-type scores were
#   sorted from high to low within each signature before the cell-type labels
#   were attached, so every signature's highest score was labelled "ATI", the
#   second "ATII", and so on. Scores now stay attached to their cell type;
#   the score values themselves are unchanged. The 5-compartment scores were
#   not affected.
#
# Inputs (in this folder):
#   - expression matrix.xlsx
#       Imputed, normalized log2 protein intensities (identical to sheet
#       "data_imp" of Script_3's protein intensity workbook): 7 annotation
#       columns + S1-S96; 3,600 proteins.
#   - sample metadata.xlsx
#       Sample.ID, Sample.group, PID, Age ("Young" / "Aged"), Region, Sex.
#   - Raredon_rat_cell_annotation.csv  (output of step 1)
#   - GSM3926541_mRat.* and GSM3926542_fRat.* raw count files (GEO GSE133747)
#
# Outputs (written to this folder):
#   - region_specific_protein_summary.xlsx
#       8 sheets (Young / Aged x Bottom / Central / Peripheral / Top vs all):
#       p.value, FDR, Log2FC for all 3,600 proteins.
#   - region_specific_upregulated_proteins.xlsx
#       8 sheets: region-specific proteins used as signatures.
#   - cell_type_signature_scores.xlsx
#       Signature scores: "cell types" (21 x 6) and "compartments" (5 x 6).
#   - p1 enrichment.pdf / p2 depletion.pdf
#       Young Central and Young Peripheral signatures: cell types with a
#       positive (enrichment) or negative (depletion) score.
#   - p1 enrichment bubble plot full.pdf / p2 depletion bubble plot full.pdf
#       Same for all 6 signatures.
#   - compartment enrichment.pdf / compartment depletion.pdf
#       5-compartment scores for all 6 signatures.
#   - Script_9_step2_region_protein_signature_in_cell_types.RData
#       Result tables, scores and plot objects (not the large count matrix).
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# =============================================================================

# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages.
# readxl: read_excel(); rio: export(); dplyr: filter / pull;
# Matrix + Seurat: count matrices and log-normalization;
# reshape2: melt(); ggplot2 + patchwork: bubble plots.
library(readxl)
library(rio)
library(dplyr)
library(Matrix)
library(Seurat)
library(reshape2)
library(ggplot2)
library(patchwork)


# Step 1: Region-specific proteins ---------------------------------------------

meta <- read_excel("sample metadata.xlsx")
expr <- as.data.frame(read_excel("expression matrix.xlsx"))
rownames(expr) <- expr$Uniprot.ID

# Young and Aged samples; the column order must match the metadata order.
meta_young <- meta %>% dplyr::filter(Age == "Young")
meta_aged <- meta %>% dplyr::filter(Age == "Aged")
expr_young <- expr[, colnames(expr) %in% meta_young$Sample.ID]
expr_aged <- expr[, colnames(expr) %in% meta_aged$Sample.ID]
stopifnot(identical(colnames(expr_young), meta_young$Sample.ID),
          identical(colnames(expr_aged), meta_aged$Sample.ID))

# One region vs the other three regions, Wilcoxon rank-sum test per protein.
wilcox_region_vs_rest <- function(target_region, expr_mat, metadata) {
  wilcox_p <- numeric(nrow(expr_mat))
  logFC <- numeric(nrow(expr_mat))
  for (i in seq_len(nrow(expr_mat))) {
    x <- as.numeric(expr_mat[i, metadata$Region == target_region])
    y <- as.numeric(expr_mat[i, metadata$Region != target_region])
    wilcox_p[i] <- wilcox.test(x = x, y = y)$p.value
    logFC[i] <- mean(x) - mean(y)
  }
  data.frame(
    Uniprot.ID = rownames(expr_mat),
    p.value = wilcox_p,
    FDR = p.adjust(wilcox_p, method = "BH"),
    Log2FC = logFC
  )
}

regions <- c("Bottom", "Central", "Peripheral", "Top")
region_specific_feature <- c(
  setNames(lapply(regions, wilcox_region_vs_rest, expr_mat = expr_young, metadata = meta_young),
           paste0("Young_", regions, "_vs_all")),
  setNames(lapply(regions, wilcox_region_vs_rest, expr_mat = expr_aged, metadata = meta_aged),
           paste0("Aged_", regions, "_vs_all"))
)

export(region_specific_feature, "region_specific_protein_summary.xlsx")


# Step 2: Region-specific protein signatures -----------------------------------

log2fc_cutoff <- 0.25
fdr_cutoff <- 0.05

# scRNA-seq gene list (identical for the male and female files).
sc_genes <- read.delim("GSM3926542_fRat.genes.tsv.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)$x

# Region-specific proteins -> keep those with a gene symbol ("-" = none)
# that is present in the scRNA-seq gene list.
select_signature <- function(res) {
  rsp <- res[res$Log2FC > log2fc_cutoff & res$FDR < fdr_cutoff, ]
  rsp <- cbind("Gene.symbol" = expr$Gene.symbol[match(rsp$Uniprot.ID, expr$Uniprot.ID)], rsp)
  rsp <- rsp[rsp$Gene.symbol != "-", ]
  rsp[rsp$Gene.symbol %in% sc_genes, ]
}

region_specific_proteins <- lapply(region_specific_feature, select_signature)
names(region_specific_proteins) <- sub("_vs_all", "", names(region_specific_feature))
sapply(region_specific_proteins, nrow)
# Reported run: Young_Bottom 0, Young_Central 105, Young_Peripheral 92,
#   Young_Top 0, Aged_Bottom 1, Aged_Central 17, Aged_Peripheral 3, Aged_Top 1

export(region_specific_proteins, "region_specific_upregulated_proteins.xlsx")

# Signatures used below (Young Bottom and Young Top have no proteins).
signatures <- list(
  "Young Central"    = region_specific_proteins$Young_Central$Gene.symbol,
  "Young Peripheral" = region_specific_proteins$Young_Peripheral$Gene.symbol,
  "Aged Bottom"      = region_specific_proteins$Aged_Bottom$Gene.symbol,
  "Aged Central"     = region_specific_proteins$Aged_Central$Gene.symbol,
  "Aged Peripheral"  = region_specific_proteins$Aged_Peripheral$Gene.symbol,
  "Aged Top"         = region_specific_proteins$Aged_Top$Gene.symbol
)


# Step 3: Log-normalized scRNA-seq expression with cell-type labels ------------

m_counts <- readMM("GSM3926541_mRat.raw.data.mtx.gz")
f_counts <- readMM("GSM3926542_fRat.raw.data.mtx.gz")
gene_names <- make.unique(sc_genes)
rownames(m_counts) <- gene_names
colnames(m_counts) <- read.delim("GSM3926541_mRat.barcodes.tsv.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)$x
rownames(f_counts) <- gene_names
colnames(f_counts) <- read.delim("GSM3926542_fRat.barcodes.tsv.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)$x

# Log-normalize the combined counts (scale factor 10,000).
query_seurat <- CreateSeuratObject(counts = cbind(m_counts, f_counts))
query_seurat <- NormalizeData(query_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
normalized_count <- GetAssayData(query_seurat, assay = "RNA", layer = "data")
rm(m_counts, f_counts, query_seurat)

# Cell-type label of each cell, matched by barcode.
cell_annotation <- read.csv("Raredon_rat_cell_annotation.csv")
cell_type <- cell_annotation$cell_type[match(colnames(normalized_count), cell_annotation$Cell)]
stopifnot(!anyNA(cell_type))


# Step 4: Signature scores per cell type and compartment -----------------------

# The 21 cell types, grouped by compartment (Cycling and Mixed clusters are
# not scored but are included in the z-scores).
celltype_group <- list(
  Epithelial = c("ATI", "ATII", "Airway epithelial", "Sox9+ epithelial"),
  Stromal = c("Fibroblast", "Smooth muscle", "Mesothelial"),
  Endothelial = c("Capillary endothelial", "Vascular endothelial", "VCAM1+ endothelial"),
  Myeloid = c("Alveolar macrophage", "Interstitial macrophage", "Mrc1+ macrophage", "Monocyte",
              "Dendritic cell", "pDC/DC", "Neutrophil", "Mast cell"),
  Lymphoid = c("T cell", "B cell", "NK cell")
)
query_cell_types <- unlist(celltype_group, use.names = FALSE)

# Gene-wise z-score of the signature genes across all cells.
signature_z <- function(signature_genes) {
  test <- as.matrix(normalized_count[signature_genes, , drop = FALSE])
  t(scale(t(test)))
}

# Mean z-score over the signature genes and the cells of each group.
# Scores are returned in the order of `groups`, named by group (the original
# version sorted the scores, which detached them from their labels).
score_groups <- function(signature_genes, groups) {
  z_test <- signature_z(signature_genes)
  sapply(groups, function(g) mean(z_test[, cell_type %in% g, drop = FALSE]))
}

# 21 cell types x 6 signatures.
query.summary <- as.data.frame(sapply(signatures, score_groups, groups = setNames(as.list(query_cell_types), query_cell_types)))

# 5 compartments x 6 signatures.
large_query_summary <- as.data.frame(sapply(signatures, score_groups, groups = celltype_group))

# Cells per compartment (reported run: Epithelial 3,333, Stromal 702,
# Endothelial 1,650, Myeloid 3,147, Lymphoid 1,823).
sapply(celltype_group, function(x) sum(cell_type %in% x))

export(
  list("cell types" = cbind(Cell.type = rownames(query.summary), query.summary),
       "compartments" = cbind(Compartment = rownames(large_query_summary), large_query_summary)),
  "cell_type_signature_scores.xlsx"
)


# Step 5: Bubble plots ---------------------------------------------------------
# Rows are grouped by compartment (colored band on the left). Enrichment
# plots show positive scores (size = score), depletion plots negative scores
# (size = |score|); other cells are left empty.

compartment_colors <- c(Epithelial = "#9EF5D3", Stromal = "#F5F59D", Endothelial = "#F2A7C7",
                        Myeloid = "#94B3F5", Lymphoid = "#BDE3FA")
row_compartment <- setNames(rep(names(celltype_group), lengths(celltype_group)), query_cell_types)

bubble_plot <- function(score_table, sig_names, direction, size_max, row_groups) {
  d <- score_table[, sig_names, drop = FALSE]
  d$Row <- rownames(d)
  d <- melt(d, id.vars = "Row", variable.name = "Signature", value.name = "Score")
  d$Row <- factor(d$Row, levels = rev(names(row_groups)))
  d$Size <- if (direction == "enrichment") ifelse(d$Score > 0, d$Score, NA) else ifelse(d$Score < 0, -d$Score, NA)

  bubbles <- ggplot(d, aes(x = Signature, y = Row)) +
    geom_point(aes(size = Size, fill = Signature), shape = 21, na.rm = TRUE) +
    scale_size_continuous(limits = c(0, size_max), range = c(1, 10),
                          name = if (direction == "enrichment") "Enrichment score" else "Depletion score") +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme(
      legend.key = element_blank(),
      axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
      axis.ticks.y = element_blank()
    )

  # Colored compartment band with labels, aligned with the bubble rows.
  band <- data.frame(Row = factor(names(row_groups), levels = rev(names(row_groups))), Group = row_groups)
  band_labels <- aggregate(as.numeric(Row) ~ Group, data = band, FUN = mean)
  names(band_labels)[2] <- "y"
  side <- ggplot(band, aes(x = 1, y = Row, fill = Group)) +
    geom_tile(width = 1, height = 1) +
    geom_text(data = band_labels, aes(x = 1, y = y, label = Group), inherit.aes = FALSE,
              angle = 90, fontface = "bold", size = 4) +
    scale_fill_manual(values = compartment_colors, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void()

  side + bubbles + plot_layout(widths = c(0.08, 1))
}

young <- c("Young Central", "Young Peripheral")
all_sigs <- names(signatures)

# Young signatures (size scale 0-0.2).
p11 <- bubble_plot(query.summary, young, "enrichment", 0.2, row_compartment)
p22 <- bubble_plot(query.summary, young, "depletion", 0.2, row_compartment)

# All 6 signatures.
p1 <- bubble_plot(query.summary, all_sigs, "enrichment", max(query.summary), row_compartment)
p2 <- bubble_plot(query.summary, all_sigs, "depletion", max(-query.summary), row_compartment)

# 5 compartments, all 6 signatures.
compartments <- setNames(names(celltype_group), names(celltype_group))
p_large_enrich <- bubble_plot(large_query_summary, all_sigs, "enrichment", max(large_query_summary), compartments)
p_large_deplete <- bubble_plot(large_query_summary, all_sigs, "depletion", max(-large_query_summary), compartments)

p11
p22

ggsave("p1 enrichment.pdf", p11, width = 11, height = 8.5)
ggsave("p2 depletion.pdf", p22, width = 11, height = 8.5)
ggsave("p1 enrichment bubble plot full.pdf", p1, width = 11, height = 8.5)
ggsave("p2 depletion bubble plot full.pdf", p2, width = 11, height = 8.5)
ggsave("compartment enrichment.pdf", p_large_enrich, width = 11, height = 6)
ggsave("compartment depletion.pdf", p_large_deplete, width = 11, height = 6)


# Step 6: Save results ---------------------------------------------------------

save(
  region_specific_feature, region_specific_proteins, signatures,
  query.summary, large_query_summary, cell_type,
  p11, p22, p1, p2, p_large_enrich, p_large_deplete,
  file = "Script_9_step2_region_protein_signature_in_cell_types.RData"
)
