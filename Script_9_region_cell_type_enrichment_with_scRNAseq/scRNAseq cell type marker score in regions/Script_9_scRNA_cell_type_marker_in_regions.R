setwd(dirname(rstudioapi::getSourceEditorContext()$path)) ## set working dir as where the script locates


library(Matrix)
library(Seurat)

library(readxl)
library(dplyr)
library(openxlsx)

meta <- read_excel("sample metadata.xlsx")
expr <- read_excel("expression matrix.xlsx")
expr<- as.data.frame(expr)
rownames(expr) <- expr$Uniprot.ID

# load data matrix-----------
m_counts <- readMM("GSM3926541_mRat.raw.data.mtx.gz")

m_genes <- read.delim(
  "GSM3926541_mRat.genes.tsv.gz",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)


m_barcodes <- read.delim(
  "GSM3926541_mRat.barcodes.tsv.gz",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

f_counts <- readMM("GSM3926542_fRat.raw.data.mtx.gz")

f_genes <- read.delim(
  "GSM3926542_fRat.genes.tsv.gz",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

f_barcodes <- read.delim(
  "GSM3926542_fRat.barcodes.tsv.gz",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)


dim(m_counts)
dim(f_counts)


gene_names <- make.unique(m_genes$x)

rownames(m_counts) <- gene_names
colnames(m_counts) <- m_barcodes$x

rownames(f_counts) <- gene_names
colnames(f_counts) <- f_barcodes$x


# build seurat objects------
m_rat <- CreateSeuratObject(
  counts = m_counts,
  project = "mRat",
  min.cells = 0,
  min.features = 0
)
m_rat$Sex <- "Male"

f_rat <- CreateSeuratObject(
  counts = f_counts,
  project = "fRat",
  min.cells = 0,
  min.features = 0
)
f_rat$Sex <- "Female"

# mt percentage label ----------

# mitochondrial genes
mt_genes <- grep(
  "^Mt-",
  rownames(m_rat),
  value = TRUE
)

mt_genes
length(mt_genes) #11

m_rat[["percent.mt"]] <- PercentageFeatureSet(
  m_rat,
  features = mt_genes
)

f_rat[["percent.mt"]] <- PercentageFeatureSet(
  f_rat,
  features = mt_genes
)

# 

m_rat <- NormalizeData(m_rat)
m_rat <- FindVariableFeatures(
  m_rat,
  selection.method = "vst",
  nfeatures = 2000
)

f_rat <- NormalizeData(f_rat)
f_rat <- FindVariableFeatures(
  f_rat,
  selection.method = "vst",
  nfeatures = 2000
)


m_top <- head(VariableFeatures(m_rat), 1000)
f_top <- head(VariableFeatures(f_rat), 1000)

integration_features <- intersect(m_top, f_top)

length(integration_features)

rat_list <- list(mRat = m_rat, fRat = f_rat)

anchors <- FindIntegrationAnchors(
  object.list = rat_list,
  anchor.features = integration_features,
  reduction = "cca",
  dims = 1:30
)

rat <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)



DefaultAssay(rat) <- "integrated"

rat <- ScaleData(
  rat,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)


rat <- RunPCA(
  rat,
  features = integration_features,
  npcs = 30,
  verbose = FALSE
)

ElbowPlot(rat, ndims = 30)

DimPlot(
  rat,
  reduction = "pca",
  group.by = "orig.ident"
)





rat <- FindNeighbors(
  rat,
  dims = 1:20
)

rat <- FindClusters(
  rat,
  resolution = 0.5
)

rat <- RunTSNE(
  rat,
  dims = 1:20
)

rat <- RunUMAP(
  rat,
  dims = 1:20
)


DimPlot(
  rat,
  reduction = "tsne",
  label = TRUE
)

DimPlot(
  rat,
  reduction = "umap",
  label = TRUE
)

table(Idents(rat))



DefaultAssay(rat) <- "RNA"

rat <- JoinLayers(
  rat,
  assay = "RNA"
)

Layers(rat[["RNA"]])


markers <- FindAllMarkers(
  rat,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)


top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)

top10

write.csv(
  markers,
  "Raredon_rat_FindAllMarkers.csv",
  row.names = FALSE
)

write.csv(
  top10,
  "Raredon_rat_top10_markers_per_cluster.csv",
  row.names = FALSE
)



DefaultAssay(rat) <- "RNA"

marker_panel <- c(
  # epithelial
  "Ager", "Aqp5", "Pdpn",
  "Sftpc", "Sftpb", "Slc34a2",
  "Scgb1a1", "Scgb3a1",
  "Pifo", "Cdhr3",
  "Sox9", "Pou2f3",
  
  # endothelial
  "Car4", "Apln",
  "Plvap", "Gpihbp1",
  "Vcam1", "Vwf",
  
  # mesenchymal
  "Pdgfra", "Tcf21",
  "Col13a1", "Col14a1",
  "Acta2", "Myh11",
  "Msln", "Wt1",
  
  # immune
  "Cd3d", "Cd3e",
  "Ms4a1", "Cd19",
  "Nkg7", "Ncr1",
  "Pparg", "Fabp4",
  "C1qa", "Folr2",
  "S100a8", "S100a9",
  "Cma1", "Tpsb2",
  "Flt3", "Siglech"
)

DotPlot(
  rat,
  features = marker_panel,
  group.by = "seurat_clusters"
) + RotatedAxis()


cluster_labels <- c(
  "0"  = "ATII",
  "1"  = "ATI",
  "2"  = "Vascular endothelial",
  "3"  = "Alveolar macrophage",
  "4"  = "T cell",
  "5"  = "Monocyte",
  "6"  = "Dendritic cell",
  "7"  = "B cell",
  "8"  = "NK cell",
  "9"  = "Sox9+ epithelial",
  "10" = "Fibroblast",
  "11" = "VCAM1+ endothelial",
  "12" = "Airway epithelial",
  "13" = "Smooth muscle",
  "14" = "Capillary endothelial",
  "15" = "Neutrophil",
  "16" = "Interstitial macrophage",
  "17" = "Mast cell",
  "18" = "Mesothelial",
  "19" = "Cycling",
  "20" = "Mixed alveolar epithelial",
  "21" = "Mrc1+ macrophage",
  "22" = "pDC/DC",
  "23" = "Cycling",
  "24" = "Mixed/doublet"
)

rat$cell_type <- unname(
  cluster_labels[as.character(rat$seurat_clusters)]
)


DimPlot(
  rat,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)



# ---------------------------------------------
# Find markers for annotated cell types
# ---------------------------------------------

query_cell_types <- c(
  "ATI",
  "ATII",
  "Airway epithelial",
  "Sox9+ epithelial",
  "Capillary endothelial",
  "Vascular endothelial",
  "VCAM1+ endothelial",
  "Fibroblast",
  "Smooth muscle",
  "Mesothelial",
  "Alveolar macrophage",
  "Interstitial macrophage",
  "Mrc1+ macrophage",
  "Monocyte",
  "Dendritic cell",
  "pDC/DC",
  "Neutrophil",
  "Mast cell",
  "T cell",
  "B cell",
  "NK cell"
)

rat_marker <- subset(
  rat,
  subset = cell_type %in% query_cell_types
)

DefaultAssay(rat_marker) <- "RNA"

Idents(rat_marker) <- "cell_type"

table(Idents(rat_marker))


celltype_markers <- FindAllMarkers(
  rat_marker,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25
)


dim(celltype_markers)

table(celltype_markers$cluster)

head(celltype_markers)

write.csv(
  celltype_markers,
  "Raredon_rat_celltype_FindAllMarkers.csv",
  row.names = FALSE
)


celltype_markers_specific <- celltype_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.5,
    pct.1 >= 0.25,
    (pct.1 - pct.2) >= 0.15
  )


table(celltype_markers_specific$cluster)

celltype_markers_protein <- celltype_markers_specific %>%
  filter(gene %in% expr$Gene.symbol)

table(celltype_markers_protein$cluster)

#VCAM1+ endothelial     Alveolar macrophage                     ATI    Vascular endothelial 
#136                     471                     328                      87 
#T cell                 NK cell                Monocyte                  B cell 
#66                      78                     146                      34 
#ATII               Mast cell       Airway epithelial          Dendritic cell 
#167                      45                     362                     231 
#Sox9+ epithelial              Neutrophil        Mrc1+ macrophage Interstitial macrophage 
#293                      47                      45                     131 
#Capillary endothelial             Mesothelial              Fibroblast           Smooth muscle 
#153                     483                     200                      99 
#pDC/DC 
#279


celltype_marker_top30 <- celltype_markers_protein %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 30,
    with_ties = FALSE
  ) %>%
  ungroup()

celltype_marker_top30 %>%
  dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2) %>%
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

write.csv(
  celltype_marker_top30,
  "Raredon_rat_celltype_top30_proteomics_markers.csv",
  row.names = FALSE
)

expr1 <- expr[expr$Gene.symbol != "-", ] 
expr_mat <- as.data.frame(expr1[, -c(1:7)])
rownames(expr_mat) <- expr1$Gene.symbol

expr_z <- t(scale(t(expr_mat)))

marker_list <- split(
  celltype_marker_top30$gene,
  celltype_marker_top30$cluster
)

celltype_score <- sapply(
  marker_list,
  function(markers) {
    genes <- intersect(markers, rownames(expr_z))
    colMeans(expr_z[genes, , drop = FALSE], na.rm = TRUE)
  }
)

celltype_score <- as.data.frame(celltype_score)
celltype_score$Sample.ID <- rownames(celltype_score)

celltype_score <- merge(
  meta,
  celltype_score,
  by = "Sample.ID"
)



library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)



celltype_cols <- names(marker_list)

length(celltype_cols)


celltype_lmm_results <- list()

for (ct in celltype_cols) {
  
  dat <- celltype_score %>%
    dplyr::select(PID, Age, Sex, Region, all_of(ct)) %>%
    dplyr::rename(Score = all_of(ct))
  
  fit <- lmer(
    Score ~ Age + Sex + Region + (1 | PID),
    data = dat
  )
  
  emm <- emmeans(
    fit,
    ~ Region
  )
  
  region_pairs <- contrast(
    emm,
    method = "pairwise",
    adjust = "none"
  )
  
  tmp <- as.data.frame(region_pairs)
  
  tmp$Cell.type <- ct
  
  celltype_lmm_results[[ct]] <- tmp
}

celltype_region_results <- bind_rows(
  celltype_lmm_results
)


celltype_region_results <- celltype_region_results %>%
  group_by(contrast) %>%
  mutate(
    FDR = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup()


celltype_region_results %>%
  dplyr::select(
    Cell.type,
    contrast,
    estimate,
    SE,
    df,
    p.value,
    FDR
  ) %>%
  arrange(FDR)




celltype_region_omnibus <- list()

for (ct in celltype_cols) {
  
  dat <- celltype_score %>%
    dplyr::select(PID, Age, Sex, Region, all_of(ct)) %>%
    dplyr::rename(Score = all_of(ct))
  
  fit <- lmer(
    Score ~ Age + Sex + Region + (1 | PID),
    data = dat
  )
  
  aov_tab <- anova(fit)
  
  tmp <- data.frame(
    Cell.type = ct,
    F.value = aov_tab["Region", "F value"],
    df1 = aov_tab["Region", "NumDF"],
    df2 = aov_tab["Region", "DenDF"],
    p.value = aov_tab["Region", "Pr(>F)"]
  )
  
  celltype_region_omnibus[[ct]] <- tmp
}

celltype_region_omnibus <- bind_rows(
  celltype_region_omnibus
)


celltype_region_omnibus <- celltype_region_omnibus %>%
  mutate(
    FDR = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(FDR)


celltype_region_omnibus <- celltype_region_omnibus %>%
  mutate(
    Significant = FDR < 0.05
  )

celltype_region_omnibus %>%
  filter(FDR < 0.05)



# ============================================================
# Post-hoc pairwise comparisons for Airway epithelial only
# ============================================================

fit_airway <- lmer(
  `Airway epithelial` ~ Age + Sex + Region + (1 | PID),
  data = celltype_score
)

emm_airway <- emmeans(
  fit_airway,
  ~ Region
)

pair_airway <- pairs(
  emm_airway,
  adjust = "tukey"
)

pair_airway_df <- as.data.frame(pair_airway)

pair_airway_df



pair_airway_df <- pair_airway_df %>%
  mutate(
    p.signif = case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001  ~ "***",
      p.value < 0.01   ~ "**",
      p.value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )
  )

pair_airway_df


celltype_score$Region <- factor(
  celltype_score$Region,
  levels = c("Bottom", "Central", "Peripheral", "Top")
)


library(ggplot2)

p_airway <- ggplot(
  celltype_score,
  aes(
    x = Region,
    y = `Airway epithelial`
  )
) +
  geom_line(
    aes(group = PID),
    alpha = 0.20,
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.35
  ) +
  geom_point(
    size = 2,
    alpha = 0.75
  ) +
  theme_bw() +
  labs(
    x = "Lung region",
    y = "Airway epithelial marker score"
  ) +
  theme(
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

p_airway


library(ggpubr)

stat_airway <- pair_airway_df %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".* - ", "", contrast)
  )

stat_airway_sig <- stat_airway %>%
  filter(p.value < 0.05)


stat_airway_sig$y.position <- seq(
  1.65,
  by = 0.18,
  length.out = nrow(stat_airway_sig)
)


p_airway <- p_airway +
  stat_pvalue_manual(
    stat_airway_sig,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  ) +
  coord_cartesian(
    ylim = c(-1.2, 2.5),
    clip = "off"
  )

p_airway

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(tidyr)
library(pheatmap)


# ============================================================
# 1. LMM + emmeans: estimate regional mean after Age/Sex adjustment
# ============================================================

emm_list <- list()

for (ct in celltype_cols) {
  
  dat <- celltype_score %>%
    dplyr::select(
      PID,
      Age,
      Sex,
      Region,
      all_of(ct)
    ) %>%
    dplyr::rename(
      Score = all_of(ct)
    )
  
  fit <- lmer(
    Score ~ Age + Sex + Region + (1 | PID),
    data = dat
  )
  
  emm <- emmeans(
    fit,
    ~ Region
  )
  
  emm_df <- as.data.frame(emm)
  
  emm_df$Cell.type <- ct
  
  emm_list[[ct]] <- emm_df
}


emm_region <- bind_rows(emm_list)



# ============================================================
# 2. Convert emmeans to heatmap matrix
# ============================================================

heatmap_df <- emm_region %>%
  dplyr::select(
    Cell.type,
    Region,
    emmean
  ) %>%
  pivot_wider(
    names_from = Region,
    values_from = emmean
  )


heatmap_mat <- as.data.frame(heatmap_df)

rownames(heatmap_mat) <- heatmap_mat$Cell.type

heatmap_mat$Cell.type <- NULL



# ============================================================
# 3. Reorder cell types
# ============================================================

cell_order <- c(
  "ATI",
  "ATII",
  "Airway epithelial",
  "Sox9+ epithelial",
  
  "Fibroblast",
  "Smooth muscle",
  "Mesothelial",
  
  "Capillary endothelial",
  "Vascular endothelial",
  "VCAM1+ endothelial",
  
  "Alveolar macrophage",
  "Interstitial macrophage",
  "Mrc1+ macrophage",
  "Monocyte",
  "Dendritic cell",
  "pDC/DC",
  "Neutrophil",
  "Mast cell",
  
  "T cell",
  "B cell",
  "NK cell"
)

heatmap_mat <- heatmap_mat[cell_order, ]



# ============================================================
# 4. Add FDR significance mark
# ============================================================

fdr_vec <- celltype_region_omnibus$FDR[
  match(
    rownames(heatmap_mat),
    celltype_region_omnibus$Cell.type
  )
]


rownames(heatmap_mat) <- ifelse(
  fdr_vec < 0.05,
  paste0(rownames(heatmap_mat), " *"),
  rownames(heatmap_mat)
)



# ============================================================
# 5. Plot heatmap
# ============================================================
library(pheatmap)
pheatmap(
  as.matrix(heatmap_mat),
  
  # show relative regional pattern within each cell type
  scale = "none",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  fontsize_row = 10,
  fontsize_col = 12,
  
  border_color = NA,
  
  angle_col = 45
)


write.csv(
  emm_region,
  "celltype_region_emmeans.csv",
  row.names = FALSE
)

write.csv(
  celltype_region_omnibus,
  "celltype_region_LMM_FDR.csv",
  row.names = FALSE
)
