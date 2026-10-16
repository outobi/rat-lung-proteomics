# =============================================================================
# Script 9 (scRNAseq cell type marker score in regions), step 1:
#   Rat lung scRNA-seq reference and cell-type marker genes
# =============================================================================
#
# Purpose:
#   Build an annotated rat lung single-cell reference from the public
#   Raredon et al. dataset (GEO GSM3926541 male rat, GSM3926542 female rat)
#   and find marker genes for each annotated cell type. The marker table
#   written here is the only input that step 2
#   (Script_9_step2_cell_type_marker_score_in_regions.R) needs from the
#   scRNA-seq data, so step 2 can be rerun without repeating this step.
#
# Workflow:
#   1. Read the raw count matrices (5,029 male and 5,924 female cells), add
#      percent mitochondrial reads, log-normalize, and find 2,000 variable
#      genes per sample.
#   2. Integrate male and female with CCA anchors (features = genes among the
#      top 1,000 variable genes of both samples; dims 1:30), scale while
#      regressing out nCount_RNA and percent.mt, PCA (30 PCs).
#   3. Cluster (dims 1:20, resolution 0.5; 25 clusters), t-SNE and UMAP.
#   4. Annotate clusters manually from canonical markers (DotPlot of
#      marker_panel) into 21 cell types; clusters labelled "Cycling",
#      "Mixed alveolar epithelial" or "Mixed/doublet" are excluded.
#   5. FindAllMarkers per cell type (RNA assay, Wilcoxon, positive markers,
#      min.pct 0.10, logfc.threshold 0.25).
#
# IMPORTANT - reproducibility of the cluster labels:
#   cluster_labels assigns cell types by cluster NUMBER. Cluster numbering
#   depends on the Seurat version and random seeds. The reported analysis
#   used Seurat 5.4.0 / SeuratObject 5.3.0 (R 4.4.3). If this step is rerun
#   with another version, check the DotPlot and UMAP before trusting the
#   labels; the reported marker table is supplied, so step 2 does not need
#   this step to be rerun.
#
# Inputs (in this folder; downloaded from GEO, series GSE133747):
#   - GSM3926541_mRat.raw.data.mtx.gz, GSM3926541_mRat.genes.tsv.gz,
#     GSM3926541_mRat.barcodes.tsv.gz   (male rat)
#   - GSM3926542_fRat.raw.data.mtx.gz, GSM3926542_fRat.genes.tsv.gz,
#     GSM3926542_fRat.barcodes.tsv.gz   (female rat)
#
# Outputs (written to this folder):
#   - Raredon_rat_celltype_FindAllMarkers.csv  (input of step 2)
#       Marker genes per annotated cell type (gene, cluster = cell type,
#       avg_log2FC, pct.1, pct.2, p_val, p_val_adj).
#   - UMAP cell type.pdf
#       UMAP colored by annotated cell type.
#   - Marker panel dotplot by cluster.pdf
#       Canonical marker expression per cluster, used for the annotation.
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# Runtime: several minutes (integration and marker search); needs several GB
# of memory.
# =============================================================================

# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages.
# Matrix: readMM(); Seurat: single-cell workflow; dplyr: pipes;
# ggplot2: saving plots.
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)


# Step 1: Read the raw count matrices ------------------------------------------

m_counts <- readMM("GSM3926541_mRat.raw.data.mtx.gz")
m_genes <- read.delim("GSM3926541_mRat.genes.tsv.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
m_barcodes <- read.delim("GSM3926541_mRat.barcodes.tsv.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

f_counts <- readMM("GSM3926542_fRat.raw.data.mtx.gz")
f_genes <- read.delim("GSM3926542_fRat.genes.tsv.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
f_barcodes <- read.delim("GSM3926542_fRat.barcodes.tsv.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

dim(m_counts) # genes x 5,029 cells
dim(f_counts) # genes x 5,924 cells

# Male and female files share the same gene list; no barcodes are shared.
stopifnot(identical(m_genes$x, f_genes$x))
gene_names <- make.unique(m_genes$x)

rownames(m_counts) <- gene_names
colnames(m_counts) <- m_barcodes$x
rownames(f_counts) <- gene_names
colnames(f_counts) <- f_barcodes$x


# Step 2: Seurat objects, QC metric, normalization -----------------------------

m_rat <- CreateSeuratObject(counts = m_counts, project = "mRat", min.cells = 0, min.features = 0)
m_rat$Sex <- "Male"

f_rat <- CreateSeuratObject(counts = f_counts, project = "fRat", min.cells = 0, min.features = 0)
f_rat$Sex <- "Female"

# Percent mitochondrial reads (11 "Mt-" genes).
mt_genes <- grep("^Mt-", rownames(m_rat), value = TRUE)
length(mt_genes) # 11

m_rat[["percent.mt"]] <- PercentageFeatureSet(m_rat, features = mt_genes)
f_rat[["percent.mt"]] <- PercentageFeatureSet(f_rat, features = mt_genes)

m_rat <- NormalizeData(m_rat)
m_rat <- FindVariableFeatures(m_rat, selection.method = "vst", nfeatures = 2000)

f_rat <- NormalizeData(f_rat)
f_rat <- FindVariableFeatures(f_rat, selection.method = "vst", nfeatures = 2000)


# Step 3: Integration, PCA, clustering, UMAP -----------------------------------

# Integration features: genes in the top 1,000 variable genes of both samples.
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

rat <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(rat) <- "integrated"

rat <- ScaleData(rat, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
rat <- RunPCA(rat, features = integration_features, npcs = 30, verbose = FALSE)
ElbowPlot(rat, ndims = 30)

rat <- FindNeighbors(rat, dims = 1:20)
rat <- FindClusters(rat, resolution = 0.5) # 25 clusters (0-24) in the reported run
rat <- RunTSNE(rat, dims = 1:20)
rat <- RunUMAP(rat, dims = 1:20)

table(Idents(rat))

# Use the RNA assay (joined layers) for marker expression.
DefaultAssay(rat) <- "RNA"
rat <- JoinLayers(rat, assay = "RNA")


# Step 4: Manual cell-type annotation ------------------------------------------

# Canonical markers used to annotate the clusters.
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

p_dot <- DotPlot(rat, features = marker_panel, group.by = "seurat_clusters") + RotatedAxis()
p_dot

# Cluster number -> cell type (valid for the reported clustering only; see
# the reproducibility note in the header).
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
stopifnot(all(levels(rat$seurat_clusters) %in% names(cluster_labels)))

rat$cell_type <- unname(cluster_labels[as.character(rat$seurat_clusters)])

p_umap <- DimPlot(rat, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE)
p_umap

ggsave("Marker panel dotplot by cluster.pdf", p_dot, width = 16, height = 8.5)
ggsave("UMAP cell type.pdf", p_umap, width = 11, height = 8.5)


# Step 5: Marker genes per annotated cell type ---------------------------------

# The 21 cell types used downstream (Cycling and Mixed clusters excluded).
query_cell_types <- c(
  "ATI", "ATII", "Airway epithelial", "Sox9+ epithelial",
  "Capillary endothelial", "Vascular endothelial", "VCAM1+ endothelial",
  "Fibroblast", "Smooth muscle", "Mesothelial",
  "Alveolar macrophage", "Interstitial macrophage", "Mrc1+ macrophage",
  "Monocyte", "Dendritic cell", "pDC/DC", "Neutrophil", "Mast cell",
  "T cell", "B cell", "NK cell"
)

rat_marker <- subset(rat, subset = cell_type %in% query_cell_types)
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

write.csv(celltype_markers, "Raredon_rat_celltype_FindAllMarkers.csv", row.names = FALSE)

# The Seurat object is not saved (several GB); record the package versions.
sessionInfo()
