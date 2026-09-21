setwd(dirname(rstudioapi::getSourceEditorContext()$path)) ## set working dir as where the script locates


library(Matrix)
library(Seurat)


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


# region speicific proteins-------
library(readxl)
library(dplyr)
library(openxlsx)

meta <- read_excel("sample metadata.xlsx")
expr <- read_excel("expression matrix.xlsx")
expr<- as.data.frame(expr)
rownames(expr) <- expr$Uniprot.ID

# Young / Aged sample IDs
young_ids <- meta %>%
  filter(Age == "Young") %>%
  pull(Sample.ID)

aged_ids <- meta %>%
  filter(Age == "Aged") %>%
  pull(Sample.ID)

length(young_ids)
length(aged_ids)


expr_young <- expr[, colnames(expr) %in% c(young_ids)]
expr_aged  <- expr[, colnames(expr) %in% c(aged_ids)]


meta_young <- meta %>%
  dplyr::filter(Age == "Young")

meta_aged <- meta %>%
  dplyr::filter(Age == "Aged")

all(colnames(expr_young) == meta_young$Sample.ID)
# TRUE
all(colnames(expr_aged) == meta_aged$Sample.ID)
# TRUE


## region vs all other regions based on Wilcoxon rank-sum test
wilcox_region_vs_rest <- function(target_region,
                                  expr_mat,
                                  metadata) {
  
  wilcox_p <- c()
  logFC <- c()
  
  for (i in seq_len(nrow(expr_mat))) {
    
    x <- as.numeric(
      expr_mat[i, metadata$Region == target_region]
    )
    
    y <- as.numeric(
      expr_mat[i, metadata$Region != target_region]
    )
    
    wilcox_p[i] <- wilcox.test(
      x = x,
      y = y
    )$p.value
    
    logFC[i] <- mean(x) - mean(y)
  }
  
  result <- data.frame(
    Uniprot.ID = rownames(expr_mat),
    p.value = wilcox_p,
    FDR = p.adjust(wilcox_p, method = "BH"),
    Log2FC = logFC
  )
  
  return(result)
}



young_bottom <- wilcox_region_vs_rest(
  target_region = "Bottom",
  expr_mat = expr_young,
  metadata = meta_young
)

young_central <- wilcox_region_vs_rest(
  target_region = "Central",
  expr_mat = expr_young,
  metadata = meta_young
)

young_peripheral <- wilcox_region_vs_rest(
  target_region = "Peripheral",
  expr_mat = expr_young,
  metadata = meta_young
)

young_top <- wilcox_region_vs_rest(
  target_region = "Top",
  expr_mat = expr_young,
  metadata = meta_young
)



aged_bottom <- wilcox_region_vs_rest(
  target_region = "Bottom",
  expr_mat = expr_aged,
  metadata = meta_aged
)

aged_central <- wilcox_region_vs_rest(
  target_region = "Central",
  expr_mat = expr_aged,
  metadata = meta_aged
)

aged_peripheral <- wilcox_region_vs_rest(
  target_region = "Peripheral",
  expr_mat = expr_aged,
  metadata = meta_aged
)

aged_top <- wilcox_region_vs_rest(
  target_region = "Top",
  expr_mat = expr_aged,
  metadata = meta_aged
)


region_specific_feature <- list(
  "Young_Bottom_vs_all"      = young_bottom,
  "Young_Central_vs_all"     = young_central,
  "Young_Peripheral_vs_all"  = young_peripheral,
  "Young_Top_vs_all"         = young_top,
  
  "Aged_Bottom_vs_all"       = aged_bottom,
  "Aged_Central_vs_all"      = aged_central,
  "Aged_Peripheral_vs_all"   = aged_peripheral,
  "Aged_Top_vs_all"          = aged_top
)

export(region_specific_feature, "region_specific_protein_summary.xlsx")



log2fc_cutoff <- 0.25
fdr_cutoff <- 0.05


young_bottom_RSP <- young_bottom[
  young_bottom$Log2FC > log2fc_cutoff &
    young_bottom$FDR < fdr_cutoff, ]

young_central_RSP <- young_central[
  young_central$Log2FC > log2fc_cutoff &
    young_central$FDR < fdr_cutoff, ]

young_peripheral_RSP <- young_peripheral[
  young_peripheral$Log2FC > log2fc_cutoff &
    young_peripheral$FDR < fdr_cutoff, ]

young_top_RSP <- young_top[
  young_top$Log2FC > log2fc_cutoff &
    young_top$FDR < fdr_cutoff, ]


aged_bottom_RSP <- aged_bottom[
  aged_bottom$Log2FC > log2fc_cutoff &
    aged_bottom$FDR < fdr_cutoff, ]

aged_central_RSP <- aged_central[
  aged_central$Log2FC > log2fc_cutoff &
    aged_central$FDR < fdr_cutoff, ]

aged_peripheral_RSP <- aged_peripheral[
  aged_peripheral$Log2FC > log2fc_cutoff &
    aged_peripheral$FDR < fdr_cutoff, ]

aged_top_RSP <- aged_top[
  aged_top$Log2FC > log2fc_cutoff &
    aged_top$FDR < fdr_cutoff, ]

c(
  Young_Bottom     = nrow(young_bottom_RSP),
  Young_Central    = nrow(young_central_RSP),
  Young_Peripheral = nrow(young_peripheral_RSP),
  Young_Top        = nrow(young_top_RSP),
  
  Aged_Bottom      = nrow(aged_bottom_RSP),
  Aged_Central     = nrow(aged_central_RSP),
  Aged_Peripheral  = nrow(aged_peripheral_RSP),
  Aged_Top         = nrow(aged_top_RSP)
)

#Young_Bottom    Young_Central Young_Peripheral        Young_Top      Aged_Bottom 
#0              113              103                0                1 
#Aged_Central  Aged_Peripheral         Aged_Top 
#19                3                1 

intersect(young_central_RSP$gene.symbol, young_peripheral_RSP$gene.symbol) # 0

# Young Central, n = 112, remove one protein segment,no gene symbol
file2 <- expr[match(young_central_RSP$Uniprot.ID, expr$Uniprot.ID),]
young_central_RSP <- cbind("Gene.symbol" = file2$Gene.symbol, young_central_RSP)
young_central_RSP <- young_central_RSP[young_central_RSP$Gene.symbol != "-",]

# Young Peripheral, n = 103
file2 <- expr[match(young_peripheral_RSP$Uniprot.ID, expr$Uniprot.ID), ]
young_peripheral_RSP <- cbind("Gene.symbol" = file2$Gene.symbol,young_peripheral_RSP)
young_peripheral_RSP <- young_peripheral_RSP[young_peripheral_RSP$Gene.symbol != "-", ]


# Aged Bottom, n = 1
file2 <- expr[match(aged_bottom_RSP$Uniprot.ID, expr$Uniprot.ID), ]
aged_bottom_RSP <- cbind("Gene.symbol" = file2$Gene.symbol,aged_bottom_RSP)
aged_bottom_RSP <- aged_bottom_RSP[aged_bottom_RSP$Gene.symbol != "-", ]


# Aged Central, n = 19
file2 <- expr[match(aged_central_RSP$Uniprot.ID, expr$Uniprot.ID), ]
aged_central_RSP <- cbind("Gene.symbol" = file2$Gene.symbol,aged_central_RSP)
aged_central_RSP <- aged_central_RSP[aged_central_RSP$Gene.symbol != "-", ]


# Aged Peripheral, n = 3
file2 <- expr[match(aged_peripheral_RSP$Uniprot.ID, expr$Uniprot.ID), ]
aged_peripheral_RSP <- cbind("Gene.symbol" = file2$Gene.symbol,aged_peripheral_RSP)
aged_peripheral_RSP <- aged_peripheral_RSP[aged_peripheral_RSP$Gene.symbol != "-", ]


# Aged Top, n = 1
file2 <- expr[match(aged_top_RSP$Uniprot.ID, expr$Uniprot.ID), ]
aged_top_RSP <- cbind("Gene.symbol" = file2$Gene.symbol,aged_top_RSP)
aged_top_RSP <- aged_top_RSP[aged_top_RSP$Gene.symbol != "-", ]


# Young Central, n = 105 overlap with scRNAseq f_gene list
young_central_RSP <- young_central_RSP[
  young_central_RSP$Gene.symbol %in% f_genes$x, ]

# Young Peripheral n = 92 overlap with scRNAseq f_gene list
young_peripheral_RSP <- young_peripheral_RSP[
  young_peripheral_RSP$Gene.symbol %in% f_genes$x, ]

# Aged Bottom, n = 1 overlap with scRNAseq f_gene list
aged_bottom_RSP <- aged_bottom_RSP[
  aged_bottom_RSP$Gene.symbol %in% f_genes$x, ]

# Aged Central, n = 17 overlap with scRNAseq f_gene list
aged_central_RSP <- aged_central_RSP[
  aged_central_RSP$Gene.symbol %in% f_genes$x, ]

# Aged Peripheral, n = 3 overlap with scRNAseq f_gene list
aged_peripheral_RSP <- aged_peripheral_RSP[
  aged_peripheral_RSP$Gene.symbol %in% f_genes$x, ]

# Aged Top, n = 1 overlap with scRNAseq f_gene list
aged_top_RSP <- aged_top_RSP[
  aged_top_RSP$Gene.symbol %in% f_genes$x, ]


region_specific_proteins <- list(
  "Young_Bottom"      = young_bottom_RSP,
  "Young_Central"     = young_central_RSP,
  "Young_Peripheral"  = young_peripheral_RSP,
  "Young_Top"         = young_top_RSP,
  
  "Aged_Bottom"       = aged_bottom_RSP,
  "Aged_Central"      = aged_central_RSP,
  "Aged_Peripheral"   = aged_peripheral_RSP,
  "Aged_Top"          = aged_top_RSP
)

export(region_specific_proteins, "region_specific_upregulated_proteins.xlsx")



# query

# combine male + female raw counts
rat_counts <- cbind(m_counts, f_counts)

# metadata
rat_meta <- rat@meta.data
query_seurat <- CreateSeuratObject(counts = rat_counts, meta.data = rat_meta)

query_seurat <- NormalizeData(
  object = query_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

normalized_count <- GetAssayData(
  query_seurat,
  assay = "RNA",
  layer = "data"
)



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

query_signature <- function(signature_genes,
                            normalized_count,
                            cell_type_annotation,
                            query_cell_types) {
  
 
  # extract signature expression
  test <- as.matrix(
    normalized_count[signature_genes, , drop = FALSE]
  )
  
  # gene-wise z-score
  gene_mean <- apply(test, 1, mean)
  gene_sd <- apply(test, 1, sd)
  
  z_test <- sweep(test, 1, gene_mean, "-")
  z_test <- sweep(z_test, 1, gene_sd, "/")
  
  # calculate score for each cell type
  score <- numeric(length(query_cell_types))
  
  for (i in seq_along(query_cell_types)) {
    
    cells <- which(
      cell_type_annotation == query_cell_types[i]
    )
    
    cell_data <- z_test[, cells, drop = FALSE]
    
    score[i] <- mean(cell_data)
  }
  
  result <- data.frame(
    Cell.type = query_cell_types,
    Score = score,
    n_genes = length(signature_genes)
  )
  
  result <- result[
    order(result$Score, decreasing = TRUE),
  ]
  
  return(result)
}

young_central_score <- query_signature(
  signature_genes = young_central_RSP$Gene.symbol,
  normalized_count = normalized_count,
  cell_type_annotation = rat$cell_type,
  query_cell_types = query_cell_types
)


young_peripheral_score <- query_signature(
  young_peripheral_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  query_cell_types
)

aged_bottom_score <- query_signature(
  aged_bottom_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  query_cell_types
)

aged_central_score <- query_signature(
  aged_central_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  query_cell_types
)

aged_peripheral_score <- query_signature(
  aged_peripheral_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  query_cell_types
)

aged_top_score <- query_signature(
  aged_top_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  query_cell_types
)



# plot

library(reshape2)
library(ggplot2)

query.summary <- cbind(
  young_central_score$Score,
  young_peripheral_score$Score,
  aged_bottom_score$Score,
  aged_central_score$Score,
  aged_peripheral_score$Score,
  aged_top_score$Score
)

query.summary <- as.data.frame(query.summary)

colnames(query.summary) <- c(
  "Young Central",
  "Young Peripheral",
  "Aged Bottom",
  "Aged Central",
  "Aged Peripheral",
  "Aged Top"
)

rownames(query.summary) <- query_cell_types


new_order <- c(
  # epithelial
  "ATI",
  "ATII",
  "Airway epithelial",
  "Sox9+ epithelial",
  
  # mesenchymal
  "Fibroblast",
  "Smooth muscle",
  "Mesothelial",
  
  # endothelial
  "Capillary endothelial",
  "Vascular endothelial",
  "VCAM1+ endothelial",
  
  # myeloid
  "Alveolar macrophage",
  "Interstitial macrophage",
  "Mrc1+ macrophage",
  "Monocyte",
  "Dendritic cell",
  "pDC/DC",
  "Neutrophil",
  "Mast cell",
  
  # lymphoid
  "T cell",
  "B cell",
  "NK cell"
)

query.summary <- query.summary[new_order, ]

query.plot <- query.summary

query.plot <- cbind(
  query.plot,
  Cell.type = rownames(query.plot)
)

query.plot <- query.plot[
  nrow(query.plot):1,
]

pcm <- melt(
  query.plot,
  id.vars = "Cell.type"
)

pcm$Cell.type <- factor(
  pcm$Cell.type,
  levels = unique(pcm$Cell.type)
)


p1 <- ggplot(
  pcm,
  aes(
    y = Cell.type,
    x = variable
  )
) +
  geom_point(
    aes(
      size = value,
      fill = variable
    ),
    shape = 21
  ) +
  scale_size_continuous(
    limits = c(0, max(pcm$value)),
    range = c(1,10)
  ) +
  theme(
    legend.key = element_blank(),
    
    axis.text.x =
      element_text(
        colour="black",
        size=12,
        face="bold",
        angle=45,
        hjust=1
      ),
    
    axis.text.y =
      element_text(
        colour="black",
        size=11,
        face="bold"
      ),
    
    axis.title =
      element_text(
        face="bold",
        size=13
      )
  ) +
  labs(
    x="Regional signatures",
    y="Cell Types"
  )

p1


pcm_young <- pcm %>%
  dplyr::filter(
    variable %in% c("Young Central", "Young Peripheral")
  )

pcm_young$variable <- factor(
  pcm_young$variable,
  levels = c("Young Central", "Young Peripheral")
)




p11 <- ggplot(
  pcm_young,
  aes(
    y = Cell.type,
    x = variable
  )
) +
  geom_point(
    aes(
      size = value,
      fill = variable
    ),
    shape = 21
  ) +
  scale_size_continuous(
    limits = c(0, 0.2),
    range = c(1,10),
   
  ) +
  theme(
    legend.key = element_blank(),
    
    axis.text.x =
      element_text(
        colour="black",
        size=12,
        face="bold",
        angle=45,
        hjust=1
      ),
    
    axis.text.y =
      element_text(
        colour="black",
        size=11,
        face="bold"
      ),
    
    axis.title =
      element_text(
        face="bold",
        size=13
      )
  ) +
  labs(
    x="Regional signatures",
    y="Cell Types"
  )

p11


pcm_depletion <- pcm[pcm$value < 0, ]
p2 <- ggplot(
  pcm_depletion,
  aes(
    y = Cell.type,
    x = variable
  )
) +
  geom_point(
    aes(
      size = abs(value),
      fill = variable
    ),
    shape = 21
  ) +
  scale_size_continuous(
    limits = c(0, max(abs(pcm_depletion$value))),
    range = c(1,10)
  ) +
  theme(
    legend.key = element_blank(),
    
    axis.text.x =
      element_text(
        colour = "black",
        size = 12,
        face = "bold",
        angle = 45,
        hjust = 1
      ),
    
    axis.text.y =
      element_text(
        colour = "black",
        size = 11,
        face = "bold"
      ),
    
    axis.title =
      element_text(
        face = "bold",
        size = 13
      )
  ) +
  labs(
    x = "Regional signatures",
    y = "Cell Types"
  )

p2

pcm_depletion_young <- pcm %>%
  dplyr::filter(
    variable %in% c("Young Central", "Young Peripheral"),
    value < 0
  )

pcm_depletion_young$variable <- factor(
  pcm_depletion_young$variable,
  levels = c("Young Central", "Young Peripheral")
)

p22 <- ggplot(
  pcm_depletion_young,
  aes(
    y = Cell.type,
    x = variable
  )
) +
  geom_point(
    aes(
      size = abs(value),
      fill = variable
    ),
    shape = 21
  ) +  scale_y_discrete(drop = FALSE) +
  scale_size_continuous(
    limits = c(0, 0.2),
    range = c(1, 10)
  ) +
  theme(
    legend.key = element_blank(),
    
    axis.text.x = element_text(
      colour = "black",
      size = 12,
      face = "bold",
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      colour = "black",
      size = 11,
      face = "bold"
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 13
    )
  ) +
  labs(
    x = "Regional signatures",
    y = "Cell Types"
  )

p22





# larger cell type enrichment

celltype_group <- list(
  Epithelial = c(
    "ATI",
    "ATII",
    "Airway epithelial",
    "Sox9+ epithelial"
  ),
  
  Endothelial = c(
    "Capillary endothelial",
    "Vascular endothelial",
    "VCAM1+ endothelial"
  ),
  
  Stromal = c(
    "Fibroblast",
    "Smooth muscle",
    "Mesothelial"
  ),
  
  Myeloid = c(
    "Alveolar macrophage",
    "Interstitial macrophage",
    "Mrc1+ macrophage",
    "Monocyte",
    "Dendritic cell",
    "pDC/DC",
    "Neutrophil",
    "Mast cell"
  ),
  
  Lymphoid = c(
    "T cell",
    "B cell",
    "NK cell"
  )
)


query_large_celltype <- function(signature_genes,
                                 normalized_count,
                                 cell_type_annotation,
                                 group_list){
  
  signature_genes <- intersect(
    signature_genes,
    rownames(normalized_count)
  )
  
  test <- as.matrix(
    normalized_count[signature_genes, , drop=FALSE]
  )
  
  # gene-wise z score
  z_test <- t(scale(t(test)))
  
  scores <- numeric(length(group_list))
  
  for(i in seq_along(group_list)){
    
    cells <- which(
      cell_type_annotation %in% group_list[[i]]
    )
    
    scores[i] <- mean(
      z_test[, cells]
    )
  }
  
  data.frame(
    Cell.type = names(group_list),
    Score = scores
  )
}


large_scores_young_central <- query_large_celltype(
  young_central_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  celltype_group
)

large_scores_young_peripheral <- query_large_celltype(
  young_peripheral_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  celltype_group
)

large_scores_aged_bottom <- query_large_celltype(
  aged_bottom_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  celltype_group
)

large_scores_aged_central <- query_large_celltype(
  aged_central_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  celltype_group
)

large_scores_aged_peripheral <- query_large_celltype(
  aged_peripheral_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  celltype_group
)

large_scores_aged_top <- query_large_celltype(
  aged_top_RSP$Gene.symbol,
  normalized_count,
  rat$cell_type,
  celltype_group
)


sapply(celltype_group, function(x)
  +     sum(rat$cell_type %in% x)
  + )
#Epithelial Endothelial     Stromal     Myeloid    Lymphoid 
#3333        1650         702        3147        1823 


large_query_summary <- data.frame(
  Cell.type = large_scores_young_central$Cell.type,
  
  `Young Central` =
    large_scores_young_central$Score,
  
  `Young Peripheral` =
    large_scores_young_peripheral$Score,
  
  `Aged Bottom` =
    large_scores_aged_bottom$Score,
  
  `Aged Central` =
    large_scores_aged_central$Score,
  
  `Aged Peripheral` =
    large_scores_aged_peripheral$Score,
  
  `Aged Top` =
    large_scores_aged_top$Score
)

rownames(large_query_summary) <- large_query_summary$Cell.type



plot_df <- large_query_summary

plot_df <- melt(
  plot_df,
  id.vars = "Cell.type"
)

colnames(plot_df)[2:3] <- c("Region", "Score")


plot_enrich <- plot_df[
  plot_df$Score > 0,
]

plot_enrich$Cell.type <- factor(
  plot_enrich$Cell.type,
  levels = rev(unique(plot_enrich$Cell.type))
)

library(ggplot2)

p_large_enrich <- ggplot(
  plot_enrich,
  aes(
    x = Region,
    y = Cell.type
  )
) +
  geom_point(
    aes(
      size = Score,
      fill = Region
    ),
    shape = 21
  ) +
  scale_size_continuous(
    range = c(2,10)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    axis.title = element_blank()
  )

p_large_enrich


plot_deplete <- plot_df[
  plot_df$Score < 0,
]

plot_deplete$Cell.type <- factor(
  plot_deplete$Cell.type,
  levels = rev(unique(plot_deplete$Cell.type))
)


p_large_deplete <- ggplot(
  plot_deplete,
  aes(
    x = Region,
    y = Cell.type
  )
) +
  geom_point(
    aes(
      size = abs(Score),
      fill = Region
    ),
    shape = 21
  ) +
  scale_size_continuous(
    range = c(2,10)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    axis.title = element_blank()
  )

p_large_deplete








