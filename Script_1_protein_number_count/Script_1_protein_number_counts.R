# =============================================================================
# Script 1: Protein group counts per sample, split by matrisome division
# =============================================================================
#
# Purpose:
#   Count the number of detected protein groups in each of the 96 lung samples,
#   for (1) all proteins, (2) all matrisome proteins, (3) core matrisome, and
#   (4) matrisome-associated proteins, then arrange the counts by experimental
#   group (16 groups x 6 samples).
#
# Inputs (all in this folder):
#   - report.pg_matrix.tsv
#       DIA-NN protein-group matrix. Columns 1-4 are Protein.Group,
#       Protein.Names, Genes, First.Protein.Description; columns 5-100 are
#       intensities for samples S1-S96. Missing values = not detected.
#   - 20220113_uniprot_annotation.RData
#       Object `anno_mapping`: a list of UniProt annotation tables
#       (human / mouse / rat). The rat table is used here.
#   - Matrisome_annotation_output.xlsx
#       Returned by the external MatrisomeAnalyzer web tool (see Step 3).
#   - Sample_index.xlsx
#       Columns Index, Sample.ID (S1-S96), Sample.group (16 groups, e.g. ABF).
#
# Outputs (written to this folder):
#   - Matrisome_annotation_input.xlsx
#       Rat UniProt annotation for the detected proteins; uploaded to
#       MatrisomeAnalyzer (Step 3).
#   - Protein_counts_summary.xlsx  (used by Script_2 and Script_3)
#       Sheets "Total_proteins", "Matrisome_proteins", "Core_matrisome",
#       "Matrisome_associated": 6 rows (replicates) x 16 columns (groups).
#       Counts are stored as numbers.
#       Sheet "Raw_table_annotation": 7 annotation columns + S1-S96 intensities.
#   - Matrisome_annotation_table.xlsx
#       Same content as the "Raw_table_annotation" sheet (single-sheet copy).
#
# Run in RStudio: setwd() and View() below require an interactive session.
# =============================================================================

# Load required packages.
# rio: import(); readxl: read_excel(); writexl: write_xlsx().
# tidyr and dplyr are loaded but not used in this script.
library(rio)
library(readxl)
library(tidyr)
library(writexl)
library(dplyr)


# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Step 1: Read the DIA-NN protein-group matrix ---------------------------------

data1 <- import("report.pg_matrix.tsv")   # 3,868 protein groups x 100 columns

# Use the first UniProt accession of each protein group as the row name
# (e.g. "P12345;Q67890" -> "P12345"). Only 8 groups contain more than one
# accession, and the first accessions are unique.
rownames(data1) <- gsub(";.*","",data1$Protein.Group) ## regex

# Drop the 4 descriptive columns; keep the 96 intensity columns S1-S96.
data1 <-  data1[,-(1:4)]

View(head(data1))


# Step 2: ID mapping and create annotation table (feature metadata) ------------

load("20220113_uniprot_annotation.RData")
anno_table <- anno_mapping[[match("rat",names(anno_mapping))]]

# Look up each detected accession in the rat UniProt table, in data1 row order.
anno_table <- anno_table[match(rownames(data1),anno_table$Uniprot.ID),]

# Drop accessions not found in the 2022 UniProt table.
# 3,859 of 3,868 protein groups are annotated. The 9 unmatched groups
# (Sdsl, Oxsr1, Micu3, Mbnl1, Rbpms, Stub1, Ppm1k, Fcgr3a, CYP2J4) are removed
# from all downstream tables, including the total protein counts.
anno_table <- subset(anno_table, !is.na(anno_table$Uniprot.ID))      #3859 protein IDs are annotated from 3868
rm(anno_mapping)

# Keep the same 3,859 protein groups in data1, in the same order as anno_table.
data1 <- data1[match(anno_table$Uniprot.ID, rownames(data1)),]

View(head(anno_table))

# Columns: Uniprot.ID, Uniprot.name, Gene.symbol, Description, Mol.Weight.
write_xlsx(
  anno_table,
  "Matrisome_annotation_input.xlsx"
)


# Step 3: External matrisome annotation (manual, outside R) --------------------
# Upload Matrisome_annotation_input.xlsx to MatrisomeAnalyzer:
#   https://izzilab-matrisomeanalyzer.share.connect.posit.cloud/
# Settings:
#   - Identifier column: gene symbol (Gene.symbol)
#   - Ortholog mapping: mouse
# Save the returned table as "Matrisome_annotation_output.xlsx".
# It contains the 5 input columns plus:
#   - "Annotated Matrisome Division": Core matrisome / Matrisome-associated /
#     Non-matrisome (60 / 125 / 3,674 protein groups)
#   - "Annotated Matrisome Category": e.g. Collagens, ECM Glycoproteins
# Gene symbols come from the 2022 UniProt table; symbols such as "-",
# "RGD:...", "ENSRNOG..." or "LOC..." cannot be matched by the tool and are
# returned as Non-matrisome.


# Step 4: Count detected protein groups per sample -----------------------------

# Import the externally annotated table.
table1 <- read_excel("Matrisome_annotation_output.xlsx")

# Check that the annotation rows are in the same order as the intensity rows.
all(table1$Uniprot.ID == rownames(data1)) # True

# Combine annotations (columns 1-7) and intensities (columns 8-103).
table1 <- cbind(table1, data1)


# Count non-missing intensities for each sample.
# Columns 1-7 contain annotations; columns 8-103 represent 96 samples (S1-S96).
# Each non-NA entry is counted as detected (the matrix contains no zero values).
j = 1
id <- c()

for (i in 8:103) {
  id[j] <- sum(!is.na(table1[[i]]))
  j = j + 1
}


# Subset protein groups by matrisome division.
# table2: all matrisome (Core matrisome + Matrisome-associated)
# table3: Core matrisome only
# table4: Matrisome-associated only
table2 <- table1[table1$`Annotated Matrisome Division` != "Non-matrisome", ]
table3 <- table2[table2$`Annotated Matrisome Division` == "Core matrisome", ]
table4 <- table2[table2$`Annotated Matrisome Division` == "Matrisome-associated", ]


# Count detected matrisome protein groups in each sample.
matrisome.id <- c()
j = 1

for (i in 8:103) {
  matrisome.id[j] <- sum(!is.na(table2[[i]]))
  j = j + 1
}


# Count detected core matrisome protein groups in each sample.
core.id <- c()
j = 1

for (i in 8:103) {
  core.id[j] <- sum(!is.na(table3[[i]]))
  j = j + 1
}


# Count detected matrisome-associated protein groups in each sample.
matrisome.asso.id <- c()
j = 1

for (i in 8:103) {
  matrisome.asso.id[j] <- sum(!is.na(table4[[i]]))
  j = j + 1
}


# Step 5: Attach experimental groups and reshape to group-wise tables ----------

# Combine the four counts into a sample-level summary (96 rows).
# Intensity columns are assumed to be ordered S1 through S96.
# data.frame() keeps the count columns numeric (cbind() would convert them to text).
summary.id <- data.frame("Sample.ID" = paste0("S", 1:96), id, matrisome.id, core.id, matrisome.asso.id)
sample.id <- as.data.frame(paste0("S", 1:96))   # not used below


# Read sample identifiers and experimental group assignments.
# Sample.group codes: 1st letter = age (A/Y), 2nd = lung region (B/T/C/P),
# 3rd = sex (F/M); 6 samples per group.
sample.table <- read_excel("Sample_index.xlsx")


# Reorder samples to match Sample_index.xlsx, then attach group labels.
# Columns: Sample.ID, id, matrisome.id, core.id, matrisome.asso.id, Sample.Group
summary.id <- summary.id[match(sample.table$Sample.ID, summary.id$Sample.ID), ]
summary.id <- cbind(summary.id, sample.table$Sample.group)
colnames(summary.id)[6] <- "Sample.Group"


# Create a separate count-and-group table for each protein category.
id.table <- summary.id[c(2, 6)]                  # total proteins
matrisome.id.table <- summary.id[c(3, 6)]        # all matrisome
core.id.table <- summary.id[c(4, 6)]             # core matrisome
matrisome.asso.id.table <- summary.id[c(5, 6)]   # matrisome-associated


# Arrange total protein counts with one column per experimental group.
# The study design assumes 16 groups with six samples per group.
# Group columns follow their first appearance in Sample_index.xlsx.
id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))  # Initialize six sample rows.

# Preserve sample-index order within each group when populating columns.
for (group in unique(id.table$Sample.Group)) {
  group_values <- id.table$id[id.table$Sample.Group == group]
  id.table.wide[[group]] <- group_values
}


# Arrange all matrisome counts in the same group-wise format.
matrisome.id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))

for (group in unique(matrisome.id.table$Sample.Group)) {
  group_values <- matrisome.id.table$matrisome.id[matrisome.id.table$Sample.Group == group]
  matrisome.id.table.wide[[group]] <- group_values
}


# Arrange core matrisome counts in the same group-wise format.
core.id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))

for (group in unique(core.id.table$Sample.Group)) {
  group_values <- core.id.table$core.id[core.id.table$Sample.Group == group]
  core.id.table.wide[[group]] <- group_values
}


# Arrange matrisome-associated counts in the same group-wise format.
matrisome.asso.id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))

for (group in unique(matrisome.asso.id.table$Sample.Group)) {
  group_values <- matrisome.asso.id.table$matrisome.asso.id[matrisome.asso.id.table$Sample.Group == group]
  matrisome.asso.id.table.wide[[group]] <- group_values
}


# Step 6: Prepare the annotated intensity table --------------------------------

# Retain input intensities and annotations, with standardized category names.
# Column renaming (Script_2 and Script_3 use these new names):
#   "Annotated Matrisome Division" -> Matrisome.Category
#     (Core matrisome / Matrisome-associated)
#   "Annotated Matrisome Category" -> Matrisome.Sub.category
#     (e.g. Collagens, ECM Glycoproteins)
raw_table_annotation <- table1
colnames(raw_table_annotation)[6] <- 'Matrisome.Category'
colnames(raw_table_annotation)[7] <- 'Matrisome.Sub.category'


# Clear category and subcategory labels for non-matrisome entries in the export
# (Non-matrisome rows get empty strings in both columns).
# This loop assumes the annotation table contains 3,859 protein-group rows.
for (i in 1:3859) {
  if (raw_table_annotation$Matrisome.Category[i] == "Non-matrisome") {
    raw_table_annotation$Matrisome.Category[i] <- ""
    raw_table_annotation$Matrisome.Sub.category[i] <- ""
  }
}


# Step 7: Write output workbooks -----------------------------------------------

# Collect count tables and the annotated intensity table as separate worksheets.
# The four count sheets are used for the protein-count bar plots.

tables_list <- list(
  "Total_proteins" = id.table.wide,
  "Matrisome_proteins" = matrisome.id.table.wide,
  "Core_matrisome" = core.id.table.wide,
  "Matrisome_associated" = matrisome.asso.id.table.wide,
  "Raw_table_annotation" = raw_table_annotation
)


# Single-sheet copy of the annotated intensity table.
write_xlsx(raw_table_annotation, "Matrisome_annotation_table.xlsx")
# Save the data table for Script_2 and _3
# (both read the "Raw_table_annotation" sheet).
write_xlsx(tables_list, "Protein_counts_summary.xlsx")
