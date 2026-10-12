# =============================================================================
# Script 2: Peptide counts per sample, split by matrisome division
# =============================================================================
#
# Purpose:
#   Count the number of detected peptides in each of the 96 lung samples,
#   for (1) all peptides, (2) all matrisome peptides, (3) core matrisome, and
#   (4) matrisome-associated peptides, then arrange the counts by experimental
#   group (16 groups x 6 samples). Each peptide inherits the matrisome
#   annotation of its protein group from Script_1.
#
# Inputs (all in this folder):
#   - Peptide.xlsx
#       DIA-NN peptide matrix, 53,032 rows (one row per peptide) x 100 columns.
#       Columns 1-4 are Protein.Group, Protein.Ids, Protein.Names, Genes;
#       columns 5-100 are intensities named by raw file path, in run order
#       AGT_018_1 ... AGT_018_96 (= samples S1-S96). Missing values = not
#       detected.
#   - Protein_counts_summary.xlsx  (output of Script_1)
#       Only sheet "Raw_table_annotation" is used: Uniprot.ID (column 1),
#       Matrisome.Category (column 6; Core matrisome / Matrisome-associated,
#       blank for Non-matrisome) and Matrisome.Sub.category (column 7).
#   - Sample_metadata.xlsx
#       96 rows in S1-S96 order. Columns Sample.ID, Sample.group, PID,
#       Investigator Animal Number, Age, Region, Sex.
#
# Output (written to this folder):
#   - Peptide_counts_summary.xlsx
#       Sheets "Total_proteins", "Matrisome_proteins", "Core_matrisome",
#       "Matrisome_associated": 6 rows (replicates) x 16 columns (groups).
#       Counts are stored as numbers. Despite the sheet names, all values
#       are peptide counts (names kept to match Script_1's layout).
#       Not read by any downstream script; used for the peptide-count plots.
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# =============================================================================

# Load required packages.
# readxl: read_excel(); rio: export().
# tidyr is loaded but not used in this script.
library(readxl)
library(tidyr)
library(rio)


# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Step 1: Read the DIA-NN peptide matrix ---------------------------------------

# 53,032 peptides x 100 columns (first worksheet).
peptide <- read_excel("Peptide.xlsx")


# Insert two empty annotation columns after the four identifier columns.
# New layout: columns 1-4 identifiers, 5-6 annotations, 7-102 samples S1-S96.
peptide <- cbind(
  peptide[, 1:4, drop = FALSE],
  Matrisome.Category = "",
  Matrisome.Sub.Category = "",
  peptide[, 5:ncol(peptide), drop = FALSE]
)


# Step 2: Transfer protein-level matrisome annotation to peptides --------------

# Read the protein-level annotation table produced by Script_1.
# read_excel() reads the blank labels of Non-matrisome proteins as NA.
raw_table <- read_excel("Protein_counts_summary.xlsx",
  sheet = "Raw_table_annotation")
colnames(raw_table)[7] <- "Matrisome.Sub.Category"   # capital "C" to match peptide


# Keep Uniprot.ID, Matrisome.Category (division) and Matrisome.Sub.Category.
raw_table1 <- raw_table[c(1, 6, 7)]
protein_list <- raw_table1$Uniprot.ID


# Transfer annotations by exact matching of Protein.Group to Uniprot.ID.
# The loop assumes 53,032 peptide rows; each Uniprot.ID is unique.
# Matching results (peptide rows):
#   - Core matrisome:        1,358
#   - Matrisome-associated:  1,651
#   - Non-matrisome:        49,565  -> labels become NA
#   - Not found:               458  -> labels stay "" (159 protein groups,
#     e.g. multi-accession groups such as "P0DP29;P0DP30;P0DP31", which are
#     not split, and the 9 proteins dropped in Script_1). Matching by the
#     first accession would add only Non-matrisome peptides.
i = 1   # not needed; reset by the for loop

for (i in 1:53032) {
  if (peptide$Protein.Group[i] %in% protein_list) {
    j <- which(protein_list %in% peptide$Protein.Group[i])
    peptide$Matrisome.Category[i] <- raw_table1$Matrisome.Category[j]
    peptide$Matrisome.Sub.Category[i] <- raw_table1$Matrisome.Sub.Category[j]
  }
}
table1 <- peptide


# Step 3: Count detected peptides per sample -----------------------------------

# Count non-missing intensities for each sample.
# Columns 1-6 contain identifiers and annotations; columns 7-102 represent
# 96 samples (S1-S96).
# Each row is one peptide; each non-NA value counts as one detected peptide
# (the matrix contains no zero values).
j = 1
id <- c()

for (i in 7:102) {
  id[j] <- sum(!is.na(table1[[i]]))
  j = j + 1
}


# Subset peptides by matrisome division.
# table2: all matrisome (Core matrisome + Matrisome-associated)
# table3: Core matrisome only
# table4: Matrisome-associated only
# Both NA (Non-matrisome) and "" (not found) are excluded from table2.
table2 <- table1[!is.na(table1$Matrisome.Category) & table1$Matrisome.Category != "", ]
table3 <- table2[table2$Matrisome.Category == "Core matrisome", ]
table4 <- table2[table2$Matrisome.Category == "Matrisome-associated", ]


# Count detected matrisome peptides in each sample.
matrisome.id <- c()
j = 1

for (i in 7:102) {
  matrisome.id[j] <- sum(!is.na(table2[[i]]))
  j = j + 1
}


# Count detected core matrisome peptides in each sample.
core.id <- c()
j = 1

for (i in 7:102) {
  core.id[j] <- sum(!is.na(table3[[i]]))
  j = j + 1
}


# Count detected matrisome-associated peptides in each sample.
matrisome.asso.id <- c()
j = 1

for (i in 7:102) {
  matrisome.asso.id[j] <- sum(!is.na(table4[[i]]))
  j = j + 1
}


# Step 4: Attach experimental groups and reshape to group-wise tables ----------

# Combine the four counts into a sample-level summary (96 rows).
# Sample columns are ordered by raw file AGT_018_1 ... AGT_018_96, i.e. S1-S96.
# data.frame() keeps the count columns numeric (cbind() would convert them to text).
summary.id <- data.frame("Sample.ID" = paste0("S", 1:96), id, matrisome.id, core.id, matrisome.asso.id)
sample.id <- as.data.frame(paste0("S", 1:96))   # not used below


# Read sample identifiers and experimental group assignments.
# Sample.group codes: 1st letter = age (A = Aged, Y = Young),
# 2nd = lung region (B = Bottom, T = Top, C = Central, P = Peripheral),
# 3rd = sex (F/M); 6 samples per group.
sample.table <- read_excel("Sample_metadata.xlsx")


# Attach group labels by row position (no matching by Sample.ID).
# Sample_metadata.xlsx rows are in the same S1-S96 order as summary.id.
# Columns: Sample.ID, id, matrisome.id, core.id, matrisome.asso.id, Sample.Group
summary.id <- cbind(summary.id, sample.table$Sample.group)
colnames(summary.id)[6] <- "Sample.Group"


# Create a separate count-and-group table for each peptide category.
id.table <- summary.id[c(2, 6)]                  # total peptides
matrisome.id.table <- summary.id[c(3, 6)]        # all matrisome
core.id.table <- summary.id[c(4, 6)]             # core matrisome
matrisome.asso.id.table <- summary.id[c(5, 6)]   # matrisome-associated


# Arrange total peptide counts with one column per experimental group.
# The study design assumes 16 groups with six samples per group.
# Group columns follow their first appearance in Sample_metadata.xlsx
# (YTM, YCM, YPM, YBM, ATF, ...), which differs from the column order in
# Script_1's Protein_counts_summary.xlsx.
id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))  # Initialize six sample rows.

# Preserve S1-S96 order within each group when populating columns.
for (group in unique(id.table$Sample.Group)) {
  group_values <- id.table$id[id.table$Sample.Group == group]
  id.table.wide[[group]] <- group_values
}


# Arrange all matrisome peptide counts in the same group-wise format.
matrisome.id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))

for (group in unique(matrisome.id.table$Sample.Group)) {
  group_values <- matrisome.id.table$matrisome.id[matrisome.id.table$Sample.Group == group]
  matrisome.id.table.wide[[group]] <- group_values
}


# Arrange core matrisome peptide counts in the same group-wise format.
core.id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))

for (group in unique(core.id.table$Sample.Group)) {
  group_values <- core.id.table$core.id[core.id.table$Sample.Group == group]
  core.id.table.wide[[group]] <- group_values
}


# Arrange matrisome-associated peptide counts in the same group-wise format.
matrisome.asso.id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))

for (group in unique(matrisome.asso.id.table$Sample.Group)) {
  group_values <- matrisome.asso.id.table$matrisome.asso.id[matrisome.asso.id.table$Sample.Group == group]
  matrisome.asso.id.table.wide[[group]] <- group_values
}


# Step 5: Write output workbook ------------------------------------------------

# Collect the four count tables as separate worksheets.
# the values in this script are peptide counts.
tables_list <- list(
  "Total_proteins" = id.table.wide,
  "Matrisome_proteins" = matrisome.id.table.wide,
  "Core_matrisome" = core.id.table.wide,
  "Matrisome_associated" = matrisome.asso.id.table.wide
)


# Export the peptide count summary for plotting the peptide-count bar plots in figure S2.
export(tables_list, "Peptide_counts_summary.xlsx")
