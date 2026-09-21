# Load required packages.
library(readxl)
library(tidyr)
library(rio)


# Summarize detected peptide entries per sample using protein-level matrisome annotations.
# Inputs: peptide.xlsx, Protein_counts_summary.xlsx, and Sample_metadata.xlsx.
# Output: Peptide_counts_summary.xlsx.


# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Read the peptide-level input table from the first worksheet.
peptide <- read_excel("Peptide.xlsx")


# Insert empty matrisome annotation columns after the four identifier columns.
peptide <- cbind(
  peptide[, 1:4, drop = FALSE],
  Matrisome.Category = "",
  Matrisome.Sub.Category = "",
  peptide[, 5:ncol(peptide), drop = FALSE]
)


# Read the protein-level annotations used to classify peptide entries.
raw_table <- read_excel("Protein_counts_summary.xlsx",
  sheet = "Raw_table_annotation")
colnames(raw_table)[7] <- "Matrisome.Sub.Category"


# Select protein identifiers, matrisome divisions, and matrisome categories.
raw_table1 <- raw_table[c(1, 6, 7)]
protein_list <- raw_table1$Uniprot.ID


# Transfer annotations by exact matching of Protein.Group to Uniprot.ID.
# The loop assumes 53,032 peptide entries and one reference match per identifier.
# Unmatched entries retain empty strings; grouped identifiers are not split.
i = 1

for (i in 1:53032) {
  if (peptide$Protein.Group[i] %in% protein_list) {
    j <- which(protein_list %in% peptide$Protein.Group[i])
    peptide$Matrisome.Category[i] <- raw_table1$Matrisome.Category[j]
    peptide$Matrisome.Sub.Category[i] <- raw_table1$Matrisome.Sub.Category[j]
  }
}
table1 <- peptide


# Count detected peptide entries per sample using non-missing intensities.
# Columns 7-102 are assumed to contain the 96 sample intensity columns.
# Each non-NA value, including zero, contributes one count; rows are not deduplicated.
j = 1
id <- c()

for (i in 7:102) {
  id[j] <- sum(!is.na(table1[[i]]))
  j = j + 1
}


# Select entries with non-NA matrisome labels, then extract the two divisions.
# The current selection retains empty strings as well as assigned labels.
table2 <- table1[!is.na(table1$Matrisome.Category), ]
table3 <- table2[table2$Matrisome.Category == "Core matrisome", ]
table4 <- table2[table2$Matrisome.Category == "Matrisome-associated", ]


# Count detected entries in the non-NA annotation subset for each sample.
matrisome.id <- c()
j = 1

for (i in 7:102) {
  matrisome.id[j] <- sum(!is.na(table2[[i]]))
  j = j + 1
}


# Count detected core matrisome peptide entries for each sample.
core.id <- c()
j = 1

for (i in 7:102) {
  core.id[j] <- sum(!is.na(table3[[i]]))
  j = j + 1
}


# Count detected matrisome-associated peptide entries for each sample.
matrisome.asso.id <- c()
j = 1

for (i in 7:102) {
  matrisome.asso.id[j] <- sum(!is.na(table4[[i]]))
  j = j + 1
}


# Combine the four count vectors into a sample-level summary.
# Sample columns are assumed to be ordered S1 through S96.
summary.id <- as.data.frame(cbind("Sample.ID" = paste0("S", 1:96), id, matrisome.id, core.id, matrisome.asso.id))
sample.id <- as.data.frame(paste0("S", 1:96))


# Read experimental group labels from the sample metadata.
sample.table <- read_excel("Sample_metadata.xlsx")


# Attach group labels by row position.
# Metadata rows are assumed to follow the same S1-S96 sample order.
summary.id <- cbind(summary.id, sample.table$Sample.group)
colnames(summary.id)[6] <- "Sample.Group"


# Create a separate count-and-group table for each annotation subset.
id.table <- summary.id[c(2, 6)]
matrisome.id.table <- summary.id[c(3, 6)]
core.id.table <- summary.id[c(4, 6)]
matrisome.asso.id.table <- summary.id[c(5, 6)]


# Arrange total peptide counts with one column per experimental group.
# The study design assumes 16 groups with six samples per group.
id.table.wide <- data.frame(matrix(nrow = 6, ncol = 0))

for (group in unique(id.table$Sample.Group)) {
  group_values <- id.table$id[id.table$Sample.Group == group]
  id.table.wide[[group]] <- group_values
}


# Arrange counts from the non-NA annotation subset in the same group-wise format.
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


# Collect the four count tables as separate Excel worksheets.
# The existing worksheet names containing "proteins" are retained for compatibility;
# the values in this script count peptide-table entries.
tables_list <- list(
  "Total_proteins" = id.table.wide,
  "Matrisome_proteins" = matrisome.id.table.wide,
  "Core_matrisome" = core.id.table.wide,
  "Matrisome_associated" = matrisome.asso.id.table.wide
)


# Export the peptide count summary for downstream plotting and inspection.
export(tables_list, "Peptide_counts_summary.xlsx")
