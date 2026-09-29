# Load required packages.
library(readxl)
library(tidyr)
library(writexl)
library(dplyr)

# Summarize detected protein groups per sample by matrisome division.
# Inputs: Matrisome_annotation_table.xlsx and Sample_index.xlsx.
# Outputs: an Excel summary workbook and an R workspace.
# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


data1 <- import("report.pg_matrix.tsv")
rownames(data1) <- gsub(";.*","",data1$Protein.Group) ## regex
data1 <-  data1[,-(1:4)]

View(head(data1))


# ID mapping and create annatation table (feature metadata) ----------------------------------

load("20220113_uniprot_annotation.RData")
anno_table <- anno_mapping[[match("rat",names(anno_mapping))]]
anno_table <- anno_table[match(rownames(data1),anno_table$Uniprot.ID),]
anno_table <- subset(anno_table, !is.na(anno_table$Uniprot.ID))      #3859 protein IDs are annotated from 3868
rm(anno_mapping)
data1 <- data1[match(anno_table$Uniprot.ID, rownames(data1)),] 

View(head(anno_table))

write_xlsx(
  anno_table,
  "Matrisome_annotation_input.xlsx"
)


# External annotation step:
# Matrisome annotations are from https://izzilab-matrisomeanalyzer.share.connect.posit.cloud/
# Identifier column is gene symbol
# Ortholog mapping using mouse species
# Save the returned table as "Matrisome_annotation_output.xlsx"


# Import the externally annotated table.
table1 <- read_excel("Matrisome_annotation_output.xlsx")
all(table1$Uniprot.ID == rownames(data1)) # True
table1 <- cbind(table1, data1)


# Count non-missing intensities for each sample.
# Columns 1-7 contain annotations; columns 8-103 represent 96 samples.
# Each non-NA entry is counted as detected, including zero values.
j = 1
id <- c()

for (i in 8:103) {
  id[j] <- sum(!is.na(table1[[i]]))
  j = j + 1
}


# Separate matrisome entries into core and matrisome-associated divisions.
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


# Combine the four counts into a sample-level summary.
# Intensity columns are assumed to be ordered S1 through S96.
summary.id <- as.data.frame(cbind("Sample.ID" = paste0("S", 1:96), id, matrisome.id, core.id, matrisome.asso.id))
sample.id <- as.data.frame(paste0("S", 1:96))


# Read sample identifiers and experimental group assignments.
sample.table <- read_excel("Sample_index.xlsx")


# Match the sample-index order before attaching experimental group labels.
summary.id <- summary.id[match(sample.table$Sample.ID, summary.id$Sample.ID), ]
summary.id <- cbind(summary.id, sample.table$Sample.group)
colnames(summary.id)[6] <- "Sample.Group"


# Create a separate count-and-group table for each protein category.
id.table <- summary.id[c(2, 6)]
matrisome.id.table <- summary.id[c(3, 6)]
core.id.table <- summary.id[c(4, 6)]
matrisome.asso.id.table <- summary.id[c(5, 6)]


# Arrange total protein counts with one column per experimental group.
# The study design assumes 16 groups with six samples per group.
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


# Retain input intensities and annotations, with standardized category names.
raw_table_annotation <- table1
colnames(raw_table_annotation)[6] <- 'Matrisome.Category'
colnames(raw_table_annotation)[7] <- 'Matrisome.Sub.category'


# Clear category and subcategory labels for non-matrisome entries in the export.
# This loop assumes the annotation table contains 3,859 protein-group rows.
for (i in 1:3859) {
  if (raw_table_annotation$Matrisome.Category[i] == "Non-matrisome") {
    raw_table_annotation$Matrisome.Category[i] <- ""
    raw_table_annotation$Matrisome.Sub.category[i] <- ""
  }
}


# Collect count tables and the annotated intensity table as separate worksheets.
# make bar plot for figure xxx

tables_list <- list(
  "Total_proteins" = id.table.wide,
  "Matrisome_proteins" = matrisome.id.table.wide,
  "Core_matrisome" = core.id.table.wide,
  "Matrisome_associated" = matrisome.asso.id.table.wide,
  "Raw_table_annotation" = raw_table_annotation
)


write_xlsx(raw_table_annotation, "Matrisome_annotation_table.xlsx")
# Save the data table for Script_2 and _3
write_xlsx(tables_list, "Protein_counts_summary.xlsx")

