# Load required packages.
library(dplyr)
library(rio)
library(readxl)
library(NormalyzerDE)
library(missForest)
library(imputeLCMD)
library(ggplot2)
library(ggbiplot)

# Preprocess protein-group intensities and generate quality-control plots.
# Inputs: DIA-NN intensity matrix, sample metadata, and matrisome annotations.
# Outputs: QC/PCA PDFs, detection counts, intensity workbooks, and RData files.
# Configure the study and output log.
# Set the working directory to this script's folder in RStudio.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

study_number <- "UW.Madison.rat.lung.proteomics"
data_file <- "report.pg_matrix.tsv"
search_method <- "DIANN"
sample_info <- "Sample_metadata.xlsx"
species <- "rat"
sample_order <- c("YTM", "YCM", "YPM", "YBM", "YTF", "YCF", "YPF", "YBF", "ATM", "ACM", "APM", "ABM", "ATF", "ACF", "APF", "ABF")
# to save space and keep consistency, Y means peripubertal, while A means young adult in all the following analysis

log_file <- paste0(study_number, " analysis log.txt")
cat(paste("Study number:", study_number), file = log_file, sep = "\n")

# load data1
data1 <- readRDS("expression_matrix.rds")

# Read the seven annotation columns from the protein-count summary workbook.
matrisome_anno_table <- read_excel("Protein_counts_summary.xlsx",
  sheet = "Raw_table_annotation")
matrisome_anno_table <- matrisome_anno_table[c(1:7)]


# Read the first sample-metadata worksheet and standardize its column names.
info_table <- import(sample_info)

colnames(info_table) <- c("Sample.ID", "Sample.Group", "Sample.PID", "Sample.Animal.Number", "Sample.Age", "Sample.Region", "Sample.Sex")

# Define the experimental-group order and encode age, region, and sex as factors.
info_table$Sample.Group <- factor(info_table$Sample.Group, levels = sample_order)
info_table$Sample.Age <- factor(info_table$Sample.Age)
info_table$Sample.Region <- factor(info_table$Sample.Region)
info_table$Sample.Sex <- factor(info_table$Sample.Sex)

cat(paste0(Sys.time(), " Total protein groups identified: ", nrow(data1), "\n"), file = log_file, append = TRUE)


# Build a sample-by-protein detection matrix: non-missing = 1, missing = 0.
# The loops assume 96 samples and 3,859 protein groups.
work1 <- as.data.frame(t(data1))

i = 1
j = 1

for (i in 1:96) {
  for (j in 1:3859) {
    if (!is.na(work1[i, j])) {work1[i, j] <- 1}
    else {work1[i, j] <- 0}
  }
}


# Apply the global detection criterion: at least 48 of 96 samples (50%).
work2 <- colSums(work1)
col_indices <- which(work2 >= 48)
work3 <- work2[work2 >= 48]


# Count observed values within each of the 16 experimental groups.
# Grouping assumes metadata rows are aligned with the samples in work1.
sample_groups <- split(work1, info_table$Sample.Group)
sample_group_sums <- lapply(sample_groups, colSums)

sample_group_number <- as.data.frame(cbind(sample_group_sums[[1]], sample_group_sums[[2]], sample_group_sums[[3]], sample_group_sums[[4]],
  sample_group_sums[[5]], sample_group_sums[[6]], sample_group_sums[[7]], sample_group_sums[[8]],
  sample_group_sums[[9]], sample_group_sums[[10]], sample_group_sums[[11]], sample_group_sums[[12]],
  sample_group_sums[[13]], sample_group_sums[[14]], sample_group_sums[[15]], sample_group_sums[[16]]

))
sample_group_number <- cbind(rownames(sample_group_number), sample_group_number)

colnames(sample_group_number) <- c("Uniprot", "YTM", "YCM", "YPM", "YBM",
  "YTF", "YCF", "YPF", "YBF",
  "ATM", "ACM", "APM", "ABM",
  "ATF", "ACF", "APF", "ABF"

)


# Export the group-level detection counts with protein identifiers.
export(sample_group_number, "Real_values_number_in_sample_group.xlsx")

sample_group_number <- sample_group_number[-1]


# Apply the conditional criterion: at least five of six samples in any group.
common_id <- c()
j = 1
i = 1

for (i in 1:3859) {
  if (any(sample_group_number[i, ] >= 5)) {
    common_id <- c(common_id, rownames(sample_group_number)[i])
  }
}


# Retain proteins meeting either the global or conditional detection criterion.
common_id1 <- union(common_id, names(work3))

data2 <- data1[rownames(data1) %in% common_id1, ] # 3600 proteins after filter


# Record retained protein counts and missing-value fractions.
cat(paste0(Sys.time(), " Protein groups passed missing value filtering: ", nrow(data2), "\n"), file = log_file, append = TRUE)
cat(paste0(Sys.time(), " Rate of missing values before filtering: ", sum(is.na(data1)) / ncol(data1) / nrow(data1), "\n"), file = log_file, append = TRUE)
cat(paste0(Sys.time(), " Rate of missing values after filtering: ", sum(is.na(data2)) / ncol(data2) / nrow(data2), "\n"), file = log_file, append = TRUE)

View(head(data2))


# Transform absolute intensities to log2 values, then apply median normalization.
# The data2 object is overwritten with log2-transformed filtered intensities.
data2 <- log2(abs(data2))

data3 <- medianNormalization(as.matrix(data2), noLogTransform = T)
data3 <- as.data.frame(data3)
rownames(data3) <- rownames(data2)

View(head(data3))


# Summarize missingness per protein and detection counts per sample for QC.
na.count1 <- rowSums(is.na(data1))
pro.count1 <- colSums(!is.na(data1))
na.count2 <- rowSums(is.na(data2))
pro.count2 <- colSums(!is.na(data2))


# Save QC plots showing missingness, detection counts, and intensity distributions
# before and after filtering or normalization.
pdf(paste0(study_number, "_QC plots ", Sys.Date(), ".pdf"))

hist(na.count1, freq = F, xlab = "NA count", ylab = "% protein id", main = "Missing value/protein before filtering", col = 4, xlim = c(0, ncol(data1)))
abline(v = ncol(data1) * 0.5, col = "red", lty = 3, lwd = 2)

barplot(pro.count1, las = 2, xlab = "", ylab = "# of protein identified", main = "protein IDs in each sample before filtering", col = 4)

hist(na.count2, freq = F, xlab = "NA count", ylab = "% protein id", main = "Missing value/protein after filtering", col = 4, xlim = c(0, ncol(data1)))
abline(v = ncol(data1) * 0.5, col = "red", lty = 3, lwd = 2)

barplot(pro.count2, las = 2, xlab = "", ylab = "# of protein identified", main = "protein IDs in each sample after filtering", col = 4)

boxplot(data2, col = 'red', boxlwd = 1.5, medlwd = 1.5, whisklwd = 1.5, staplelwd = 1.5, outlwd = 1.5, cex.axis = 0.6, font.axis = 2, xlab = "Sample", ylab = "Log2 protein intensity", main = "Boxplot before normalization")

boxplot(data3, col = 'red', boxlwd = 1.5, medlwd = 1.5, whisklwd = 1.5, staplelwd = 1.5, outlwd = 1.5, cex.axis = 0.6, font.axis = 2, xlab = "Sample", ylab = "Log2 protein intensity", main = "Boxplot after normalization")

dev.off()


# Imputation stage 1: apply random forest imputation to the normalized matrix.
# Extract the completed matrix and record the out-of-bag error.
data4 <- missForest(data3)
data5 <- data4$ximp

cat(paste0(Sys.time(), " Imputation has been performed using Random Forest, OOB error of RF = ", data4$OOBerror, "\n"), file = log_file, append = TRUE)

View(head(data5))


# Build a mask for proteins with no observed values in an experimental group.
# The mask assumes 3,859 protein rows and sample columns S1 through S96.
work4 <- matrix(1, nrow = 3859, ncol = 96)
colnames(work4) <- paste0("S", 1:96)

for (group in info_table$Sample.Group) {
  zero_proteins <- which(sample_group_number[[group]] == 0)
  group_samples <- info_table$Sample.ID[info_table$Sample.Group == group]
  work4[zero_proteins, group_samples] <- 0
}


# Subset the mask to retained proteins and inspect row alignment.
work4 <- as.data.frame(work4)

rownames(work4) <- rownames(data1)
work5 <- work4[rownames(work4) %in% rownames(data5), ]
all(rownames(work5) == rownames(data5))


# Apply the mask to the random forest estimates and convert zero entries to NA.
# This restores missingness for group-wise absence before the second imputation.
data6 <- data5 * work5
data6[data6 == 0] <- NA


# Imputation stage 2: impute remaining missing values using QRILC
# (quantile regression imputation of left-censored data), with tune.sigma = 0.5.
data7 <- impute.QRILC(data6, tune.sigma = 0.5)[[1]]


# Count entries changed by the second imputation for each protein and sample.
data5_7 <- data5 - data7

row_zeros <- apply(data5_7, 1, function(x) sum(x != 0, na.rm = TRUE))
col_zeros <- apply(data5_7, 2, function(x) sum(x != 0, na.rm = TRUE))
col_zeros

# Random forest and QRILC imputation were performed without setting a random seed
# in the original analysis. Re-running these steps may therefore produce
# different imputed values and downstream results.
# The imputed matrix used in the reported analyses is provided as an intermediate
# data file to enable reproduction of downstream analyses from that matrix.
cat(paste0(Sys.time(), " Imputation has been performed using QRILC on the basis of RF, protein IDs changed:", sum(col_zeros), "\n"), file = log_file, append = TRUE)


# Load the external ggbiplot helper required for the PCA visualizations.
source("T:/Proteomics_Data_Analysis/SOP_Liang/Proteomics data analysis_R template/ggbiplot.r")


# Run PCA with samples as observations and save plots colored by sample group,


# Display PCA by age, with point shapes indicating sex.  Figure 2a
g2 <- ggbiplot(prcomp(t(data7)), var.axes = F, labels = NULL, labels.size = 3, groups = info_table$Sample.Age, obs.scale = 1, var.scale = 1, ellipse = T, circle = T)
g2 <- g2 + scale_shape_discrete(solid = F) + geom_point(aes(shape = info_table$Sample.Sex))
g2 <- g2 + theme_bw(base_size = 16)

plot(g2)


# Display PCA by sex, with point shapes indicating region.  Figure 5a
g3 <- ggbiplot(prcomp(t(data7)), var.axes = F, labels = NULL, labels.size = 3, groups = info_table$Sample.Sex, obs.scale = 1, var.scale = 1, ellipse = T, circle = T)
g3 <- g3 + scale_shape_discrete(solid = F) + geom_point(aes(shape = info_table$Sample.Region))
g3 <- g3 + theme_bw(base_size = 16)

plot(g3)



# Display PCA by region, with point shapes indicating age.  Figure 6a
g5 <- ggbiplot(prcomp(t(data7)), var.axes = F, labels = NULL, labels.size = 3, groups = info_table$Sample.Region, obs.scale = 1, var.scale = 1, ellipse = T, circle = T)
g5 <- g5 + scale_shape_discrete(solid = F) + geom_point(aes(shape = info_table$Sample.Age))
g5 <- g5 + theme_bw(base_size = 16)

plot(g5)


# Export the selected regional, age, and sex PCA views to a combined PDF.
pdf("PCA.pdf")
plot(g2)
plot(g3)
plot(g5)
dev.off()


# Collect original, filtered/log2-transformed, normalized, and imputed matrices.
data_table <- list("data_ori" = data1,
  "data_filtered" = data2,
  "data_nor" = data3,
  "data_imp" = data7
)


# Save the analysis matrices, annotation object, and sample metadata.
save(data_table, matrisome_anno_table, info_table, file = paste0(study_number, "_data tables_", Sys.Date(), ".RData"))


# Inspect annotation ordering and protein-row consistency across matrices.
# These expressions report logical results without stopping execution.
matrisome_anno_table1 <- matrisome_anno_table[order(matrisome_anno_table$Uniprot.ID), ]
all(matrisome_anno_table1[[1]] == matrisome_anno_table[[1]])
remove(matrisome_anno_table1)
all(rownames(data1) == matrisome_anno_table$Uniprot.ID)
all(rownames(data2) == rownames(data3))
all(rownames(data2) == rownames(data5))
all(rownames(data2) == rownames(data7))


# Match annotations to protein identifiers separately for each matrix
# and prepare the four annotated worksheets for export.
data_table_export <- list("data_ori" = cbind(matrisome_anno_table[match(rownames(data1), matrisome_anno_table$Uniprot.ID), ], data1),
  "data_filtered" = cbind(matrisome_anno_table[match(rownames(data2), matrisome_anno_table$Uniprot.ID), ], data2),
  "data_nor" = cbind(matrisome_anno_table[match(rownames(data3), matrisome_anno_table$Uniprot.ID), ], data3),
  "data_imp" = cbind(matrisome_anno_table[match(rownames(data7), matrisome_anno_table$Uniprot.ID), ], data7)
)


# Export annotated intensity tables and record completion in the log.
export(data_table_export, file = paste0(study_number, "_", search_method, "_protein intensity_", Sys.Date(), ".xlsx"))
cat(paste0(Sys.time(), " Excel file of protein intensity data tables is created", "\n"), file = log_file, append = TRUE)


# Save the current R workspace, including intermediate objects.
save.image("Script_3_data_preprocessing_data_summary.RData")
