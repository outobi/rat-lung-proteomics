# =============================================================================
# Script 3: Protein intensity preprocessing (filtering, normalization,
#           imputation) and QC / PCA plots
# =============================================================================
#
# Purpose:
#   Starting from the annotated protein-group intensity table of Script_1,
#   (1) filter protein groups by detection rate, (2) log2-transform and
#   median-normalize, (3) impute missing values in two stages (random forest,
#   then QRILC for group-wise absence), and (4) draw QC and PCA plots.
#   The imputed matrix is the input for all differential-expression scripts.
#
# Inputs (all in this folder):
#   - Protein_counts_summary.xlsx  (output of Script_1)
#       Sheet "Raw_table_annotation": columns 1-7 are annotations
#       (Uniprot.ID, Uniprot.name, Gene.symbol, Description, Mol.Weight,
#       Matrisome.Category, Matrisome.Sub.category); columns 8-103 are raw
#       intensities for samples S1-S96. 3,859 protein groups.
#   - Sample_metadata.xlsx
#       96 rows in S1-S96 order. Columns Sample.ID, Sample.group, PID,
#       Investigator Animal Number, Age, Region, Sex.
#       Age labels "Young" / "Aged" denote 7-week-old (7w) / 7-11-month-old (7-11m) rats.
#
# Outputs (written to this folder):
#   - UW.Madison.rat.lung.proteomics analysis log.txt
#       Run log (overwritten on each run).
#   - Real_values_number_in_sample_group.xlsx  (used by Script_8)
#       3,859 protein groups x 16 groups: number of samples (0-6) with an
#       observed intensity in each group.
#   - UW.Madison.rat.lung.proteomics_QC plots.pdf
#       Missing-value histograms, per-sample protein counts, and intensity
#       boxplots before/after filtering and normalization.
#   - PCA.pdf
#       PCA of the imputed matrix (Figures 2a, 5a, 6a).
#   - UW.Madison.rat.lung.proteomics_data_tables.RData
#       (used by Script_4, Script_6 and Script_8)
#       data_table (list of data_ori, data_filtered, data_nor, data_imp),
#       matrisome_anno_table, info_table.
#   - UW.Madison.rat.lung.proteomics_DIANN_protein_intensity.xlsx
#       (used by Script_4, Script_6 and Script_8)
#       Sheets data_ori, data_filtered, data_nor, data_imp, each with the
#       7 annotation columns followed by S1-S96.
#   - Script_3_data_preprocessing_data_summary.RData
#       Full R workspace, including intermediate objects.
#
# Reproducibility:
#   Output file names carry no date, so re-running this script overwrites
#   the files used in the reported analyses (originally generated on
#   2026-09-15). Filtering and normalization are deterministic, but
#   imputation is random (no seed set; see Step 5), so re-running gives
#   slightly different imputed values and downstream results. Keep a copy
#   of the reported files before re-running.
#
# Run in RStudio: setwd() and View() below require an interactive session.
# =============================================================================

# Load required packages.
# rio: import() / export(); readxl: read_excel();
# NormalyzerDE: medianNormalization(); missForest: missForest();
# imputeLCMD: impute.QRILC(); ggplot2 + ggbiplot: PCA plots.
library(rio)
library(readxl)
library(NormalyzerDE)
library(missForest)
library(imputeLCMD)
library(ggplot2)
library(ggbiplot)


# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Step 0: Study settings and run log -------------------------------------------

study_number <- "UW.Madison.rat.lung.proteomics"   # prefix of output file names
data_file <- "report.pg_matrix.tsv"                # not used below
search_method <- "DIANN"                           # used in the Excel file name
sample_info <- "Sample_metadata.xlsx"
species <- "rat"                                   # not used below

# Order of the 16 experimental groups (factor levels and output columns).
# Sample.group codes: 1st letter = age (Y = 7w, A = 7-11m),
# 2nd = lung region (T = Top, C = Central, P = Peripheral, B = Bottom),
# 3rd = sex (M/F); 6 samples per group.
sample_order <- c("YTM", "YCM", "YPM", "YBM", "YTF", "YCF", "YPF", "YBF", "ATM", "ACM", "APM", "ABM", "ATF", "ACF", "APF", "ABF")
# Y = young (7-week-old, 7w) and A = aged (7-11-month-old, 7-11m) in all following analyses.

# Start a new log file (overwrites any existing log).
log_file <- paste0(study_number, " analysis log.txt")
cat(paste("Study number:", study_number), file = log_file, sep = "\n")



# Step 1: Read intensities, annotations and sample metadata --------------------

# Read the annotated intensity table from Script_1.
# Converting to a data.frame and setting Uniprot.ID as row names is required:
# all later filtering and matching steps use these row names.
# (The reported data tables RData holds matrisome_anno_table as a tibble without
# row names; its contents are the same.)
matrisome_anno_table <- as.data.frame(read_excel("Protein_counts_summary.xlsx",
  sheet = "Raw_table_annotation"))
rownames(matrisome_anno_table) <- matrisome_anno_table$Uniprot.ID

# Split into intensities (data1: 3,859 proteins x 96 samples, raw values)
# and the 7 annotation columns.
data1 <- matrisome_anno_table[-c(1:7)]
matrisome_anno_table <- matrisome_anno_table[c(1:7)]




# Read the first sample-metadata worksheet and standardize its column names.
# Rows are in S1-S96 order, matching the columns of data1.
info_table <- import(sample_info)

colnames(info_table) <- c("Sample.ID", "Sample.Group", "Sample.PID", "Sample.Animal.Number", "Sample.Age", "Sample.Region", "Sample.Sex")

# Define the experimental-group order and encode age, region, and sex as factors.
info_table$Sample.Group <- factor(info_table$Sample.Group, levels = sample_order)
info_table$Sample.Age <- factor(info_table$Sample.Age)
info_table$Sample.Region <- factor(info_table$Sample.Region)
info_table$Sample.Sex <- factor(info_table$Sample.Sex)

cat(paste0(Sys.time(), " Total protein groups identified: ", nrow(data1), "\n"), file = log_file, append = TRUE)


# Step 2: Filter protein groups by detection rate ------------------------------

# Build a sample-by-protein detection matrix: non-missing = 1, missing = 0.
# work1: 96 rows (samples S1-S96) x 3,859 columns (Uniprot.ID).
# The loops assume 96 samples and 3,859 protein groups.
work1 <- as.data.frame(t(data1))

i = 1   # not needed; reset by the for loops
j = 1

for (i in 1:96) {
  for (j in 1:3859) {
    if (!is.na(work1[i, j])) {work1[i, j] <- 1}
    else {work1[i, j] <- 0}
  }
}


# Criterion 1 (global): detected in at least 48 of 96 samples (50%).
# work2: detection count per protein; work3: proteins passing (3,432).
work2 <- colSums(work1)
col_indices <- which(work2 >= 48)   # not used below
work3 <- work2[work2 >= 48]


# Count observed values within each of the 16 experimental groups.
# Grouping assumes metadata rows are aligned with the samples in work1.
# split() orders the groups by the factor levels, i.e. sample_order.
sample_groups <- split(work1, info_table$Sample.Group)
sample_group_sums <- lapply(sample_groups, colSums)

# sample_group_number: 3,859 proteins x 16 groups (values 0-6), in sample_order.
sample_group_number <- as.data.frame(cbind(sample_group_sums[[1]], sample_group_sums[[2]], sample_group_sums[[3]], sample_group_sums[[4]],
  sample_group_sums[[5]], sample_group_sums[[6]], sample_group_sums[[7]], sample_group_sums[[8]],
  sample_group_sums[[9]], sample_group_sums[[10]], sample_group_sums[[11]], sample_group_sums[[12]],
  sample_group_sums[[13]], sample_group_sums[[14]], sample_group_sums[[15]], sample_group_sums[[16]]

))
sample_group_number <- cbind(rownames(sample_group_number), sample_group_number)

# Column names follow sample_order.
colnames(sample_group_number) <- c("Uniprot", "YTM", "YCM", "YPM", "YBM",
  "YTF", "YCF", "YPF", "YBF",
  "ATM", "ACM", "APM", "ABM",
  "ATF", "ACF", "APF", "ABF"

)


# Export the group-level detection counts with protein identifiers.
export(sample_group_number, "Real_values_number_in_sample_group.xlsx")

# Drop the Uniprot column; row names still hold the Uniprot.ID.
sample_group_number <- sample_group_number[-1]


# Criterion 2 (conditional): detected in at least 5 of 6 samples in any
# group (3,596 proteins).
common_id <- c()
j = 1
i = 1

for (i in 1:3859) {
  if (any(sample_group_number[i, ] >= 5)) {
    common_id <- c(common_id, rownames(sample_group_number)[i])
  }
}


# Retain proteins meeting either the global or conditional detection criterion.
# data2: 3,600 proteins x 96 samples (raw intensities, original row order).
common_id1 <- union(common_id, names(work3))

data2 <- data1[rownames(data1) %in% common_id1, ] # 3600 proteins after filter


# Record retained protein counts and missing-value fractions
# (13.3% before and 8.5% after filtering in the reported run).
cat(paste0(Sys.time(), " Protein groups passed missing value filtering: ", nrow(data2), "\n"), file = log_file, append = TRUE)
cat(paste0(Sys.time(), " Rate of missing values before filtering: ", sum(is.na(data1)) / ncol(data1) / nrow(data1), "\n"), file = log_file, append = TRUE)
cat(paste0(Sys.time(), " Rate of missing values after filtering: ", sum(is.na(data2)) / ncol(data2) / nrow(data2), "\n"), file = log_file, append = TRUE)

View(head(data2))


# Step 3: Log2 transformation and median normalization -------------------------

# Transform intensities to log2 values, then apply median normalization.
# The data2 object is overwritten with log2-transformed filtered intensities.
# abs() has no effect here (all intensities are positive).
data2 <- log2(abs(data2))

# data3: median-normalized log2 intensities (input values already log2).
data3 <- medianNormalization(as.matrix(data2), noLogTransform = T)
data3 <- as.data.frame(data3)
rownames(data3) <- rownames(data2)

View(head(data3))


# Step 4: QC plots -------------------------------------------------------------

# Summarize missingness per protein and detection counts per sample for QC.
# Suffix 1 = before filtering (data1), 2 = after filtering (data2).
na.count1 <- rowSums(is.na(data1))
pro.count1 <- colSums(!is.na(data1))
na.count2 <- rowSums(is.na(data2))
pro.count2 <- colSums(!is.na(data2))


# Save QC plots showing missingness, detection counts, and intensity distributions
# before and after filtering or normalization.
# Red dotted line in the histograms = 50% of samples (48 NA).
pdf(paste0(study_number, "_QC plots.pdf"))

hist(na.count1, freq = F, xlab = "NA count", ylab = "% protein id", main = "Missing value/protein before filtering", col = 4, xlim = c(0, ncol(data1)))
abline(v = ncol(data1) * 0.5, col = "red", lty = 3, lwd = 2)

barplot(pro.count1, las = 2, xlab = "", ylab = "# of protein identified", main = "protein IDs in each sample before filtering", col = 4)

hist(na.count2, freq = F, xlab = "NA count", ylab = "% protein id", main = "Missing value/protein after filtering", col = 4, xlim = c(0, ncol(data1)))
abline(v = ncol(data1) * 0.5, col = "red", lty = 3, lwd = 2)

barplot(pro.count2, las = 2, xlab = "", ylab = "# of protein identified", main = "protein IDs in each sample after filtering", col = 4)

boxplot(data2, col = 'red', boxlwd = 1.5, medlwd = 1.5, whisklwd = 1.5, staplelwd = 1.5, outlwd = 1.5, cex.axis = 0.6, font.axis = 2, xlab = "Sample", ylab = "Log2 protein intensity", main = "Boxplot before normalization")

boxplot(data3, col = 'red', boxlwd = 1.5, medlwd = 1.5, whisklwd = 1.5, staplelwd = 1.5, outlwd = 1.5, cex.axis = 0.6, font.axis = 2, xlab = "Sample", ylab = "Log2 protein intensity", main = "Boxplot after normalization")

dev.off()


# Step 5: Two-stage imputation of missing values -------------------------------

# Imputation stage 1: apply random forest imputation to the normalized matrix.
# All missing values are imputed (29,324 cells); observed values are unchanged.
# Extract the completed matrix and record the out-of-bag error.
data4 <- missForest(data3)
data5 <- data4$ximp

cat(paste0(Sys.time(), " Imputation has been performed using Random Forest, OOB error of RF = ", data4$OOBerror, "\n"), file = log_file, append = TRUE)

View(head(data5))


# Build a mask for proteins with no observed values in an experimental group.
# work4: 3,859 proteins x 96 samples; 1 = keep, 0 = protein not detected in
# any of the 6 samples of that sample's group.
# The mask assumes 3,859 protein rows and sample columns S1 through S96.
# (The loop visits each group 6 times; the result is the same.)
work4 <- matrix(1, nrow = 3859, ncol = 96)
colnames(work4) <- paste0("S", 1:96)

for (group in info_table$Sample.Group) {
  zero_proteins <- which(sample_group_number[[group]] == 0)
  group_samples <- info_table$Sample.ID[info_table$Sample.Group == group]
  work4[zero_proteins, group_samples] <- 0
}


# Subset the mask to the 3,600 retained proteins and check row alignment.
work4 <- as.data.frame(work4)

rownames(work4) <- rownames(data1)
work5 <- work4[rownames(work4) %in% rownames(data5), ]
all(rownames(work5) == rownames(data5))   # TRUE


# Apply the mask to the random forest estimates and convert zero entries to NA.
# This restores missingness for group-wise absence before the second imputation,
# so that proteins absent from a whole group get low (left-censored) values
# instead of random forest estimates.
data6 <- data5 * work5
data6[data6 == 0] <- NA


# Imputation stage 2: impute remaining missing values using QRILC
# (quantile regression imputation of left-censored data), with tune.sigma = 0.5.
# data7: final imputed matrix (3,600 proteins x 96 samples, no NA).
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
# The logged number (2,652 in the reported run) is the number of intensity
# values (matrix cells) re-imputed by QRILC.
cat(paste0(Sys.time(), " Imputation has been performed using QRILC on the basis of RF, protein intensity numbers changed:", sum(col_zeros), "\n"), file = log_file, append = TRUE)



# Step 6: PCA plots ------------------------------------------------------------

# Run PCA on the imputed matrix with samples as observations (prcomp(t(data7))),
# using ggbiplot() from the ggbiplot package. Ellipses and colors show the
# grouping factor; point shapes show a second factor.


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


# Export the age, sex and region PCA views to a combined PDF.
pdf("PCA.pdf")
plot(g2)
plot(g3)
plot(g5)
dev.off()


# Step 7: Save data tables -----------------------------------------------------

# Collect original, filtered/log2-transformed, normalized, and imputed matrices.
# data_ori:      raw intensities, 3,859 proteins
# data_filtered: log2 intensities, 3,600 proteins (with NA)
# data_nor:      median-normalized log2 intensities, 3,600 proteins (with NA)
# data_imp:      imputed normalized log2 intensities, 3,600 proteins (no NA)
data_table <- list("data_ori" = data1,
  "data_filtered" = data2,
  "data_nor" = data3,
  "data_imp" = data7
)


# Save the analysis matrices, annotation object, and sample metadata
# (read by Script_4, Script_6 and Script_8).
save(data_table, matrisome_anno_table, info_table, file = paste0(study_number, "_data_tables.RData"))


# Inspect annotation ordering and protein-row consistency across matrices.
# These expressions report logical results without stopping execution
# (all TRUE in the reported run).
matrisome_anno_table1 <- matrisome_anno_table[order(matrisome_anno_table$Uniprot.ID), ]
all(matrisome_anno_table1[[1]] == matrisome_anno_table[[1]])   # rows sorted by Uniprot.ID
remove(matrisome_anno_table1)
all(rownames(data1) == matrisome_anno_table$Uniprot.ID)
all(rownames(data2) == rownames(data3))
all(rownames(data2) == rownames(data5))
all(rownames(data2) == rownames(data7))


# Match annotations to protein identifiers separately for each matrix
# and prepare the four annotated worksheets for export
# (7 annotation columns followed by S1-S96).
data_table_export <- list("data_ori" = cbind(matrisome_anno_table[match(rownames(data1), matrisome_anno_table$Uniprot.ID), ], data1),
  "data_filtered" = cbind(matrisome_anno_table[match(rownames(data2), matrisome_anno_table$Uniprot.ID), ], data2),
  "data_nor" = cbind(matrisome_anno_table[match(rownames(data3), matrisome_anno_table$Uniprot.ID), ], data3),
  "data_imp" = cbind(matrisome_anno_table[match(rownames(data7), matrisome_anno_table$Uniprot.ID), ], data7)
)


# Export annotated intensity tables (read by Script_4, Script_6 and Script_8)
# and record completion in the log.
export(data_table_export, file = paste0(study_number, "_", search_method, "_protein_intensity.xlsx"))
cat(paste0(Sys.time(), " Excel file of protein intensity data tables is created", "\n"), file = log_file, append = TRUE)


# Save the current R workspace.
save.image("Script_3_data_preprocessing_data_summary.RData")
