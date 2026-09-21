# Load required packages.
library(ggplot2)
library(ggrepel)
library(readxl)
library(UpSetR)
library(dplyr)

# Optional: set the working directory to this script's folder in RStudio.
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Female-enriched protein sets: Figure S10a --------------------------------------
female_group <- read_excel(
  "Female vs male DEP in 8 pairwise groups.xlsx",
  sheet = "8 groups up in female"
)

# Extract protein lists for the eight age-by-region groups.
Bottom_Aged <- female_group$bottom_aged
Central_Aged <- female_group$central_aged
Peri_Aged <- female_group$peri_aged
Top_Aged <- female_group$top_aged

Bottom_Young <- female_group$bottom_young
Central_Young <- female_group$central_young
Peri_Young <- female_group$peri_young
Top_Young <- female_group$top_young

# Remove missing entries from each protein list.
Bottom_Aged <- Bottom_Aged[!is.na(Bottom_Aged)]
Central_Aged <- Central_Aged[!is.na(Central_Aged)]
Peri_Aged <- Peri_Aged[!is.na(Peri_Aged)]
Top_Aged <- Top_Aged[!is.na(Top_Aged)]

Bottom_Young <- Bottom_Young[!is.na(Bottom_Young)]
Central_Young <- Central_Young[!is.na(Central_Young)]
Peri_Young <- Peri_Young[!is.na(Peri_Young)]
Top_Young <- Top_Young[!is.na(Top_Young)]

# Build a binary membership matrix for the eight age-by-region groups.
# Collect unique proteins across all eight groups.
all_proteins_female <- unique(c(
  Bottom_Aged, Central_Aged, Peri_Aged, Top_Aged,
  Bottom_Young, Central_Young, Peri_Young, Top_Young
))

# Initialize the matrix with proteins as rows and groups as columns.
binary_matrix_female <- matrix(0, nrow = length(all_proteins_female), ncol = 8)
rownames(binary_matrix_female) <- all_proteins_female
colnames(binary_matrix_female) <- c(
  "Bottom_Aged", "Central_Aged", "Peripheral_Aged", "Top_Aged",
  "Bottom_Young", "Central_Young", "Peripheral_Young", "Top_Young"
)

# Mark membership in each group with 1.
binary_matrix_female[Bottom_Aged, "Bottom_Aged"] <- 1
binary_matrix_female[Central_Aged, "Central_Aged"] <- 1
binary_matrix_female[Peri_Aged, "Peripheral_Aged"] <- 1
binary_matrix_female[Top_Aged, "Top_Aged"] <- 1
binary_matrix_female[Bottom_Young, "Bottom_Young"] <- 1
binary_matrix_female[Central_Young, "Central_Young"] <- 1
binary_matrix_female[Peri_Young, "Peripheral_Young"] <- 1
binary_matrix_female[Top_Young, "Top_Young"] <- 1

# Convert to a data frame for UpSetR.
upset_data_female <- as.data.frame(binary_matrix_female)

# Make Figure S10a.
upset_plot_female <- upset(
  upset_data_female,
  nsets = 8,
  order.by = "freq"
)

# Save Figure S10a without an extra blank page.
pdf("Figure_S10a_UpSet_up_in_female.pdf", width = 10, height = 7)
print(upset_plot_female, newpage = FALSE)
dev.off()


# Export female-enriched protein lists ------------------------------------------

# Use the membership data frame to identify exact intersections.
protein_df <- upset_data_female

# Identify observed membership patterns, excluding empty patterns.
intersection_patterns <- unique(protein_df)
intersection_patterns <- intersection_patterns[rowSums(intersection_patterns) > 0, ]

# Initialize containers for intersection protein lists and labels.
intersection_protein_lists <- list()
intersection_names <- c()

# Match proteins to each exact membership pattern.
for (i in seq_len(nrow(intersection_patterns))) {
  pattern <- intersection_patterns[i, ]

  # Collect proteins with the same membership across all eight groups.
  matching_proteins <- c()
  for (j in seq_len(nrow(protein_df))) {
    if (identical(as.numeric(protein_df[j, ]), as.numeric(pattern))) {
      matching_proteins <- c(matching_proteins, rownames(protein_df)[j])
    }
  }

  # Label the active groups for this exact intersection.
  active_groups <- names(pattern)[pattern == 1]
  if (length(active_groups) == 1) {
    intersection_name <- paste("Only", active_groups)
  } else {
    intersection_name <- paste(active_groups, collapse = " ∩ ")
  }

  intersection_protein_lists[[i]] <- matching_proteins
  intersection_names[i] <- intersection_name
}

# Sort intersections by decreasing protein count.
intersection_sizes <- sapply(intersection_protein_lists, length)
sorted_indices <- order(intersection_sizes, decreasing = TRUE)

# Export exact intersection labels, counts, and protein lists.
output_file <- "protein_lists_by_intersection_female.txt"
sink(output_file)
cat("PROTEIN LISTS FOR EACH INTERSECTION GROUP IN UPSET PLOT (FEMALE DATA)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

for (i in sorted_indices) {
  cat("INTERSECTION GROUP", which(sorted_indices == i), ":", intersection_names[i], "\n")
  cat("Number of proteins:", length(intersection_protein_lists[[i]]), "\n")
  cat("Proteins:\n")
  if (length(intersection_protein_lists[[i]]) > 0) {
    cat(intersection_protein_lists[[i]], sep = "\n")
  } else {
    cat("(No proteins in this intersection)\n")
  }
  cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")
}
sink()

# Export the complete protein list for each individual group.
output_file2 <- "protein_lists_by_individual_groups_female.txt"
sink(output_file2)
categories <- colnames(protein_df)
for (i in seq_along(categories)) {
  proteins_in_category <- rownames(protein_df)[protein_df[[categories[i]]] == 1]
  cat("Group:", categories[i], "\n")
  cat("Number of proteins:", length(proteins_in_category), "\n")
  cat("Proteins:\n")
  cat(proteins_in_category, sep = "\n")
  cat("\n\n")
}
sink()

cat("Intersection protein lists saved to", output_file, "\n")
cat("Individual group lists saved to", output_file2, "\n")


# Male-enriched protein sets: Figure S10b ----------------------------------------
male_group <- read_excel(
  "Female vs male DEP in 8 pairwise groups.xlsx",
  sheet = "8 groups up in male"
)

# Extract protein lists for the eight age-by-region groups.
Bottom_Aged <- male_group$bottom_aged
Central_Aged <- male_group$central_aged
Peri_Aged <- male_group$peri_aged
Top_Aged <- male_group$top_aged

Bottom_Young <- male_group$bottom_young
Central_Young <- male_group$central_young
Peri_Young <- male_group$peri_young
Top_Young <- male_group$top_young

# Remove missing entries from each protein list.
Bottom_Aged <- Bottom_Aged[!is.na(Bottom_Aged)]
Central_Aged <- Central_Aged[!is.na(Central_Aged)]
Peri_Aged <- Peri_Aged[!is.na(Peri_Aged)]
Top_Aged <- Top_Aged[!is.na(Top_Aged)]

Bottom_Young <- Bottom_Young[!is.na(Bottom_Young)]
Central_Young <- Central_Young[!is.na(Central_Young)]
Peri_Young <- Peri_Young[!is.na(Peri_Young)]
Top_Young <- Top_Young[!is.na(Top_Young)]

# Build a binary membership matrix for the eight age-by-region groups.
# Collect unique proteins across all eight groups.
all_proteins_male <- unique(c(
  Bottom_Aged, Central_Aged, Peri_Aged, Top_Aged,
  Bottom_Young, Central_Young, Peri_Young, Top_Young
))

# Initialize the matrix with proteins as rows and groups as columns.
binary_matrix_male <- matrix(0, nrow = length(all_proteins_male), ncol = 8)
rownames(binary_matrix_male) <- all_proteins_male
colnames(binary_matrix_male) <- c(
  "Bottom_Aged", "Central_Aged", "Peripheral_Aged", "Top_Aged",
  "Bottom_Young", "Central_Young", "Peripheral_Young", "Top_Young"
)

# Mark membership in each group with 1.
binary_matrix_male[Bottom_Aged, "Bottom_Aged"] <- 1
binary_matrix_male[Central_Aged, "Central_Aged"] <- 1
binary_matrix_male[Peri_Aged, "Peripheral_Aged"] <- 1
binary_matrix_male[Top_Aged, "Top_Aged"] <- 1
binary_matrix_male[Bottom_Young, "Bottom_Young"] <- 1
binary_matrix_male[Central_Young, "Central_Young"] <- 1
binary_matrix_male[Peri_Young, "Peripheral_Young"] <- 1
binary_matrix_male[Top_Young, "Top_Young"] <- 1

# Convert to a data frame for UpSetR.
upset_data_male <- as.data.frame(binary_matrix_male)

# Make Figure S10b.
upset_plot_male <- upset(
  upset_data_male,
  nsets = 8,
  order.by = "freq"
)

# Save Figure S10b without an extra blank page.
pdf("Figure_S10b_UpSet_up_in_male.pdf", width = 10, height = 7)
print(upset_plot_male, newpage = FALSE)
dev.off()


# Export male-enriched protein lists --------------------------------------------

# Use the membership data frame to identify exact intersections.
protein_df <- upset_data_male

# Identify observed membership patterns, excluding empty patterns.
intersection_patterns <- unique(protein_df)
intersection_patterns <- intersection_patterns[rowSums(intersection_patterns) > 0, ]

# Initialize containers for intersection protein lists and labels.
intersection_protein_lists <- list()
intersection_names <- c()

# Match proteins to each exact membership pattern.
for (i in seq_len(nrow(intersection_patterns))) {
  pattern <- intersection_patterns[i, ]

  # Collect proteins with the same membership across all eight groups.
  matching_proteins <- c()
  for (j in seq_len(nrow(protein_df))) {
    if (identical(as.numeric(protein_df[j, ]), as.numeric(pattern))) {
      matching_proteins <- c(matching_proteins, rownames(protein_df)[j])
    }
  }

  # Label the active groups for this exact intersection.
  active_groups <- names(pattern)[pattern == 1]
  if (length(active_groups) == 1) {
    intersection_name <- paste("Only", active_groups)
  } else {
    intersection_name <- paste(active_groups, collapse = " ∩ ")
  }

  intersection_protein_lists[[i]] <- matching_proteins
  intersection_names[i] <- intersection_name
}

# Sort intersections by decreasing protein count.
intersection_sizes <- sapply(intersection_protein_lists, length)
sorted_indices <- order(intersection_sizes, decreasing = TRUE)

# Export exact intersection labels, counts, and protein lists.
output_file <- "protein_lists_by_intersection male.txt"
sink(output_file)
cat("PROTEIN LISTS FOR EACH INTERSECTION GROUP IN UPSET PLOT\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

for (i in sorted_indices) {
  cat("INTERSECTION GROUP", which(sorted_indices == i), ":", intersection_names[i], "\n")
  cat("Number of proteins:", length(intersection_protein_lists[[i]]), "\n")
  cat("Proteins:\n")
  if (length(intersection_protein_lists[[i]]) > 0) {
    cat(intersection_protein_lists[[i]], sep = "\n")
  } else {
    cat("(No proteins in this intersection)\n")
  }
  cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")
}
sink()

# Export the complete protein list for each individual group.
output_file2 <- "protein_lists_by_individual_groups male.txt"
sink(output_file2)
categories <- colnames(protein_df)
for (i in seq_along(categories)) {
  proteins_in_category <- rownames(protein_df)[protein_df[[categories[i]]] == 1]
  cat("Group:", categories[i], "\n")
  cat("Number of proteins:", length(proteins_in_category), "\n")
  cat("Proteins:\n")
  cat(proteins_in_category, sep = "\n")
  cat("\n\n")
}
sink()

cat("Intersection protein lists saved to", output_file, "\n")
cat("Individual group lists saved to", output_file2, "\n")
