# =============================================================================
# Script 10: Prepare the helper data for the two Shiny applications
# =============================================================================
#
# Purpose:
#   Package the multiple-factor differential-expression results of Script_6
#   (8 age, 8 sex and 24 region comparisons) for the two Shiny apps in
#   apps/. The apps do not refit models or recompute adjusted p-values.
#   (Young = 7-week-old, 7w; Aged = 7-11-month-old, 7-11m.)
#
# Input (in this folder):
#   - Script_6_multiple_factor_DE.RData
#       Copy of the Script_6 workspace (same results as
#       ../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.RData).
#       *_comparisons_unfiltered: model results for all 3,600 proteins
#       (no p-value adjustment); age.vs.young4, female.vs.male4,
#       region.vs.region4: Method 4 results (proteins with >= 5 observed
#       values in either group, BH-adjusted within each comparison);
#       df0: the 7 protein annotation columns.
#
# Outputs (written to this folder and copied into apps/, see Step 5):
#   - helper age.RData, helper sex.RData, helper region.RData
#       (for apps/protein_explorer)
#       df1: annotations + estimate and FDR for every comparison (one row per
#            protein with a gene symbol; FDR is NA when the protein was
#            removed by the Method 4 filter in that comparison);
#       df4: transposed display matrix (rows = estimates, then FDRs, in the
#            app's comparison order; columns = gene symbols), 4 significant
#            digits;
#       df6: protein annotations sorted by gene symbol.
#   - helper volcano plot.RData  (for apps/volcano_explorer)
#       age.vs.young, female.vs.male, region.vs.region: Method 4 result
#       tables, named by their display label and ordered as shown in the app.
#   - sessionInfo.txt
#       R and package versions of the session that ran this script.
#
# Comparison order used by both apps:
#   Age (Aged versus Young):   Bottom / Top / Central / Peripheral x Male / Female
#   Sex (Female versus Male):  Bottom / Top / Central / Peripheral x Aged / Young
#   Region (6 region pairs):   within Young Male, Young Female, Aged Male,
#                              Aged Female
#   Positive estimate = higher in the first group (Aged, Female, or the
#   first region of the pair).
#
# Run in RStudio: setwd() below requires an interactive RStudio session.
# =============================================================================

# Set the working directory to this script's folder (RStudio only).
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages.
# dplyr: arrange().
library(dplyr)


# Step 1: Load the Script_6 results --------------------------------------------

load("Script_6_multiple_factor_DE.RData")


# Step 2: Comparison order and display labels ----------------------------------

# App display label = Script_6 list name, in the order shown in the apps.
age_order <- c(
  "Bottom Male"       = "bottom.male",
  "Bottom Female"     = "bottom.female",
  "Top Male"          = "top.male",
  "Top Female"        = "top.female",
  "Central Male"      = "central.male",
  "Central Female"    = "central.female",
  "Peripheral Male"   = "peri.male",
  "Peripheral Female" = "peri.female"
)

sex_order <- c(
  "Bottom Aged"      = "bottom.age",
  "Bottom Young"     = "bottom.young",
  "Top Aged"         = "top.age",
  "Top Young"        = "top.young",
  "Central Aged"     = "central.age",
  "Central Young"    = "central.young",
  "Peripheral Aged"  = "peri.age",
  "Peripheral Young" = "peri.young"
)

region_pairs <- c(
  "Bottom versus Central"     = "bottom.vs.central",
  "Bottom versus Peripheral"  = "bottom.vs.peri",
  "Bottom versus Top"         = "bottom.vs.top",
  "Central versus Peripheral" = "central.vs.peri",
  "Central versus Top"        = "central.vs.top",
  "Peripheral versus Top"     = "peri.vs.top"
)
region_groups <- c(
  "Young Male"   = "male.young",
  "Young Female" = "female.young",
  "Aged Male"    = "male.age",
  "Aged Female"  = "female.age"
)
region_order <- unlist(lapply(names(region_groups), function(g) {
  setNames(paste(region_pairs, region_groups[[g]], sep = "."), paste0(g, ": ", names(region_pairs)))
}))

stopifnot(
  setequal(age_order, names(age.vs.young4)),
  setequal(sex_order, names(female.vs.male4)),
  setequal(region_order, names(region.vs.region4))
)


# Step 3: Protein explorer helpers ---------------------------------------------

# For each comparison, take the unfiltered model results and fill in the
# Method 4 FDR (NA = protein removed by the filter in that comparison).
fill_fdr <- function(unfiltered, filtered) {
  unfiltered$p.value.adjusted <- filtered$p.value.adjusted[match(unfiltered$Uniprot.ID, filtered$Uniprot.ID)]
  unfiltered
}

make_protein_helper <- function(unfiltered_list, filtered_list, order, file) {
  tables <- lapply(order, function(k) fill_fdr(unfiltered_list[[k]], filtered_list[[k]]))

  # df1: annotations + (estimate, FDR) for each comparison, in app order.
  df1 <- cbind(df0, do.call(cbind, unname(lapply(tables, function(t) t[c("estimate", "p.value.adjusted")]))))
  df1 <- df1[df1$Gene.symbol != "-", ] # proteins with a gene symbol only
  rownames(df1) <- df1$Gene.symbol

  # df4: rows = all estimates, then all FDRs; columns = gene symbols.
  stats <- df1[-(1:ncol(df0))]
  n <- length(order)
  df4 <- as.data.frame(t(stats[c(seq(1, 2 * n, 2), seq(2, 2 * n, 2))]))
  df4 <- signif(df4, 4)

  # df6: protein annotations sorted by gene symbol (English collation).
  df6 <- df0 %>% arrange(Gene.symbol, .locale = "en")

  save(df1, df4, df6, file = file)
}

make_protein_helper(age_comparisons_unfiltered, age.vs.young4, age_order, "helper age.RData")
make_protein_helper(sex_comparisons_unfiltered, female.vs.male4, sex_order, "helper sex.RData")
make_protein_helper(region_comparisons_unfiltered, region.vs.region4, region_order, "helper region.RData")


# Step 4: Volcano explorer helper ----------------------------------------------

# Method 4 result tables, named by display label and in app order, so that
# the app finds each comparison by name (not by list position).
age.vs.young <- setNames(age.vs.young4[age_order],
                         paste0(names(age_order), ": Aged versus Young"))
female.vs.male <- setNames(female.vs.male4[sex_order],
                           paste0(names(sex_order), ": Female versus Male"))
region.vs.region <- setNames(region.vs.region4[region_order], names(region_order))

save(age.vs.young, female.vs.male, region.vs.region, file = "helper volcano plot.RData")


# Step 5: Copy the helper files into the apps ---------------------------------

# Each app loads its own copy, so it can run or be deployed on its own.
file.copy(c("helper age.RData", "helper sex.RData", "helper region.RData"), "apps/protein_explorer", overwrite = TRUE)
file.copy("helper volcano plot.RData", "apps/volcano_explorer", overwrite = TRUE)


# Record the R and package versions of this session.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
