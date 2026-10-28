# =============================================================================
# Script 03: Supplementary Table S10 - batch-effect assessment
# =============================================================================
#
# Purpose:
#   Build Supplementary Table S10 (batch-effect assessment of sample
#   preparation day and LC-MS/MS injection order) from the results of
#   Script_01, in the same layout as the other supplementary tables
#   (Table S1): per worksheet, a bold title in row 1, italic notes in rows
#   2-3, the column header in row 5 (bold white text on blue, centred) and
#   the data from row 6; Arial 12, centred cells, no gridlines, frozen header.
#   Labels follow the manuscript: 7-11-month-old (7-11m) / 7-week-old (7w),
#   group codes such as 7-11m-F-B, "estimate" (model-derived log2
#   difference) and differentially expressed proteins (DEPs).
#
# Worksheets:
#   Sample preparation       96 samples: annotation, preparation day (1-4)
#                            and LC-MS/MS injection position (1-96)
#   Design balance           sample counts per preparation day for age, sex
#                            and region; chi-square / Fisher P values
#                            (preparation day) and Kruskal-Wallis P values
#                            (injection order)
#   Sensitivity summary      per comparison and pooled: DEPs in the original
#                            and adjusted model, retained / lost / gained,
#                            direction concordance, Pearson r, median change
#   Sensitivity per protein  per protein x comparison: estimate, FDR and DEP
#                            call in both models, with gene symbol and
#                            matrisome category
#   Matrisome DEP status     the 38 matrisome proteins of Fig. 4h x the 8
#                            comparisons: "estimate (FDR)" in both models
#                            ("filtered" if removed by the detection filter);
#                            DEP cells in bold, cells whose significance group
#                            changed highlighted in orange
#
# Inputs:
#   - results/batch_results.RData  (output of Script_01): cmp, workdata.
#   - data/sample_metadata.xlsx: sample annotation and technical variables.
#   - data/UW.Madison.rat.lung.proteomics_data_tables.RData  (Script_3):
#       matrisome_anno_table for gene symbols and matrisome categories.
#   - data/38 matrisome DEPs with unique features in 8 pairwise comparisons.txt
#       Uniprot.ID / Gene.symbol of the Fig. 4h proteins, in heatmap row order.
#
# Outputs:
#   - tables/Supplementary Table 10 Batch effect assessment of sample
#     preparation day and LC-MS-MS injection order.xlsx
#   - results/03_supplementary_table_workspace.RData  (full R workspace)
#
# How to run:
#   After Script_01, from the batch_effect_analysis folder:
#       Rscript code/03_supplementary_table.R
#   or Source it in RStudio. Optional arguments: [indir] [resultsdir] [outdir].
#   Runs in a few seconds.
#
# Reproducibility:
#   The Fisher exact test uses Monte Carlo P values (1e5 replicates); the
#   seed is set before each test (same values as in Script_01).
# =============================================================================

# Load required packages.
# dplyr: table handling; readxl: read_excel(); openxlsx: formatted workbook.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
})


# Step 0: Working directory and settings ---------------------------------------

# Locate the batch_effect_analysis folder (parent of code/), whether the
# script is run with Rscript, sourced in RStudio, or run line by line in
# RStudio. All paths below are relative to that folder.
find_root <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) return(normalizePath(file.path(dirname(f), "..")))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(normalizePath(file.path(dirname(p), "..")))
  }
  getwd()
}
setwd(find_root())
if (!dir.exists("data")) stop("Cannot find the data/ folder. Set the working directory to batch_effect_analysis/, ",
                              "e.g. setwd('.../batch effect visit/batch_effect_analysis'). Current: ", getwd())

# Command-line arguments (none in RStudio, so the defaults are used).
args <- commandArgs(TRUE)
indir  <- if (length(args) >= 1) args[1] else "data"      # input data folder
resdir <- if (length(args) >= 2) args[2] else "results"   # Script_01 output folder
outdir <- if (length(args) >= 3) args[3] else "tables"    # table output folder
out_file <- "Supplementary Table 10 Batch effect assessment of sample preparation day and LC-MS-MS injection order.xlsx"

# Script_01 results: cmp, tech, vpart, pc_assoc, robust, overall, workdata.
load(file.path(resdir, "batch_results.RData"))

# Manuscript labels. Age levels in workdata are Young / Aged.
age_lab <- c(Young = "7-week-old (7w)", Aged = "7-11-month-old (7-11m)")
grp <- function(age, region, sex) paste(age, substr(sex, 1, 1), substr(region, 1, 1), sep = "-")

# Comparison label as in the manuscript, e.g. "7-11m-F-B versus 7w-F-B".
cmp <- cmp %>%
  mutate(Comparison = paste(grp("7-11m", Region, Sex), "versus", grp("7w", Region, Sex)))

# Fig. 4h column order: BF, BM, CF, CM, PF, PM, TF, TM.
reg <- rep(c("Bottom", "Central", "Peripheral", "Top"), each = 2)
sx  <- rep(c("Female", "Male"), 4)
comp_order <- paste(grp("7-11m", reg, sx), "versus", grp("7w", reg, sx))
comp_short <- paste(reg, sx)


# Step 1: Sheet "Sample preparation" -------------------------------------------

# One row per sample in S1-S96 order, with preparation day and injection
# position.
meta <- read_excel(file.path(indir, "sample_metadata.xlsx"), sheet = 1)
S_A <- meta %>%
  select(`Sample ID`, `Sample Group`, PID, `Investigator Animal Number`, Age, Sex, Region,
         `Preparation Day`, `LC-MS/MS Injection Order` = `Injection Order`) %>%
  arrange(as.integer(sub("S", "", `Sample ID`)))


# Step 2: Sheet "Design balance" -----------------------------------------------

# Sample counts per level of age, sex and region on each preparation day.
S_B1 <- bind_rows(lapply(c("Age", "Sex", "Region"), function(v) {
  t <- as.data.frame.matrix(table(workdata[[v]], workdata$PreparationDay))
  colnames(t) <- paste("Day", colnames(t))
  lev <- rownames(t)
  if (v == "Age") lev <- age_lab[lev]
  cbind(Factor = v, Level = unname(lev), t, Total = rowSums(t))
}))

# Association of each biological factor with preparation day (chi-square,
# Fisher exact with Monte Carlo P value) and difference in injection order
# between factor levels (Kruskal-Wallis). Sample-level tests (n = 96).
S_B2 <- bind_rows(lapply(c("Age", "Sex", "Region"), function(v) {
  tab <- table(workdata[[v]], workdata$PreparationDay)
  data.frame(
    Factor = v,
    `Chi-square P (preparation day)`     = round(suppressWarnings(chisq.test(tab)$p.value), 3),
    `Fisher exact P (preparation day)`   = round({ set.seed(1); fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value }, 3),
    `Kruskal-Wallis P (injection order)` = round(kruskal.test(workdata$InjectionOrder ~ workdata[[v]])$p.value, 3),
    check.names = FALSE)
}))


# Step 3: Sheet "Sensitivity summary" ------------------------------------------

# Per comparison and pooled. Retained = DEP in both models; Lost = original
# only; Gained = adjusted only. Same direction = % of original DEPs with the
# same sign of the estimate in the adjusted model. Pearson r over all
# proteins tested in the comparison.
summ <- function(d) {
  d %>% summarise(
    `Proteins tested` = n(),
    `DEPs, original model` = sum(DEP0),
    `DEPs, adjusted model` = sum(DEP1),
    Retained = sum(DEP0 & DEP1),
    Lost     = sum(DEP0 & !DEP1),
    Gained   = sum(!DEP0 & DEP1),
    `Retained (%)` = round(100 * sum(DEP0 & DEP1) / sum(DEP0), 1),
    `Same direction among original DEPs (%)` = round(100 * mean(sign(estimate0[DEP0]) == sign(estimate1[DEP0])), 1),
    `Pearson r of estimates` = round(cor(estimate0, estimate1), 3),
    `Median absolute change in estimate` = round(median(abs(estimate1 - estimate0)), 3))
}
S_C1 <- bind_rows(
  cmp %>% group_by(Comparison) %>% summ() %>% ungroup() %>% arrange(match(Comparison, comp_order)),
  cmp %>% summ() %>% mutate(Comparison = "All 8 comparisons", .before = 1))


# Step 4: Sheet "Sensitivity per protein" --------------------------------------

# Gene symbol and matrisome category from the Script_3 annotation table.
e <- new.env()
load(file.path(indir, "UW.Madison.rat.lung.proteomics_data_tables.RData"), envir = e)
genes <- e$matrisome_anno_table %>% select(Uniprot.ID, Gene.symbol, Matrisome.Category)

# One row per protein x comparison (proteins passing the comparison-specific
# detection filter), comparisons in Fig. 4h order.
S_C2 <- cmp %>%
  transmute(Uniprot.ID, Comparison,
            `Estimate, original model` = estimate0,
            `FDR, original model` = p.value.adjusted0,
            `DEP, original model` = DEP0,
            `Estimate, adjusted model` = estimate1,
            `FDR, adjusted model` = p.value.adjusted1,
            `DEP, adjusted model` = DEP1)
S_C2 <- left_join(genes, S_C2, by = "Uniprot.ID", multiple = "all") %>%
  filter(!is.na(Comparison)) %>%
  arrange(Uniprot.ID, match(Comparison, comp_order))


# Step 5: Sheet "Matrisome DEP status" (Fig. 4h proteins) ----------------------

# Fig. 4h protein list, in heatmap row order.
mat38 <- read.delim(file.path(indir, "38 matrisome DEPs with unique features in 8 pairwise comparisons.txt"))
stopifnot(nrow(mat38) == 38)

# Cell content: "estimate (FDR)" for every tested protein, e.g.
# "-1.19 (3.5e-06)"; "filtered" if removed by the detection filter (blank in
# Fig. 4h). The significance group (higher in 7-11m, higher in 7w, not DEP)
# is shown by formatting in Step 6: DEP cells in bold, cells whose group
# differs between the two models highlighted.
cell <- function(est, fdr, dep) {
  ifelse(is.na(dep), "filtered", sprintf("%.2f (%s)", est, formatC(fdr, digits = 2, format = "g")))
}
grid <- expand.grid(Uniprot.ID = mat38$Uniprot.ID, Comparison = comp_order, stringsAsFactors = FALSE) %>%
  left_join(cmp %>% select(Uniprot.ID, Comparison, estimate0, p.value.adjusted0, DEP0,
                           estimate1, p.value.adjusted1, DEP1),
            by = c("Uniprot.ID", "Comparison")) %>%
  mutate(orig = cell(estimate0, p.value.adjusted0, DEP0),
         adj  = cell(estimate1, p.value.adjusted1, DEP1),
         changed = !is.na(DEP0) & ((DEP0 != DEP1) | (DEP0 & sign(estimate0) != sign(estimate1))))

# Wide table: one row per protein (gene symbols in upper case, as in
# Fig. 4h); for each comparison an "original" and an "adjusted" column.
S_C3 <- data.frame(Uniprot.ID = mat38$Uniprot.ID, Gene.symbol = toupper(mat38$Gene.symbol), check.names = FALSE)
for (k in seq_along(comp_order)) {
  g <- grid[grid$Comparison == comp_order[k], ]
  g <- g[match(mat38$Uniprot.ID, g$Uniprot.ID), ]
  S_C3[[paste0(comp_short[k], ", original")]] <- g$orig
  S_C3[[paste0(comp_short[k], ", adjusted")]] <- g$adj
}
S_C3$`Comparisons with changed status` <- sapply(mat38$Uniprot.ID, function(id) sum(grid$changed[grid$Uniprot.ID == id]))

n_tested  <- sum(!is.na(grid$DEP0))
n_changed <- sum(grid$changed)


# Step 6: Write the workbook in the supplementary-table style ------------------

# Styles copied from Table S1: Arial 12, dark grey text, centred cells;
# header bold white on blue (#4472C4); bold title; italic notes.
base_font <- list(fontName = "Arial", fontSize = 12)
st_body   <- do.call(createStyle, c(base_font, list(fontColour = "#1F2937", halign = "center")))
st_header <- do.call(createStyle, c(base_font, list(fontColour = "#FFFFFF", fgFill = "#4472C4",
                     textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE)))
st_title  <- do.call(createStyle, c(base_font, list(textDecoration = "bold")))
st_note   <- do.call(createStyle, c(base_font, list(textDecoration = "italic")))
st_dep    <- do.call(createStyle, c(base_font, list(fontColour = "#1F2937", halign = "center", textDecoration = "bold")))
st_chg    <- do.call(createStyle, c(base_font, list(fontColour = "#1F2937", halign = "center", fgFill = "#FCE4B6")))
st_both   <- do.call(createStyle, c(base_font, list(fontColour = "#1F2937", halign = "center", fgFill = "#FCE4B6",
                     textDecoration = "bold")))

wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial")

# Write one table with its header row in the Table S1 style.
put_table <- function(sheet, df, row) {
  writeData(wb, sheet, df, startRow = row, headerStyle = st_header)
  addStyle(wb, sheet, st_body, rows = row + seq_len(nrow(df)), cols = seq_len(ncol(df)), gridExpand = TRUE)
  setRowHeights(wb, sheet, rows = row, heights = 21)
}

# One worksheet: title (row 1), up to two note lines (rows 2-3), blank row 4,
# table from row 5, header frozen, gridlines off.
add_sheet <- function(sheet, title, notes, df, widths) {
  addWorksheet(wb, sheet, gridLines = FALSE)
  writeData(wb, sheet, title, startRow = 1)
  addStyle(wb, sheet, st_title, rows = 1, cols = 1)
  writeData(wb, sheet, notes, startRow = 2)
  addStyle(wb, sheet, st_note, rows = 1 + seq_along(notes), cols = 1)
  put_table(sheet, df, 5)
  freezePane(wb, sheet, firstActiveRow = 6)
  setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = widths)
}

add_sheet("Sample preparation", "Sample preparation and acquisition order",
  c("96 rat lung samples; one row per sample, with group code, animal ID, age, sex, region, sample preparation day (1-4) and LC-MS/MS injection position (1-96).",
    "24 samples were prepared per day; Sample ID connects to Table S1."),
  S_A, c(16, 20, 10, 31, 25, 15, 19, 20, 30))

# Design balance: count table from row 5, test table below it.
add_sheet("Design balance", "Distribution of biological groups across sample preparation days",
  c("Number of samples per preparation day. P values test the association of age, sex and region with preparation day (chi-square; Fisher exact test, 1e5 Monte Carlo replicates)",
    "and differences in LC-MS/MS injection order between groups (Kruskal-Wallis). Tests are at sample level (n = 96)."),
  S_B1, c(14, 30, 34, 34, 36, 10, 10))
put_table("Design balance", S_B2, nrow(S_B1) + 8)

model_notes <- c(
  "Original model: Protein ~ Age + Sex + Region + Age:Sex + Age:Region + Sex:Region + (1 | Subject); adjusted model: original model + sample preparation day (categorical) + LC-MS/MS injection order (continuous, scaled).",
  "7-11-month-old minus 7-week-old estimates within each sex and region; comparison-specific detection filter and Benjamini-Hochberg correction as in the Methods; DEP: FDR < 0.05 and |estimate| > log2(1.5).")

add_sheet("Sensitivity summary",
  "7-11-month-old versus 7-week-old multiple-factor comparisons before and after adjustment for preparation day and injection order",
  model_notes, S_C1, c(28, 18, 22, 22, 12, 10, 10, 15, 26, 18, 24))

add_sheet("Sensitivity per protein",
  "Per-protein 7-11-month-old versus 7-week-old estimates before and after adjustment for preparation day and injection order",
  model_notes, S_C2, c(16, 16, 22, 28, 18, 18, 16, 18, 18, 16))

add_sheet("Matrisome DEP status",
  "Differential expression status of the 38 matrisome proteins in Fig. 4h before and after adjustment for preparation day and injection order",
  c("Cell values: estimate (FDR). Bold, differentially expressed (FDR < 0.05 and |estimate| > log2(1.5)); filtered, removed by the comparison-specific detection filter (blank in Fig. 4h).",
    sprintf("Orange, significance group (higher in 7-11m, higher in 7w, or not differentially expressed) differs between models: %d of %d tested protein-comparison pairs; unchanged in %d (%.1f%%).",
            n_changed, n_tested, n_tested - n_changed, 100 * (n_tested - n_changed) / n_tested)),
  S_C3, c(14, 14, rep(17, 16), 18))

# Matrisome sheet: bold = DEP in that model; orange = group changed between
# models (table body starts in row 6; comparison k occupies columns 2k+1
# (original) and 2k+2 (adjusted)).
for (k in seq_along(comp_order)) {
  g <- grid[grid$Comparison == comp_order[k], ]
  g <- g[match(mat38$Uniprot.ID, g$Uniprot.ID), ]
  for (j in 1:2) {
    dep <- if (j == 1) g$DEP0 else g$DEP1
    col <- 2 * k + j
    r_dep  <- which(dep %in% TRUE & !g$changed)
    r_chg  <- which(g$changed & !(dep %in% TRUE))
    r_both <- which(g$changed & dep %in% TRUE)
    if (length(r_dep))  addStyle(wb, "Matrisome DEP status", st_dep,  rows = r_dep + 5,  cols = col)
    if (length(r_chg))  addStyle(wb, "Matrisome DEP status", st_chg,  rows = r_chg + 5,  cols = col)
    if (length(r_both)) addStyle(wb, "Matrisome DEP status", st_both, rows = r_both + 5, cols = col)
  }
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, file.path(outdir, out_file), overwrite = TRUE)

# Print the key results.
print(as.data.frame(S_C1), width = 200)
print(S_B2)
cat(sprintf("Matrisome Fig. 4h: %d tested pairs, %d changed\n", n_tested, n_changed))


# Save the full R workspace, including intermediate objects.
save.image(file.path(resdir, "03_supplementary_table_workspace.RData"))
