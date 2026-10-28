# =============================================================================
# Script 01: Batch-effect assessment and sensitivity analysis
#            (response to Reviewer comment 12)
# =============================================================================
#
# Purpose:
#   Assess whether the two technical factors of the experiment, sample
#   preparation day and LC-MS/MS injection order, affect the reported
#   Aged vs Young differential abundance results. Five analyses:
#     (1) Design balance: are age, sex and region associated with
#         preparation day or injection order? (chi-square / Fisher,
#         Kruskal-Wallis)
#     (2) Global structure: PCA of the imputed matrix; association (R2) of
#         PC1-PC10 with each biological and technical factor.
#     (3) Protein-level test: per protein, likelihood-ratio test (LRT)
#         comparing the original model (Model 0) with a model that also
#         contains the technical covariates (Model 1).
#     (4) Variance partitioning: per protein, fraction of variance explained
#         by each factor (variancePartition-style random-effects model).
#     (5) Sensitivity analysis: Aged vs Young contrasts from Model 0 and
#         Model 1, filtered and FDR-adjusted exactly as in Script_6
#         (Method 4), and compared (DEPs retained / lost / gained,
#         direction concordance, correlation of log2 fold changes).
#
# Models (one per protein, 3,600 proteins):
#   Model 0 (original, identical to Script_6):
#       Protein ~ Age + Sex + Region + Age:Region + Sex:Region + Age:Sex
#                 + (1 | PID)
#   Model 1 (batch-adjusted):
#       Model 0 + PreparationDay (factor, 4 levels)
#               + InjectionOrderScaled (injection position 1-96, z-scored)
#   Within-day preparation order is read and reported (Spearman
#   correlations) but is not a model covariate.
#
# Contrasts, filter and multiple testing (as in Script_6, Method 4):
#   Aged - Young within each Region x Sex (8 comparisons), emmeans with
#   Kenward-Roger df. Per comparison, keep proteins with >= 5 observed
#   (non-imputed) values in at least one of the two groups, then BH-adjust
#   within that comparison. DEP: FDR < 0.05 and |log2FC| > log2(1.5).
#   Model 0 reproduces the reported 1,448 DEP calls (670 unique proteins).
#
# Inputs (folder data/):
#   - UW.Madison.rat.lung.proteomics_data_tables.RData  (output of Script_3)
#       data_table$data_imp: imputed, normalized log2 intensities
#       (3,600 proteins x S1-S96); info_table: sample annotation.
#   - sample_metadata.xlsx
#       96 rows; columns "Sample ID", "Preparation Day" (1-4),
#       "Within-day Preparation Order" (1-24), "Injection Order" (1-96).
#   - Real_values_number_in_sample_group.xlsx  (output of Script_3)
#       Number of observed values (0-6) per protein in each of 16 groups.
#   - 38 matrisome DEPs with unique features in 8 pairwise comparisons.txt
#       Uniprot.ID / Gene.symbol of the 38 matrisome proteins (Figure 3h).
#
# Outputs (folder results/):
#   - batch_results.RData
#       Key objects read by Script_02 and Script_03:
#       cmp       per protein x comparison: estimate, SE, FDR and DEP call
#                 from Model 0 (suffix 0) and Model 1 (suffix 1)
#       tech      per protein: LRT and Type III F tests of technical terms
#       vpart     per protein: variance fractions
#       pc_assoc  PC1-PC10: variance explained and R2 with each factor
#       robust    per comparison: robustness summary
#       overall   all comparisons pooled
#       workdata  model table (96 samples)
#   - Batch_effect_sensitivity_analysis.xlsx
#       All intermediate tables, one worksheet each (see Step 8).
#   - Batch_effect_assessment_figures.pdf
#       Diagnostic plots (PCA by day / injection order, PC1 vs injection
#       order, variance partitioning, Model 0 vs Model 1 estimates, LRT
#       p-value histogram).
#   - 01_batch_effect_models_workspace.RData
#       Full R workspace, including intermediate objects.
#
# How to run:
#   From the batch_effect_analysis folder:
#       Rscript code/01_batch_effect_models.R
#   or open this file in RStudio and Source it. The working directory is set
#   to batch_effect_analysis/ automatically (Step 0).
#   Optional arguments: [indir] [outdir] [n_test]; n_test fits only the first
#   n proteins (quick test; FDR values are then not meaningful).
#   Runtime: ~20-30 min on one core (3,600 proteins x 2 models, plus an LRT
#   refit and a variance-partition model each). Model fitting is parallelized
#   with mclapply on macOS / Linux.
#
# Reproducibility:
#   Model fitting is deterministic. The Fisher exact test uses Monte Carlo
#   p-values; a seed is set before each test in Step 3 (same values as in
#   Script_03).
# =============================================================================

# Load required packages.
# lme4 + lmerTest: lmer(), Type III F tests; emmeans: contrasts (Kenward-Roger
# df requires pbkrtest); dplyr / tidyr: table handling; readxl: read_excel();
# openxlsx: write.xlsx(); ggplot2: diagnostic plots; parallel: mclapply().
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(parallel)
})


# Step 0: Working directory and settings ---------------------------------------

# Locate the batch_effect_analysis folder (parent of code/), whether the
# script is run with Rscript, sourced in RStudio, or run line by line in
# RStudio. All paths below (data/, results/) are relative to that folder.
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
indir  <- if (length(args) >= 1) args[1] else "data"      # input folder
outdir <- if (length(args) >= 2) args[2] else "results"   # output folder
n_test <- if (length(args) >= 3) as.integer(args[3]) else NA   # NA = all proteins

# Number of cores for model fitting (mclapply forks on macOS / Linux;
# set ncores <- 1 on Windows).
ncores <- max(1, detectCores() - 1)

# Degrees of freedom for emmeans contrasts. Kenward-Roger is the emmeans
# default when pbkrtest is installed, and was used in Script_6.
emm_options(lmer.df = "kenward-roger")


# Step 1: Read intensities and sample information ------------------------------

# data7: imputed, normalized log2 intensities from Script_3
# (3,600 proteins x 96 samples, columns S1-S96, no NA), as in Script_6.
# info: sample annotation (Sample.ID, Sample.Group, Sample.PID, Sample.Age,
# Sample.Region, Sample.Sex), rows in S1-S96 order.
e <- new.env()
load(file.path(indir, "UW.Madison.rat.lung.proteomics_data_tables.RData"), envir = e)
data7 <- as.matrix(e$data_table$data_imp)
info  <- e$info_table
stopifnot(identical(colnames(data7), info$Sample.ID))


# Step 2: Add technical variables and build the model table --------------------

# Read preparation day, within-day preparation order and injection position
# for each sample, and match them to the S1-S96 order of data7.
meta <- read_excel(file.path(indir, "sample_metadata.xlsx"), sheet = 1)
meta <- transmute(meta,
  Sample.ID      = as.character(`Sample ID`),
  PreparationDay = as.integer(`Preparation Day`),
  WithinDayOrder = as.numeric(`Within-day Preparation Order`),
  InjectionOrder = as.numeric(`Injection Order`)
)
m <- match(info$Sample.ID, meta$Sample.ID)
stopifnot(!anyNA(m))

# workdata: one row per sample, in the column order of data7.
# Age levels Young, Aged so that contrasts are Aged - Young (A = 7-11m,
# Y = 7w). PID = animal (random intercept). PreparationDay is a factor
# (4 levels); InjectionOrder is continuous.
workdata <- data.frame(
  Sample.ID      = info$Sample.ID,
  Age            = factor(info$Sample.Age, levels = c("Young", "Aged")),
  Sex            = factor(info$Sample.Sex),
  Region         = factor(info$Sample.Region),
  PID            = factor(info$Sample.PID),
  PreparationDay = factor(meta$PreparationDay[m], levels = 1:4),
  WithinDayOrder = meta$WithinDayOrder[m],
  InjectionOrder = meta$InjectionOrder[m]
)

# Centre and scale the order variables (mean 0, SD 1); the injection-order
# coefficient is then the change in log2 abundance per SD of injection
# position.
workdata$WithinDayOrderScaled <- as.numeric(scale(workdata$WithinDayOrder))
workdata$InjectionOrderScaled <- as.numeric(scale(workdata$InjectionOrder))
stopifnot(!anyNA(workdata), setequal(workdata$InjectionOrder, 1:96))

# Observed-value counts per protein and group, for the Method 4 filter
# (Step 6). Columns are group codes, e.g. ABF = Aged Bottom Female.
counts <- as.data.frame(read_excel(file.path(indir, "Real_values_number_in_sample_group.xlsx")))
rownames(counts) <- counts$Uniprot
counts$Uniprot <- NULL
stopifnot(all(rownames(data7) %in% rownames(counts)))


# Step 3: Design balance and collinearity --------------------------------------

# For each biological factor, test association with preparation day
# (chi-square and Fisher exact test on the factor x day table) and whether
# injection order differs between factor levels (Kruskal-Wallis).
# Tests are at sample level (96 samples, 4 per animal). The Fisher test uses
# Monte Carlo p-values (1e5 replicates, seed 1).
bal <- bind_rows(lapply(c("Age", "Sex", "Region"), function(v) {
  tab <- table(workdata[[v]], workdata$PreparationDay)
  data.frame(
    Variable = v,
    ChiSq_vs_PreparationDay_p      = suppressWarnings(chisq.test(tab)$p.value),
    Fisher_vs_PreparationDay_p     = { set.seed(1); fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value },
    KruskalWallis_InjectionOrder_p = kruskal.test(workdata$InjectionOrder ~ workdata[[v]])$p.value
  )
}))

# Sample counts per factor level and preparation day.
crosstabs <- bind_rows(lapply(c("Age", "Sex", "Region"), function(v) {
  t <- as.data.frame.matrix(table(workdata[[v]], workdata$PreparationDay))
  colnames(t) <- paste0("Day", colnames(t))
  cbind(Variable = v, Level = rownames(t), t)
}))

# Spearman correlations among the three technical variables
# (preparation day and injection order are strongly correlated, rho ~ 0.93,
# because acquisition followed preparation-day blocks).
tech_cor <- cor(cbind(PreparationDay = as.numeric(workdata$PreparationDay),
                      WithinDayOrder = workdata$WithinDayOrder,
                      InjectionOrder = workdata$InjectionOrder), method = "spearman")

# Variance inflation factors of the main-effect design (biological +
# technical). VIF ~ 1 for biological terms means they are not collinear
# with the technical covariates.
X <- model.matrix(~ Age + Sex + Region + PreparationDay + InjectionOrderScaled, workdata)[, -1]
vif <- sapply(colnames(X), function(cn) {
  1 / (1 - summary(lm(X[, cn] ~ X[, colnames(X) != cn]))$r.squared)
})

# Model 1 fixed-effect design must be full rank (no perfect confounding).
design_batch <- model.matrix(~ Age + Sex + Region + Age:Region + Sex:Region + Age:Sex +
                               PreparationDay + InjectionOrderScaled, workdata)
stopifnot(qr(design_batch)$rank == ncol(design_batch))


# Step 4: PCA and association of PCs with each factor --------------------------

# PCA of the imputed matrix with samples as observations, centred and not
# scaled (same as prcomp(t(data7)) in Script_3).
pc <- prcomp(t(data7), center = TRUE, scale. = FALSE)
ve <- summary(pc)$importance[2, 1:10]   # proportion of variance, PC1-PC10
pcs <- as.data.frame(pc$x[, 1:10])
pcs$Sample.ID <- rownames(pcs)
pcs <- left_join(pcs, workdata, by = "Sample.ID")

# For each PC, R2 of a linear model PC ~ factor (one factor at a time) and,
# for the technical factors, the F-test p-value.
pc_assoc <- bind_rows(lapply(paste0("PC", 1:10), function(p) {
  r2 <- function(f) summary(lm(as.formula(paste(p, "~", f)), pcs))$r.squared
  pv <- function(f) anova(lm(as.formula(paste(p, "~", f)), pcs))[1, "Pr(>F)"]
  data.frame(PC = p, VarExplained = ve[p],
    R2_Age = r2("Age"), R2_Sex = r2("Sex"), R2_Region = r2("Region"), R2_PID = r2("PID"),
    R2_PreparationDay = r2("PreparationDay"), p_PreparationDay = pv("PreparationDay"),
    R2_InjectionOrder = r2("InjectionOrder"), p_InjectionOrder = pv("InjectionOrder"))
}))


# Step 5: Per-protein models ----------------------------------------------------

# Model formulas.
# f0: original model (Script_6). f1: f0 + technical covariates.
# fvp: variance-partition model; all categorical factors as random
# intercepts (as in the variancePartition package), injection order fixed.
# See the caveats in fit_one() below.
f0  <- Protein ~ Age + Sex + Region + Age:Region + Sex:Region + Age:Sex + (1 | PID)
f1  <- update(f0, . ~ . + PreparationDay + InjectionOrderScaled)
fvp <- Protein ~ (1 | Age) + (1 | Sex) + (1 | Region) + (1 | PID) + (1 | PreparationDay) + InjectionOrderScaled

# Aged - Young contrast within each Region x Sex (8 rows per protein).
# Equivalent to emmeans(trt.vs.ctrl ~ Age | Region + Sex) in Script_6.
age_contr <- function(fit) {
  emm <- emmeans(fit, ~ Age | Region + Sex)
  as.data.frame(summary(contrast(emm, method = list("Aged - Young" = c(-1, 1)))))
}

# Fit all models for protein i and return:
#   r0 / r1: age contrasts from Model 0 / Model 1
#   tech:    LRT and F tests of the technical terms
#   vp:      variance fractions
fit_one <- function(i) {
  d <- workdata
  d$Protein <- as.numeric(data7[i, ])
  ctrl <- lmerControl(calc.derivs = FALSE)

  # Model 0 and Model 1 (REML), used for the contrasts.
  m0 <- lmerTest::lmer(f0, d, control = ctrl)
  m1 <- lmerTest::lmer(f1, d, control = ctrl)

  # Likelihood-ratio test of the technical block: both models refitted by
  # maximum likelihood (REML likelihoods are not comparable when fixed
  # effects differ). Chi-square with 4 df (PreparationDay 3 + InjectionOrder 1).
  a1 <- anova(update(m0, REML = FALSE), update(m1, REML = FALSE))

  # Type III F tests (Satterthwaite) of each technical term in Model 1.
  an <- anova(m1, type = 3)

  # Variance partitioning: random-effect variances, plus the variance of the
  # fixed-effect predictions for injection order, as fractions of the total.
  # The 7 fractions (Age, Sex, Region, PID, PreparationDay, InjectionOrder,
  # Residual) sum to 1 for each protein.
  #
  # Caveats (this is an approximation; interpret medians, not single proteins):
  #   - Age and Sex have only 2 levels and PreparationDay 4. A variance
  #     estimated from so few levels is imprecise and is often estimated as 0
  #     (boundary fit), so per-protein values are noisy and medians are lower
  #     than means.
  #   - The model has no interaction terms, so the biological fractions are
  #     marginal and differ from the effects in Model 0 / Model 1.
  #   - Injection order is a fixed effect. Its fraction is the variance of
  #     its fitted values (coefficient^2 x variance of the scaled order),
  #     added to the random-effect variances. var() removes the intercept, so
  #     the intercept does not contribute; variancePartition uses the same
  #     approach for fixed effects.

  vp <- tryCatch({
    mv <- lmer(fvp, d, REML = FALSE, control = ctrl)
    vc <- as.data.frame(VarCorr(mv))
    v  <- setNames(vc$vcov, vc$grp)
    v["InjectionOrder"] <- var(predict(mv, re.form = NA))
    v / sum(v)
  }, error = function(err) rep(NA_real_, 7))

  id <- rownames(data7)[i]
  r0 <- age_contr(m0)
  r0$Uniprot.ID <- id
  r1 <- age_contr(m1)
  r1$Uniprot.ID <- id
  tech <- data.frame(Uniprot.ID = id,
    LRT_chisq = a1$Chisq[2], LRT_df = a1$Df[2], LRT_p = a1$`Pr(>Chisq)`[2],
    F_PreparationDay_p = an["PreparationDay", "Pr(>F)"],
    F_InjectionOrder_p = an["InjectionOrderScaled", "Pr(>F)"],
    Coef_InjectionOrder_perSD = fixef(m1)["InjectionOrderScaled"],
    singular0 = isSingular(m0),   # TRUE when the animal variance is estimated as 0
    singular1 = isSingular(m1))
  vpd <- as.data.frame(t(vp))
  vpd$Uniprot.ID <- id
  list(r0 = r0, r1 = r1, tech = tech, vp = vpd)
}

# Fit all proteins in parallel. "boundary (singular) fit" messages are
# expected (animal variance estimated as 0 for some proteins) and do not
# affect the contrasts. Proteins whose fit fails are reported and dropped.
idx <- seq_len(nrow(data7))
if (!is.na(n_test)) idx <- idx[1:n_test]
t0 <- Sys.time()
res <- mclapply(idx, function(i) {
  tryCatch(fit_one(i), error = function(err) {
    message(i, ": ", conditionMessage(err))
    NULL
  })
}, mc.cores = ncores)
cat("Fitting time:", format(Sys.time() - t0), "\n")
ok <- !sapply(res, is.null)
cat("Failed:", sum(!ok), "\n")
res <- res[ok]

# Combine the per-protein results into tables.
grab  <- function(k) bind_rows(lapply(res, `[[`, k))
age0  <- grab("r0")     # Model 0 contrasts (proteins x 8 comparisons)
age1  <- grab("r1")     # Model 1 contrasts
tech  <- grab("tech")   # technical tests, one row per protein
vpart <- grab("vp")     # variance fractions, one row per protein

# BH adjustment of the technical tests across all 3,600 proteins.
tech$LRT_FDR              <- p.adjust(tech$LRT_p, "BH")
tech$F_PreparationDay_FDR <- p.adjust(tech$F_PreparationDay_p, "BH")
tech$F_InjectionOrder_FDR <- p.adjust(tech$F_InjectionOrder_p, "BH")


# Step 6: Method 4 filter, FDR and DEP calls (as in Script_6) ------------------

# For each Region x Sex comparison: keep proteins with >= 5 observed values
# in the aged or the young group (group codes A/Y + region initial + sex
# initial, e.g. ABF vs YBF), BH-adjust the contrast p-values within the
# comparison, and call DEPs (FDR < 0.05 and |log2FC| > log2(1.5)).
filter_and_adjust <- function(df) {
  df %>%
    group_by(Region, Sex) %>%
    group_modify(~ {
      ag <- paste0("A", substr(.y$Region, 1, 1), substr(.y$Sex, 1, 1))
      yg <- paste0("Y", substr(.y$Region, 1, 1), substr(.y$Sex, 1, 1))
      stopifnot(ag %in% colnames(counts), yg %in% colnames(counts))
      keep <- counts[.x$Uniprot.ID, ag] >= 5 | counts[.x$Uniprot.ID, yg] >= 5
      out <- .x[keep, , drop = FALSE]
      out$p.value.adjusted <- p.adjust(out$p.value, "BH")
      out
    }) %>%
    ungroup() %>%
    mutate(DEP = !is.na(p.value.adjusted) & p.value.adjusted < 0.05 & abs(estimate) > log2(1.5))
}
age0 <- filter_and_adjust(age0)
age1 <- filter_and_adjust(age1)


# Step 7: Compare Model 0 and Model 1 ------------------------------------------

# cmp: one row per protein x comparison (26,281 rows), with estimate, SE,
# FDR and DEP from Model 0 (suffix 0) and Model 1 (suffix 1).
sel <- function(df, s) {
  df %>%
    select(Uniprot.ID, Region, Sex, estimate, SE, p.value.adjusted, DEP) %>%
    rename_with(~ paste0(.x, s), c(estimate, SE, p.value.adjusted, DEP))
}
cmp <- sel(age0, "0") %>%
  inner_join(sel(age1, "1"), by = c("Uniprot.ID", "Region", "Sex")) %>%
  mutate(Comparison = paste0(Region, "-", Sex, ": Aged vs Young"),
         Diff1 = estimate1 - estimate0)

# Robustness summary per comparison:
#   Shared / Lost / Gained: DEP calls in both / Model 0 only / Model 1 only
#   Direction_concordance: % of Model 0 DEPs with the same sign in Model 1
#   Pearson_r: correlation of log2FC between models (all proteins)
#   Median_SE_ratio: SE(Model 1) / SE(Model 0); < 1 = more precise
summ_fun <- function(d, e, D, lab) {
  d %>%
    group_by(Comparison) %>%
    summarise(
      Model = lab, N_proteins = n(),
      DEPs_Model0 = sum(DEP0), DEPs_Adjusted = sum(.data[[D]]),
      Shared = sum(DEP0 & .data[[D]]), Lost = sum(DEP0 & !.data[[D]]), Gained = sum(!DEP0 & .data[[D]]),
      Pct_Model0_DEPs_retained = 100 * Shared / pmax(DEPs_Model0, 1),
      Direction_concordance_DEPs_pct = 100 * mean(sign(estimate0[DEP0]) == sign(.data[[e]][DEP0])),
      Pearson_r = cor(estimate0, .data[[e]]),
      Median_abs_diff_log2FC = median(abs(.data[[e]] - estimate0)),
      Median_SE_ratio = median(.data[[sub("estimate", "SE", e)]] / SE0),
      .groups = "drop")
}
robust <- summ_fun(cmp, "estimate1", "DEP1", "Model 1: + PreparationDay + InjectionOrder")

# Pooled over the 8 comparisons.
overall <- data.frame(
  Model = "Model 1", Comparisons = 8,
  Total_DEPs_Model0   = sum(cmp$DEP0),
  Total_DEPs_Adjusted = sum(cmp$DEP1),
  Shared              = sum(cmp$DEP0 & cmp$DEP1),
  Pct_retained        = 100 * sum(cmp$DEP0 & cmp$DEP1) / sum(cmp$DEP0),
  Pearson_r_all       = cor(cmp$estimate0, cmp$estimate1))

# The 38 matrisome proteins of Figure 3h (checked at the end of the script).
mat38 <- tryCatch(
  read.delim(file.path(indir, "38 matrisome DEPs with unique features in 8 pairwise comparisons.txt")),
  error = function(e) NULL)

# Number of proteins with a significant technical effect (FDR < 0.05).
tech_summary <- data.frame(
  Test = c("Joint LRT: PreparationDay + InjectionOrder (4 df)",
           "PreparationDay (Type III F, Model 1)",
           "InjectionOrder (Type III F, Model 1)"),
  N_proteins = nrow(tech),
  N_FDR_lt_0.05 = c(sum(tech$LRT_FDR < .05),
                    sum(tech$F_PreparationDay_FDR < .05),
                    sum(tech$F_InjectionOrder_FDR < .05)))
tech_summary$Pct <- 100 * tech_summary$N_FDR_lt_0.05 / tech_summary$N_proteins

# Median and mean variance fraction of each component across proteins.
vp_components <- setdiff(names(vpart), "Uniprot.ID")
vp_summary <- data.frame(
  Component  = vp_components,
  Median_pct = sapply(vp_components, function(k) 100 * median(vpart[[k]], na.rm = TRUE)),
  Mean_pct   = sapply(vp_components, function(k) 100 * mean(vpart[[k]], na.rm = TRUE)))


# Step 8: Diagnostic plots and exports -----------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Diagnostic plots (the publication figure is made by Script_02).
pdf(file.path(outdir, "Batch_effect_assessment_figures.pdf"), width = 11, height = 8)

# PCA coloured by preparation day, shapes = region.
print(ggplot(pcs, aes(PC1, PC2, colour = PreparationDay, shape = Region)) +
  geom_point(size = 2.5) + theme_bw() +
  labs(title = "PCA (all 3,600 proteins) coloured by preparation day",
       x = sprintf("PC1 (%.1f%%)", 100 * ve[1]), y = sprintf("PC2 (%.1f%%)", 100 * ve[2])))

# PCA coloured by injection order, shapes = age.
print(ggplot(pcs, aes(PC1, PC2, colour = InjectionOrder, shape = Age)) +
  geom_point(size = 2.5) + theme_bw() + scale_colour_viridis_c() +
  labs(title = "PCA coloured by LC-MS/MS injection order",
       x = sprintf("PC1 (%.1f%%)", 100 * ve[1]), y = sprintf("PC2 (%.1f%%)", 100 * ve[2])))

# PC1 against injection order (drift check).
print(ggplot(pcs, aes(InjectionOrder, PC1, colour = PreparationDay)) +
  geom_point() + theme_bw() + labs(title = "PC1 versus injection order"))

# Variance fractions per protein, components ordered by median.
vlong <- tidyr::pivot_longer(vpart, -Uniprot.ID, names_to = "Component", values_to = "Fraction")
print(ggplot(vlong, aes(reorder(Component, -Fraction, median, na.rm = TRUE), 100 * Fraction)) +
  geom_violin(fill = "grey85") + geom_boxplot(width = .12, outlier.size = .3) + theme_bw() +
  labs(x = NULL, y = "% variance explained", title = "Variance partitioning per protein"))

# Model 0 vs Model 1 Aged - Young log2FC; red = Model 0 DEP.
print(ggplot(cmp, aes(estimate0, estimate1, colour = DEP0)) +
  geom_point(alpha = .35, size = .7) + geom_abline(linetype = 2, colour = "red") +
  facet_grid(Sex ~ Region) + theme_bw() +
  scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "firebrick")) +
  labs(x = "Model 0 log2FC (Aged - Young)", y = "Model 1 (batch-adjusted) log2FC", colour = "Model 0 DEP",
       title = "Robustness of age effects to technical covariate adjustment"))

# LRT p-value distribution (a peak near 0 = technical effects present).
hist(tech$LRT_p, breaks = 50, main = "Joint LRT p-values for technical covariates", xlab = "p")
dev.off()

# All tables in one workbook, one worksheet each.
write.xlsx(list(
  README = data.frame(
    Item = c("Model 0", "Model 1", "Data", "DEP rule", "df method"),
    Description = c(deparse1(f0), deparse1(f1),
                    "data_imp (3600 proteins x 96 samples, log2, imputed)",
                    "Method 4 filter (>=5 observed in either group), BH within comparison, FDR<0.05 & |log2FC|>log2(1.5)",
                    "Kenward-Roger (emmeans)")),
  Design_balance = bal,
  Design_crosstabs = crosstabs,
  Technical_cor_spearman = cbind(Variable = rownames(tech_cor), as.data.frame(tech_cor)),
  VIF = data.frame(Term = names(vif), VIF = vif),
  PCA_association = pc_assoc,
  Variance_partition_summary = vp_summary,
  Technical_tests_summary = tech_summary,
  Robustness_summary = robust,
  Robustness_overall = overall,
  Protein_comparison = cmp,
  Technical_tests_per_protein = tech,
  Variance_partition_per_protein = vpart),
  file.path(outdir, "Batch_effect_sensitivity_analysis.xlsx"), overwrite = TRUE)

# Key objects for Script_02 (figure) and Script_03 (table).
save(cmp, tech, vpart, pc_assoc, robust, overall, workdata,
     file = file.path(outdir, "batch_results.RData"))


# Step 9: Print key results -----------------------------------------------------

# Expected in the reported run: 1,448 Model 0 DEP calls, 1,336 retained
# (92.3%), pooled r = 0.984, 2,626 proteins with LRT FDR < 0.05, and
# 106 of 111 matrisome-38 DEP calls retained.
print(bal)
print(round(vif, 2))
print(pc_assoc[1:5, ], digits = 3)
print(vp_summary, digits = 3)
print(tech_summary)
print(as.data.frame(robust), digits = 3)
print(overall, digits = 3)
if (!is.null(mat38)) {
  z <- cmp %>% filter(Uniprot.ID %in% mat38$Uniprot.ID, DEP0)
  cat("Matrisome-38 DEP calls in Model 0:", nrow(z), " retained in Model 1:", sum(z$DEP1), "\n")
}


# Save the full R workspace, including intermediate objects.
save.image(file.path(outdir, "01_batch_effect_models_workspace.RData"))
