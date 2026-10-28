# =============================================================================
# Script 02: Supplementary figure - batch-effect assessment
# =============================================================================
#
# Purpose:
#   Draw the four-panel supplementary figure that accompanies the response to
#   Reviewer comment 12, from the results of Script_01:
#     a  PCA, PC1 vs PC2, coloured by preparation day (shapes = age)
#     b  PCA, PC1 vs PC3, coloured by preparation day (shapes = age)
#        (PC3 is the component most associated with preparation day)
#     c  Variance partitioning: per-protein % variance explained by
#        preparation day, injection order, age, region, sex, animal and
#        residual (3,600 proteins; labels = medians)
#     d  Aged - Young log2 fold changes from the original model (x) and the
#        batch-adjusted model (y) for the 8 Region x Sex comparisons; red =
#        DEP in the original model; each panel shows Pearson r and the number
#        of original DEPs retained.
#
# PCA settings (as in Script_3, Figures 2a / 5a / 6a):
#   ggbiplot(prcomp(t(data7)), obs.scale = 1, var.scale = 1, ellipse = TRUE,
#            circle = TRUE, var.axes = FALSE)
#   prcomp() centres and does not scale; ellipses are 68% normal-probability
#   ellipses per group. Here groups = preparation day; open shapes = age.
#
# Inputs:
#   - results/batch_results.RData  (output of Script_01)
#       cmp (Model 0 vs Model 1 contrasts), vpart (variance fractions),
#       workdata (sample table with PreparationDay and Age).
#   - data/UW.Madison.rat.lung.proteomics_data_tables.RData  (Script_3)
#       data_table$data_imp for the PCA.
#
# Outputs:
#   - figures/Supplementary_Figure_batch_effect.pdf / .png
#       11 x 14 in; png at 300 dpi.
#   - figures/PCA_by_preparation_day.pdf
#       Script_3-style PCA, one view per page (PC1 vs PC2, PC1 vs PC3).
#   - results/02_supplementary_figure_workspace.RData
#       Full R workspace.
#
# How to run:
#   After Script_01, from the batch_effect_analysis folder:
#       Rscript code/02_supplementary_figure.R
#   or Source it in RStudio. Optional arguments: [indir] [resultsdir] [outdir].
#   Runs in a few seconds.
# =============================================================================

# Load required packages.
# ggbiplot (vqv/ggbiplot, the version used in Script_3) attaches plyr. Load it
# BEFORE dplyr, otherwise plyr::summarise / plyr::mutate mask the dplyr
# functions used below and the script fails.
# ggplot2: plots; dplyr / tidyr: table handling; patchwork: panel layout.
suppressPackageStartupMessages({
  library(ggbiplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
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
outdir <- if (length(args) >= 3) args[3] else "figures"   # figure output folder


# Step 1: Read results and intensities ------------------------------------------

# Script_01 results: cmp, tech, vpart, pc_assoc, robust, overall, workdata.
load(file.path(resdir, "batch_results.RData"))

# data7: imputed log2 intensities (3,600 proteins x S1-S96), for the PCA.
e <- new.env()
load(file.path(indir, "UW.Madison.rat.lung.proteomics_data_tables.RData"), envir = e)
data7 <- as.matrix(e$data_table$data_imp)
stopifnot(identical(colnames(data7), workdata$Sample.ID))


# Step 2: Colours and plot theme ------------------------------------------------

# Preparation-day colours (colour-blind-checked categorical palette).
# Text uses neutral ink colours, not the series colours.
day_cols <- c(`1` = "#2a78d6", `2` = "#eb6834", `3` = "#1baf7a", `4` = "#eda100")
ink  <- "#0b0b0b"   # primary text
ink2 <- "#52514e"   # secondary text, axes

# Shared theme for panels c and d (panels a and b use theme_bw, as in Script_3).
th <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(colour = ink2, linewidth = .3),
        axis.ticks = element_line(colour = ink2, linewidth = .3),
        axis.text = element_text(colour = ink2),
        axis.title = element_text(colour = ink),
        plot.title = element_text(face = "bold", size = 9.5, colour = ink),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        plot.tag = element_text(face = "bold", size = 12))


# Step 3: Panels a and b - PCA by preparation day ------------------------------

# PCA of the imputed matrix with samples as observations (as in Script_3).
pca <- prcomp(t(data7))

# Grouping (colour + ellipse) = preparation day; shape = age.
Day     <- factor(workdata$PreparationDay)
Age_lab <- factor(ifelse(workdata$Age == "Aged", "7-11m", "7w"), levels = c("7w", "7-11m"))

# ggbiplot draws filled points coloured by group; black open shapes for age
# are drawn on top. ch = the two PCs to plot.
pca_plot <- function(ch, ttl) {
  g <- ggbiplot(pca, choices = ch, var.axes = FALSE, labels = NULL, groups = Day,
                obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
  g + geom_point(aes(shape = Age_lab), colour = "black", size = 1.6) +
    scale_shape_manual(values = c(`1` = 16, `2` = 16, `3` = 16, `4` = 16, `7w` = 1, `7-11m` = 2),
                       breaks = c("7w", "7-11m"), name = "Age") +
    guides(fill = "none",
           colour = guide_legend(override.aes = list(shape = 16, size = 3, fill = NA, linetype = 0))) +
    labs(colour = "Preparation day", title = ttl) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9.5),
          plot.tag = element_text(face = "bold", size = 12))
}
pa <- pca_plot(c(1, 2), "PCA, PC1 vs PC2")
pb <- pca_plot(c(1, 3), "PCA, PC1 vs PC3")


# Step 4: Panel c - variance partitioning --------------------------------------

# Long format: one row per protein x component. Technical components
# (preparation day, injection order) are filled orange; biological
# components and the residual blue.
lab <- c(PreparationDay = "Prep.\nday", InjectionOrder = "Injection\norder", Age = "Age",
         Region = "Region", Sex = "Sex", PID = "Animal", Residual = "Resid.")
vl <- vpart %>%
  pivot_longer(-Uniprot.ID, names_to = "Component", values_to = "Fraction") %>%
  filter(!is.na(Fraction)) %>%
  mutate(Technical = ifelse(Component %in% c("PreparationDay", "InjectionOrder"),
                            "Technical", "Biological / residual"),
         Component = factor(lab[Component], levels = lab))

# Median % variance per component, printed above each violin.
med <- vl %>% group_by(Component) %>% summarise(m = 100 * median(Fraction))

pc_ <- ggplot(vl, aes(Component, 100 * Fraction, fill = Technical)) +
  geom_violin(scale = "width", colour = NA, alpha = .55) +
  geom_boxplot(width = .14, outlier.shape = NA, fill = "white", colour = ink2, linewidth = .3) +
  geom_text(data = med, aes(Component, 105, label = sprintf("%.1f%%", m)),
            inherit.aes = FALSE, size = 2.6, colour = ink2) +
  scale_fill_manual(values = c(Technical = "#eb6834", `Biological / residual` = "#2a78d6"), name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Variance explained per protein (%)", title = "Variance partitioning") +
  th + theme(legend.position = "top", legend.justification = "left")


# Step 5: Panel d - age effects, original vs adjusted model --------------------

# Facet order: 4 regions x Female, then 4 regions x Male.
cmp <- cmp %>%
  mutate(Comp = factor(paste(Region, Sex),
                       levels = paste(rep(c("Bottom", "Central", "Peripheral", "Top"), 2),
                                      rep(c("Female", "Male"), each = 4))))

# Per-panel labels: Pearson r (all proteins, unclipped values) and the number
# of original-model DEPs that remain DEPs in the adjusted model.
ann <- cmp %>%
  group_by(Comp) %>%
  summarise(r = cor(estimate0, estimate1), n0 = sum(DEP0), kept = sum(DEP0 & DEP1)) %>%
  mutate(lab = sprintf("r = %.3f\n%d/%d DEPs kept", r, kept, n0))

# Axis range fixed to +/-4 log2FC; the few points beyond are drawn on the
# panel border (number reported in the caption).
lim   <- c(-4, 4)
n_out <- sum(abs(cmp$estimate0) > 4 | abs(cmp$estimate1) > 4)
cmp   <- cmp %>% mutate(estimate0 = pmax(pmin(estimate0, 4), -4),
                        estimate1 = pmax(pmin(estimate1, 4), -4))

# Grey = not DEP, red = DEP in the original model; dashed line y = x.
pd <- ggplot(cmp, aes(estimate0, estimate1)) +
  geom_abline(slope = 1, intercept = 0, colour = ink2, linewidth = .3, linetype = 2) +
  geom_point(data = filter(cmp, !DEP0), colour = "#bdbcb6", size = .45, stroke = 0) +
  geom_point(data = filter(cmp, DEP0), colour = "#d03b3b", size = .7, stroke = 0) +
  geom_text(data = ann, aes(x = lim[1], y = lim[2], label = lab),
            hjust = 0, vjust = 1, size = 2.4, colour = ink, lineheight = .9) +
  facet_wrap(~ Comp, nrow = 2) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression("Original model log"[2]*"FC (Aged - Young)"),
       y = expression("Adjusted model log"[2]*"FC"),
       title = "Age effects before and after adjustment for preparation day and injection order",
       caption = sprintf(paste("Red: differentially abundant in the original model (FDR < 0.05, |FC| > 1.5).",
                               "Dashed line: y = x. %d protein-comparison pairs with |log2FC| > 4 are drawn",
                               "on the panel border. r computed on unclipped values."), n_out)) +
  th + theme(plot.caption = element_text(colour = ink2, size = 7.5, hjust = 0))


# Step 6: Assemble and save ----------------------------------------------------

# Row 1: panels a and b with one shared legend; row 2: panel c; row 3: panel d.
row1 <- (pa | pb) + plot_layout(guides = "collect")
fig  <- row1 / pc_ / pd +
  plot_layout(heights = c(1.05, 0.75, 1.45)) +
  plot_annotation(tag_levels = "a")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "Supplementary_Figure_batch_effect.pdf"), fig, width = 11, height = 14)   # default pdf device (cairo_pdf needs XQuartz on macOS)
ggsave(file.path(outdir, "Supplementary_Figure_batch_effect.png"), fig, width = 11, height = 14, dpi = 300)

# Script_3-style PCA (larger font), one page per view.
pdf(file.path(outdir, "PCA_by_preparation_day.pdf"))
for (ch in list(c(1, 2), c(1, 3))) plot(pca_plot(ch, NULL) + theme_bw(base_size = 16))
dev.off()
cat("saved\n")


# Save the full R workspace, including intermediate objects.
save.image(file.path(resdir, "02_supplementary_figure_workspace.RData"))

