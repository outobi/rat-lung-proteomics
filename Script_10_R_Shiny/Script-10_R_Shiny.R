# Load required packages.
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(tidyr)
library(rio)
library(readxl)
library(ggplot2)
library(gridExtra)


# Retain dependencies used by the later, unchanged analysis sections.
library(limma)
library(reshape2)
library(Rmisc)
library(stringr)
library(multcomp)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load the saved imputed matrices and sample information.
load("Script_6_multiple_factor_DE.RData")


# Prepare comparison tables for single-gene queries R Shiny app.
# Fill unfiltered result tables with Method 4 FDR vaues where available.

# Prepare age versus young results.
table <- list()
for (j in 1:8) {
  bbb <- age_comparisons_unfiltered[[j]]
  bbb$p.value.adjusted <- NA

  i = 1
  ccc <- age.vs.young4[[j]]
  for (items in ccc$Uniprot.ID) {
    if (items %in% bbb$Uniprot.ID) {
      index <- which(bbb$Uniprot.ID == items)
      bbb$p.value.adjusted[index] <- ccc$p.value.adjusted[i]
      i = i + 1
    }
  }

  table[[j]] <- bbb
}

names(table) <- names(age.vs.young4)

table1 <- list()
table1[c(1, 2, 3, 4, 5, 6, 7, 8)] <- table[c(5, 1, 8, 4, 6, 2, 7, 3)]

age_vs_young_index_all <- c(
  "Bottom Male", "Bottom Female", "Top Male", "Top Female",
  "Central Male", "Central Female", "Peripheral Male", "Peripheral Female"
)
names(table1) <- age_vs_young_index_all
table <- table1
# Combine the eight aged-versus-young comparison tables.

df1 <- cbind(
  df0,
  table[[1]][c(15, 20)],
  table[[2]][c(15, 20)],
  table[[3]][c(15, 20)],
  table[[4]][c(15, 20)],
  table[[5]][c(15, 20)],
  table[[6]][c(15, 20)],
  table[[7]][c(15, 20)],
  table[[8]][c(15, 20)]
)
df1 <- df1[df1$Gene.symbol != "-", ]
rownames(df1) <- df1$Gene.symbol
df2 <- df1[c(8:23)]

df3 <- df2[c(1, 3, 5, 7, 9, 11, 13, 15, 2, 4, 6, 8, 10, 12, 14, 16)]
df4 <- as.data.frame(t(df3))
df4 <- signif(df4, 4)


# Prepare the reference table with matrisome annotations.
df6 <- df0

df6 <- df6 %>% arrange(Gene.symbol)

save(df1, df4, df6, file = "helper age.RData")


# Prepare female-versus-male results.

table <- list()
for (j in 1:8) {
  bbb <- sex_comparisons_unfiltered[[j]]
  bbb$p.value.adjusted <- NA

  i = 1
  ccc <- female.vs.male4[[j]]
  for (items in ccc$Uniprot.ID) {
    if (items %in% bbb$Uniprot.ID) {
      index <- which(bbb$Uniprot.ID == items)
      bbb$p.value.adjusted[index] <- ccc$p.value.adjusted[i]
      i = i + 1
    }
  }

  table[[j]] <- bbb
}

names(table) <- names(female.vs.male4)

table1 <- list()
table1[c(1, 2, 3, 4, 5, 6, 7, 8)] <- table[c(1, 5, 4, 8, 2, 6, 3, 7)]

female.vs.male.index.all <- c(
  "Bottom Age",
  "Bottom Young",
  "Top Age",
  "Top Young",
  "Central Age",
  "Central Young",
  "Peripheral Age",
  "Peripheral Young"
)

names(table1) <- female.vs.male.index.all
table <- table1


# Combine the eight female versus male comparison tables.

df1 <- cbind(
  df0,
  table[[1]][c(15, 20)],
  table[[2]][c(15, 20)],
  table[[3]][c(15, 20)],
  table[[4]][c(15, 20)],
  table[[5]][c(15, 20)],
  table[[6]][c(15, 20)],
  table[[7]][c(15, 20)],
  table[[8]][c(15, 20)]
)
df1 <- df1[df1$Gene.symbol != "-", ]
rownames(df1) <- df1$Gene.symbol
df2 <- df1[c(8:23)]

df3 <- df2[c(1, 3, 5, 7, 9, 11, 13, 15, 2, 4, 6, 8, 10, 12, 14, 16)]
df4 <- as.data.frame(t(df3))
df4 <- signif(df4, 4)


# Prepare the reference table with matrisome annotations.
df6 <- df0
df6 <- df6 %>% arrange(Gene.symbol)

save(df1, df4, df6, file = "helper sex.RData")


# Prepare region-versus-region results.

table <- list()
for (j in 1:24) {
  bbb <- region_comparisons_unfiltered[[j]]
  bbb$p.value.adjusted <- NA

  i = 1
  ccc <- region.vs.region4[[j]]
  for (items in ccc$Uniprot.ID) {
    if (items %in% bbb$Uniprot.ID) {
      index <- which(bbb$Uniprot.ID == items)
      bbb$p.value.adjusted[index] <- ccc$p.value.adjusted[i]
      i = i + 1
    }
  }

  table[[j]] <- bbb
}
names(table) <- names(region.vs.region4)


# Combine the 24 region versus region comparison tables.

# Order comparisons by young male, young female, aged male, and aged female.
table1 <- table$bottom.vs.central.male.young

table2 <- table$bottom.vs.peri.male.young

table3 <- table$bottom.vs.top.male.young

table4 <- table$central.vs.peri.male.young

table5 <- table$central.vs.top.male.young

table6 <- table$peri.vs.top.male.young

table7 <- table[[13]]
table8 <- table[[14]]
table9 <- table[[15]]
table10 <- table[[16]]
table11 <- table[[17]]
table12 <- table[[18]]

table13 <- table[[7]]
table14 <- table[[8]]
table15 <- table[[9]]
table16 <- table[[10]]
table17 <- table[[11]]
table18 <- table[[12]]

table19 <- table[[1]]
table20 <- table[[2]]
table21 <- table[[3]]
table22 <- table[[4]]
table23 <- table[[5]]
table24 <- table[[6]]

df1 <- cbind(
  df0, table1[c(15, 20)], table2[c(15, 20)], table3[c(15, 20)], table4[c(15, 20)], table5[c(15, 20)], table6[c(15, 20)], table7[c(15, 20)], table8[c(15, 20)],
  table9[c(15, 20)], table10[c(15, 20)], table11[c(15, 20)], table12[c(15, 20)], table13[c(15, 20)], table14[c(15, 20)], table15[c(15, 20)], table16[c(15, 20)],
  table17[c(15, 20)], table18[c(15, 20)], table19[c(15, 20)], table20[c(15, 20)], table21[c(15, 20)], table22[c(15, 20)], table23[c(15, 20)], table24[c(15, 20)]
)

df1 <- df1[df1$Gene.symbol != "-", ]
rownames(df1) <- df1$Gene.symbol
df2 <- df1[c(8:55)]

df3 <- df2[c(seq(1, 48, 2), seq(2, 48, 2))]
df4 <- as.data.frame(t(df3))
df4 <- signif(df4, 4)

df6 <- df0
df6 <- df6 %>% arrange(Gene.symbol)

save(df1, df4, df6, file = "helper region.RData")

region.vs.region.index.all <- c(
  "Young Male: Bottom versus Central", "Young Male: Bottom versus Peripheral",
  "Young Male: Bottom versus Top", "Young Male: Central versus Peripheral",
  "Young Male: Central versus Top", "Young Male: Peripheral versus Top",
  "Young Female: Bottom versus Central", "Young Female: Bottom versus Peripheral",
  "Young Female: Bottom versus Top", "Young Female: Central versus Peripheral",
  "Young Female: Central versus Top", "Young Female: Peripheral versus Top",
  "Aged Male: Bottom versus Central", "Aged Male: Bottom versus Peripheral",
  "Aged Male: Bottom versus Top", "Aged Male: Central versus Peripheral",
  "Aged Male: Central versus Top", "Aged Male: Peripheral versus Top",
  "Aged Female: Bottom versus Central", "Aged Female: Bottom versus Peripheral",
  "Aged Female: Bottom versus Top", "Aged Female: Central versus Peripheral",
  "Aged Female: Central versus Top", "Aged Female: Peripheral versus Top"
)


# Prepare the volcano plot R Shiny app helper data--------
# Save Method 4 comparison lists for downstream volcano plots.

rm(list = ls()[!ls() %in% c("age.vs.young4", "female.vs.male4", "region.vs.region4")])
age.vs.young <- age.vs.young4
female.vs.male <- female.vs.male4
region.vs.region <- region.vs.region4
rm(age.vs.young4)
rm(female.vs.male4)
rm(region.vs.region4)

save.image("helper volcano plot.RData")
