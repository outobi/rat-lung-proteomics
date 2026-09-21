# File I/O source index

This index records read/write calls and nearby source lines. It is a navigation aid, not proof that all input files exist or that every output can be regenerated. Line numbers refer to the archived scripts. Dynamic filenames and commented calls should be interpreted in the source.

## Script_10_R_Shiny/Script-10_R_Shiny.R

- [Line 25](../Script_10_R_Shiny/Script-10_R_Shiny.R#L25): `load("Script_6_multiple_factor_DE.RData")  `
- [Line 88](../Script_10_R_Shiny/Script-10_R_Shiny.R#L88): `save(df1, df4, df6, file = "helper age.RData")  `
- [Line 157](../Script_10_R_Shiny/Script-10_R_Shiny.R#L157): `save(df1, df4, df6, file = "helper sex.RData")  `
- [Line 235](../Script_10_R_Shiny/Script-10_R_Shiny.R#L235): `save(df1, df4, df6, file = "helper region.RData")  region.vs.region.index.all <- c(`
- [Line 264](../Script_10_R_Shiny/Script-10_R_Shiny.R#L264): `save.image("helper volcano plot.RData")`

## Script_10_R_Shiny/apps/protein_explorer/app.R

- [Line 9](../Script_10_R_Shiny/apps/protein_explorer/app.R#L9): `load("helper age.RData")    # contains df4_age (previously df4) and df1 df4_age <- df4 load("helper region.RData") # contains df4_region`
- [Line 11](../Script_10_R_Shiny/apps/protein_explorer/app.R#L11): `load("helper region.RData") # contains df4_region df4_region <- df4 load("helper sex.RData")    # contains df4_sex`
- [Line 13](../Script_10_R_Shiny/apps/protein_explorer/app.R#L13): `load("helper sex.RData")    # contains df4_sex df4_sex <- df4 df4_age$A1bg`
- [Line 349](../Script_10_R_Shiny/apps/protein_explorer/app.R#L349): `write.csv(all_data, file, row.names = FALSE) } )`

## Script_10_R_Shiny/apps/volcano_explorer/app.R

- [Line 8](../Script_10_R_Shiny/apps/volcano_explorer/app.R#L8): `load("helper volcano plot.RData") # Define comparison groups `
- [Line 330](../Script_10_R_Shiny/apps/volcano_explorer/app.R#L330): `ggsave(file, plot, width = 10, height = 8, dpi = 300) } )`
- [Line 339](../Script_10_R_Shiny/apps/volcano_explorer/app.R#L339): `write.csv(get_comparison_data(), file, row.names = FALSE) } )`

## Script_1_protein_number_count/Script_1_protein_number_counts.R

- [Line 14](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L14): `data1 <- import("report.pg_matrix.tsv") rownames(data1) <- gsub(";.*","",data1$Protein.Group) ## regex data1 <-  data1[,-(1:4)]`
- [Line 23](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L23): `load("20220113_uniprot_annotation_MW_IPA_matrisome.RData") anno_table <- anno_mapping[[match("rat",names(anno_mapping))]] anno_table <- anno_table[match(rownames(data1),anno_table$Uniprot.ID),]`
- [Line 32](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L32): `write_xlsx( anno_table, "Matrisome_annotation_input.xlsx"`
- [Line 46](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L46): `table1 <- read_excel("Matrisome_annotation_output.xlsx") all(table1$Uniprot.ID == rownames(data1)) # True table1 <- cbind(table1, data1)`
- [Line 106](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L106): `sample.table <- read_excel("Sample_index.xlsx")  `
- [Line 188](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L188): `write_xlsx(tables_list, "Protein_counts_summary.xlsx")  `
- [Line 192](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L192): `saveRDS(data1, file =  "expression_matrix.rds") `

## Script_2_peptide_number_count/Script_2_peptide_number_count.R

- [Line 17](../Script_2_peptide_number_count/Script_2_peptide_number_count.R#L17): `peptide <- read_excel("Peptide.xlsx")  `
- [Line 30](../Script_2_peptide_number_count/Script_2_peptide_number_count.R#L30): `raw_table <- read_excel("Protein_counts_summary.xlsx", sheet = "Raw_table_annotation") colnames(raw_table)[7] <- "Matrisome.Sub.Category"`
- [Line 111](../Script_2_peptide_number_count/Script_2_peptide_number_count.R#L111): `sample.table <- read_excel("Sample_metadata.xlsx")  `
- [Line 176](../Script_2_peptide_number_count/Script_2_peptide_number_count.R#L176): `export(tables_list, "Peptide_counts_summary.xlsx")`

## Script_3_data_preprocessing/Script_3_data_preprocessing.R

- [Line 30](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L30): `data1 <- readRDS("expression_matrix.rds")  # Read the seven annotation columns from the protein-count summary workbook.`
- [Line 33](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L33): `matrisome_anno_table <- read_excel("Protein_counts_summary.xlsx", sheet = "Raw_table_annotation") matrisome_anno_table <- matrisome_anno_table[c(1:7)]`
- [Line 39](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L39): `info_table <- import(sample_info)  colnames(info_table) <- c("Sample.ID", "Sample.Group", "Sample.PID", "Sample.Animal.Number", "Sample.Age", "Sample.Region", "Sample.Sex")`
- [Line 95](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L95): `export(sample_group_number, "Real_values_number_in_sample_group.xlsx")  sample_group_number <- sample_group_number[-1]`
- [Line 146](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L146): `pdf(paste0(study_number, "_QC plots ", Sys.Date(), ".pdf"))  hist(na.count1, freq = F, xlab = "NA count", ylab = "% protein id", main = "Missing value/protein before filtering", col = 4, xlim = c(0, ncol(data1)))`
- [Line 222](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L222): `source("T:/Proteomics_Data_Analysis/SOP_Liang/Proteomics data analysis_R template/ggbiplot.r")  `
- [Line 254](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L254): `pdf("PCA.pdf") plot(g2) plot(g3)`
- [Line 270](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L270): `save(data_table, matrisome_anno_table, info_table, file = paste0(study_number, "_data tables_", Sys.Date(), ".RData"))  `
- [Line 294](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L294): `export(data_table_export, file = paste0(study_number, "_", search_method, "_protein intensity_", Sys.Date(), ".xlsx")) cat(paste0(Sys.time(), " Excel file of protein intensity data tables is created", "\n"), file = log_file, append = TRUE) `
- [Line 299](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L299): `save.image("Script_3_data_preprocessing_data_summary.RData")`

## Script_4_single_factor_DE/Script_4_Single_factor_DE.R

- [Line 31](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L31): `load("UW.Madison.rat.lung.proteomics_data tables_2026-09-15.RData")  `
- [Line 36](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L36): `data_imp <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein intensity_2026-09-15.xlsx", sheet = "data_imp" )`
- [Line 57](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L57): `workdata <- read_excel("Sample_metadata.xlsx")  `
- [Line 221](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L221): `export(single_factor_DEP_summary, file = "Single_factor_DEP_summary.xlsx")  age_DE <- read_excel("Single_factor_DEP_summary.xlsx", sheet = "age")`
- [Line 223](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L223): `age_DE <- read_excel("Single_factor_DEP_summary.xlsx", sheet = "age")  `
- [Line 253](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L253): `sex_DE <- read_excel("Single_factor_DEP_summary.xlsx", sheet = "sex") res1 <- sex_DE res1$LogPValue <- (-1) * log10(res1$p.value.adjusted)`
- [Line 279](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L279): `pdf("20260915_volcano_plots_all_comparisons.pdf", width = 12, height = 8)  plot(p1)`
- [Line 297](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L297): `res1 <- read_excel("Single_factor_DEP_summary.xlsx", sheet = region_comparisons[[comparison_name]]) res1$LogPValue <- (-1) * log10(res1$p.value.adjusted) res2 <- res1[res1$p.value.adjusted < 0.05, ]`
- [Line 377](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L377): `data <- read_excel("DE summary1.xlsx", sheet = region_comparisons[[comparison_name]]) counts <- count_DEPs(data) dep_counts <- rbind(`
- [Line 454](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L454): `data <- read_excel("DE summary1.xlsx", sheet = region_comparisons[[comparison_name]]) counts <- count_DEPs_FC(data, fc_cutoff) dep_counts_FC <- rbind(`
- [Line 516](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L516): `pdf("20260915_DEP_numbers_barplot.pdf", width = 10, height = 6) print(p_dep) dev.off()`
- [Line 564](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L564): `pdf("20260915_matrisome_DEP_numbers_barplot.pdf", width = 10, height = 6) print(p_dep_matrisome) dev.off()`
- [Line 586](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L586): `sample_group <- read_excel("Sample_index.xlsx") sample_group_order <- unlist(sample_group$Index) `
- [Line 641](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L641): `pdf("20260915_matrisome_DEP_heatmap_and_barplot.pdf", width = 10, height = 6) print(p8) print(p9)`
- [Line 648](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L648): `save.image("Single_factor_DEP_summary_data_summary.RData")`

## Script_5_single_factor_pathway/age_vs_young_GSEA/GSEA_aged_vs_young_single_factor.R

- [Line 19](../Script_5_single_factor_pathway/age_vs_young_GSEA/GSEA_aged_vs_young_single_factor.R#L19): `up <- read_excel( "GSEA_summary_aged_vs_young_single_factor.xlsx", sheet = "GO aged"`
- [Line 68](../Script_5_single_factor_pathway/age_vs_young_GSEA/GSEA_aged_vs_young_single_factor.R#L68): `down <- read_excel( "GSEA_summary_aged_vs_young_single_factor.xlsx", sheet = "GO young"`
- [Line 119](../Script_5_single_factor_pathway/age_vs_young_GSEA/GSEA_aged_vs_young_single_factor.R#L119): `up <- read_excel( "GSEA_summary_aged_vs_young_single_factor.xlsx", sheet = "Reactome young"`
- [Line 162](../Script_5_single_factor_pathway/age_vs_young_GSEA/GSEA_aged_vs_young_single_factor.R#L162): `pdf("GSEA_age_vs_young_single_factor.pdf", width = 11, height = 8.5) print(p1) print(p2)`

## Script_5_single_factor_pathway/age_vs_young_ORA/ORA_aged_vs_young_single_factor.R

- [Line 16](../Script_5_single_factor_pathway/age_vs_young_ORA/ORA_aged_vs_young_single_factor.R#L16): `up <- read_excel("ORA_summary.xlsx", sheet = "aged GO CC")  # Order Description based on the absolute value of Log10(FDR) in descending order`
- [Line 37](../Script_5_single_factor_pathway/age_vs_young_ORA/ORA_aged_vs_young_single_factor.R#L37): `up <- read_excel("ORA_summary.xlsx", sheet = "aged GO summary")  up <- up %>%`
- [Line 66](../Script_5_single_factor_pathway/age_vs_young_ORA/ORA_aged_vs_young_single_factor.R#L66): `up <- read_excel("ORA_summary.xlsx", sheet = "young Reactome")  # Order Description based on the absolute value of Log10(FDR) in descending order`
- [Line 89](../Script_5_single_factor_pathway/age_vs_young_ORA/ORA_aged_vs_young_single_factor.R#L89): `up <- read_excel("ORA_summary.xlsx", sheet = "young GO summary")  `
- [Line 116](../Script_5_single_factor_pathway/age_vs_young_ORA/ORA_aged_vs_young_single_factor.R#L116): `pdf("ORA_age_vs_young_single_factor.pdf", width = 11, height = 8.5) print(p1) print(p2)`

## Script_5_single_factor_pathway/female_vs_male_GSEA/GSEA_summary_female_vs_male_single_factor.R

- [Line 21](../Script_5_single_factor_pathway/female_vs_male_GSEA/GSEA_summary_female_vs_male_single_factor.R#L21): `down <- read_excel("GSEA_summary_female_vs_male_single_factor.xlsx")  `
- [Line 63](../Script_5_single_factor_pathway/female_vs_male_GSEA/GSEA_summary_female_vs_male_single_factor.R#L63): `pdf("GSEA_female_vs_male_single_factor.pdf", width = 11, height = 8.5) print(p4) dev.off()`

## Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R

- [Line 19](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L19): `up <- read_excel("GSEA_summary_central_vs_peripheral.xlsx", sheet = "GO" )`
- [Line 50](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L50): `up <- read_excel("GSEA_summary_central_vs_peripheral.xlsx", sheet = "Reactome central" )`
- [Line 76](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L76): `down <- read_excel("GSEA_summary_central_vs_peripheral.xlsx", sheet = "Reactome peripheral" )`
- [Line 102](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L102): `down <- read_excel("GSEA_summary_bottom_vs_top.xlsx", sheet = "GO" )`
- [Line 134](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L134): `down <- read_excel("GSEA_summary_bottom_vs_top.xlsx", sheet = "Reactome" )`
- [Line 164](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L164): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO bottom vs central")  # Sort terms and preserve the resulting y-axis order (bottom to top).`
- [Line 194](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L194): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "Reactome bottom vs central")  # Sort terms and preserve the resulting y-axis order (bottom to top).`
- [Line 218](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L218): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO bottom vs peri")  # Sort terms and preserve the resulting y-axis order (bottom to top).`
- [Line 248](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L248): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO central vs top")  # Sort terms and preserve the resulting y-axis order (bottom to top).`
- [Line 278](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L278): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "Reactome central vs top")  # Sort terms and preserve the resulting y-axis order (bottom to top).`
- [Line 302](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L302): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "Reactome top vs central" )`
- [Line 328](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L328): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO peri vs top" )`
- [Line 360](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L360): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "GO top vs peri" )`
- [Line 392](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L392): `up <- read_excel("Supplementary_table_region_vs_reigon_single_factor_GSEA.xlsx", sheet = "Reactome peri vs top" )`
- [Line 421](../Script_5_single_factor_pathway/region_vs_region_GSEA/GSEA_region_vs_region_single_factor.R#L421): `pdf("GSEA_region_vs_region_single_factor.pdf", width = 14, height = 8.5) for (region_plot in region_plots) { print(region_plot)`

## Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R

- [Line 30](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L30): `load("UW.Madison.rat.lung.proteomics_data tables_2026-09-15.RData")  `
- [Line 35](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L35): `data_imp <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein intensity_2026-09-15.xlsx", sheet = "data_imp" )`
- [Line 41](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L41): `workdata <- read_excel("Sample_metadata.xlsx")  sample_group_number <- read_excel("Real_values_number_in_sample_group.xlsx")`
- [Line 43](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L43): `sample_group_number <- read_excel("Real_values_number_in_sample_group.xlsx") sample_group_number2 <- sample_group_number[sample_group_number$Uniprot %in% rownames(data7), ] `
- [Line 233](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L233): `export(age.vs.young4, "age.vs.young4.xlsx") export(female.vs.male4, "female.vs.male4.xlsx") export(region.vs.region4, "region.vs.region4.xlsx")`
- [Line 234](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L234): `export(female.vs.male4, "female.vs.male4.xlsx") export(region.vs.region4, "region.vs.region4.xlsx") `
- [Line 235](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L235): `export(region.vs.region4, "region.vs.region4.xlsx")  `
- [Line 297](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L297): `export(results_list4, file = "Method 4 DE_protein_summary_all_comparisons_FDR0.05_FC1.5.xlsx")  `
- [Line 446](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L446): `pdf("DE_protein_summary_plots.pdf", width = 12, height = 18) grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2) dev.off()`
- [Line 593](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L593): `pdf("DE_matrisome_protein_summary_plots.pdf", width = 12, height = 18) grid.arrange(p7, p8, p9, p10, p11, p12, ncol = 2) dev.off()`
- [Line 623](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L623): `pdf("age_vs_young_DEP_multiple_facotor_volcano_plot.pdf", width = 12, height = 8)  for (i in 1:8) {`
- [Line 733](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L733): `export(summary_table3, "summary_table3.xlsx")  # Plot the 38 selected matrisome DEPs across eight comparisons.`
- [Line 783](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L783): `boxplot_sample_order <- read.table("Sample_group_order.txt", header = TRUE) boxplot_intensities <- boxplot_intensities[, boxplot_sample_order$Sample.ID] boxplot_groups <- c(`
- [Line 793](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L793): `pdf("All_Proteins_Boxplots.pdf", width = 11.5, height = 8)  for (i in seq_len(nrow(boxplot_intensities))) {`
- [Line 858](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L858): `ggsave("Age_vs_young_DEP_correlation_map.pdf", pp3, width = 8, height = 7)  save.image("Script_6_multiple_factor_DE.RData")`
- [Line 860](../Script_6_multiple_factor_DE/Script_6_multiple_factor_DE.R#L860): `save.image("Script_6_multiple_factor_DE.RData")`

## Script_7_mulltiple_factor_pathway/Figure S8 TOP YOUNG, female vs male ORA/bubble plot.R

- [Line 20](../Script_7_mulltiple_factor_pathway/Figure%20S8%20TOP%20YOUNG%2C%20female%20vs%20male%20ORA/bubble%20plot.R#L20): `up <- read_excel("ORA list in female vs male top young.xlsx", sheet = "GO")  # Order terms within each GO category by increasing Log10(FDR).`
- [Line 51](../Script_7_mulltiple_factor_pathway/Figure%20S8%20TOP%20YOUNG%2C%20female%20vs%20male%20ORA/bubble%20plot.R#L51): `up <- read_excel("ORA list in female vs male top young.xlsx", sheet = "Reactome")  # Order terms by increasing absolute Log10(FDR); larger values appear at the top.`
- [Line 74](../Script_7_mulltiple_factor_pathway/Figure%20S8%20TOP%20YOUNG%2C%20female%20vs%20male%20ORA/bubble%20plot.R#L74): `pdf("top young female > male, represatitive ORA and Reactome.pdf", width = 11, height = 8.5)  print(p1)`

## Script_7_mulltiple_factor_pathway/age vs young/S7_multiple_factor_pathway.R

- [Line 14](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L14): `data0 <- read_excel("summary_table3.xlsx")  data0 <- as.data.frame(data0)`
- [Line 56](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L56): `pdf("Figure_S11a_UpSet_up_in_age.pdf", width = 10, height = 7) print(upset_plot_aged, newpage = FALSE) `
- [Line 96](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L96): `pdf("Figure_S11b_UpSet_up_in_young.pdf", width = 10, height = 7) print(upset_plot_young, newpage = FALSE) `
- [Line 124](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L124): `export(test1, "female_up.xlsx") # 1. Aged > young in females, pooled across regions. export(test2, "male_up.xlsx") # 2. Aged > young in males, pooled across regions. # make venn graph in S11a manually`
- [Line 125](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L125): `export(test2, "male_up.xlsx") # 2. Aged > young in males, pooled across regions. # make venn graph in S11a manually `
- [Line 137](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L137): `export(test1, "bottom_up.xlsx") # 3. Aged > young in bottom, pooled across sexes. export(test2, "central_up.xlsx") # 4. Aged > young in central, pooled across sexes. export(test3, "peripheral_up.xlsx") # 5. Aged > young in peripheral, pooled across sexes.`
- [Line 138](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L138): `export(test2, "central_up.xlsx") # 4. Aged > young in central, pooled across sexes. export(test3, "peripheral_up.xlsx") # 5. Aged > young in peripheral, pooled across sexes. export(test4, "top_up.xlsx") # 6. Aged > young in top, pooled across sexes.`
- [Line 139](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L139): `export(test3, "peripheral_up.xlsx") # 5. Aged > young in peripheral, pooled across sexes. export(test4, "top_up.xlsx") # 6. Aged > young in top, pooled across sexes. # make venn graph in S11b manually`
- [Line 140](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L140): `export(test4, "top_up.xlsx") # 6. Aged > young in top, pooled across sexes. # make venn graph in S11b manually `
- [Line 165](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L165): `export(test3, "female_down.xlsx") # 7. Young > aged in females, pooled across regions. export(test4, "male_down.xlsx") # 8. Young > aged in males, pooled across regions. # make venn graph in S11a manually`
- [Line 166](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L166): `export(test4, "male_down.xlsx") # 8. Young > aged in males, pooled across regions. # make venn graph in S11a manually `
- [Line 178](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L178): `export(test1, "bottom_down.xlsx") # 9. Young > aged in bottom, pooled across sexes. export(test2, "central_down.xlsx") # 10. Young > aged in central, pooled across sexes. export(test3, "peripheral_down.xlsx") # 11. Young > aged in peripheral, pooled across sexes.`
- [Line 179](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L179): `export(test2, "central_down.xlsx") # 10. Young > aged in central, pooled across sexes. export(test3, "peripheral_down.xlsx") # 11. Young > aged in peripheral, pooled across sexes. export(test4, "top_down.xlsx") # 12. Young > aged in top, pooled across sexes.`
- [Line 180](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L180): `export(test3, "peripheral_down.xlsx") # 11. Young > aged in peripheral, pooled across sexes. export(test4, "top_down.xlsx") # 12. Young > aged in top, pooled across sexes. # make venn graph in S11b manually`
- [Line 181](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L181): `export(test4, "top_down.xlsx") # 12. Young > aged in top, pooled across sexes. # make venn graph in S11b manually `
- [Line 187](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L187): `BP_female_up <- read_excel( "ORA summary on sex.xlsx", sheet = "female aged"`
- [Line 193](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L193): `BP_male_up <- read_excel( "ORA summary on sex.xlsx", sheet = "male aged"`
- [Line 231](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L231): `BP_female_down <- read_excel( "ORA summary on sex.xlsx", sheet = "female young"`
- [Line 237](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L237): `BP_male_down <- read_excel( "ORA summary on sex.xlsx", sheet = "male young"`
- [Line 275](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L275): `pdf("Age_vs_young_in_female_or_male_ORA_summary.pdf", width = 10, height = 7)  print(p1)`
- [Line 289](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L289): `BP_bottom_up <- read_excel( "ORA summary on region.xlsx", sheet = "bottom up in aged"`
- [Line 295](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L295): `BP_central_up <- read_excel( "ORA summary on region.xlsx", sheet = "central up in aged"`
- [Line 301](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L301): `BP_peripheral_up <- read_excel( "ORA summary on region.xlsx", sheet = "peri up in aged"`
- [Line 307](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L307): `BP_top_up <- read_excel( "ORA summary on region.xlsx", sheet = "top up in aged"`
- [Line 376](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L376): `pdf("Age_>_young_in_regions_ORA_summary.pdf", width = 10, height = 7)  print(p11)`
- [Line 387](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L387): `BP_bottom_down <- read_excel( "ORA summary on region.xlsx", sheet = "bottom up in young"`
- [Line 393](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L393): `BP_central_down <- read_excel( "ORA summary on region.xlsx", sheet = "central up in young"`
- [Line 399](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L399): `BP_peripheral_down <- read_excel( "ORA summary on region.xlsx", sheet = "peri up in young"`
- [Line 405](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L405): `BP_top_down <- read_excel( "ORA summary on region.xlsx", sheet = "top up in young"`
- [Line 474](../Script_7_mulltiple_factor_pathway/age%20vs%20young/S7_multiple_factor_pathway.R#L474): `pdf("Age_<_young_in_regions_ORA_summary.pdf", width = 10, height = 7)  print(p55)`

## Script_7_mulltiple_factor_pathway/female vs male/upset and venngraph.R

- [Line 13](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L13): `female_group <- read_excel( "Female vs male DEP in 8 pairwise groups.xlsx", sheet = "8 groups up in female"`
- [Line 76](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L76): `pdf("Figure_S10a_UpSet_up_in_female.pdf", width = 10, height = 7) print(upset_plot_female, newpage = FALSE) dev.off()`
- [Line 124](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L124): `sink(output_file) cat("PROTEIN LISTS FOR EACH INTERSECTION GROUP IN UPSET PLOT (FEMALE DATA)\n") cat(paste(rep("=", 80), collapse = ""), "\n\n")`
- [Line 139](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L139): `sink()  # Export the complete protein list for each individual group.`
- [Line 143](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L143): `sink(output_file2) categories <- colnames(protein_df) for (i in seq_along(categories)) {`
- [Line 153](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L153): `sink()  cat("Intersection protein lists saved to", output_file, "\n")`
- [Line 160](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L160): `male_group <- read_excel( "Female vs male DEP in 8 pairwise groups.xlsx", sheet = "8 groups up in male"`
- [Line 223](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L223): `pdf("Figure_S10b_UpSet_up_in_male.pdf", width = 10, height = 7) print(upset_plot_male, newpage = FALSE) dev.off()`
- [Line 271](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L271): `sink(output_file) cat("PROTEIN LISTS FOR EACH INTERSECTION GROUP IN UPSET PLOT\n") cat(paste(rep("=", 80), collapse = ""), "\n\n")`
- [Line 286](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L286): `sink()  # Export the complete protein list for each individual group.`
- [Line 290](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L290): `sink(output_file2) categories <- colnames(protein_df) for (i in seq_along(categories)) {`
- [Line 300](../Script_7_mulltiple_factor_pathway/female%20vs%20male/upset%20and%20venngraph.R#L300): `sink()  cat("Intersection protein lists saved to", output_file, "\n")`

## Script_8_multiple_factor_age:sex_and_age:region_interaction/GSEA/Script_8_multiple_factor_age:sex_and_age:region_interaction_GSEA.R

- [Line 13](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L13): `up <- read_excel("gsea summary.xlsx", sheet = "female vs male go")  up <- up %>%`
- [Line 36](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L36): `ggsave( filename = "GSEA GO, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group.pdf", plot = p1,`
- [Line 46](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L46): `up <- read_excel("gsea summary.xlsx", sheet = "female vs male reactome")  `
- [Line 66](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L66): `ggsave( filename = "GSEA Reactome, (Aged vs Young)Female - (Aged vs Young)Male, enriched in the female group.pdf", plot = p2,`
- [Line 79](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L79): `up <- read_excel("gsea summary.xlsx", sheet = "bottom vs peri GO")  up <- up %>%`
- [Line 102](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L102): `ggsave( filename = "GSEA GO MF, (Aged vs Young)Bottom - (Aged vs Young)Peripheral, enriched in the bottom group.pdf", plot = p3,`
- [Line 114](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L114): `up <- read_excel("gsea summary.xlsx", sheet = "central vs peri GO")  up <- up %>%`
- [Line 137](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L137): `ggsave( filename = "GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Peripheral, enriched in the central group.pdf", plot = p4,`
- [Line 148](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L148): `up <- read_excel("gsea summary.xlsx", sheet = "central vs peri reactome")  `
- [Line 168](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L168): `ggsave( filename = "GSEA GO Reactome, (Aged vs Young)Central - (Aged vs Young)Peripheral, enriched in the central group.pdf", plot = p5,`
- [Line 181](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L181): `up <- read_excel("gsea summary.xlsx", sheet = "central vs top GO")  up <- up %>%`
- [Line 204](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L204): `ggsave( filename = "GSEA GO CC, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group.pdf", plot = p6,`
- [Line 215](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L215): `up <- read_excel("gsea summary.xlsx", sheet = "central vs top reactome")  `
- [Line 234](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L234): `ggsave( filename = "GSEA GO Reactome, (Aged vs Young)Central - (Aged vs Young)Top, enriched in the central group.pdf", plot = p7,`
- [Line 243](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/GSEA/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_GSEA.R#L243): `save.image("GSEA plotting summary 20260909.RData")`

## Script_8_multiple_factor_age:sex_and_age:region_interaction/Script_8_multiple_factor_age:sex_and_age:region_interaction_DE.R

- [Line 26](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L26): `load("UW.Madison.rat.lung.proteomics_data tables_2026-09-15.RData")  data_imp <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein intensity_2026-09-15.xlsx",`
- [Line 28](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L28): `data_imp <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein intensity_2026-09-15.xlsx", sheet = "data_imp" )`
- [Line 34](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L34): `workdata <- read_excel("Sample_metadata.xlsx")  sample_group_number <- read_excel("Real_values_number_in_sample_group.xlsx")`
- [Line 36](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L36): `sample_group_number <- read_excel("Real_values_number_in_sample_group.xlsx") sample_group_number2 <- sample_group_number[sample_group_number$Uniprot %in% rownames(data7), ] `
- [Line 105](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L105): `export(age_sex_interaction_df,"age sex interaction protein list.xlsx")  `
- [Line 137](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L137): `ggsave( filename = "age sex interaction volcano plot.pdf", plot = p1,`
- [Line 231](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L231): `openxlsx::write.xlsx( age_region_list, file = "age region interaction protein list.xlsx",`
- [Line 266](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L266): `ggsave( filename = "age region interaction volcano plot bottom vs central.pdf", plot = p1,`
- [Line 303](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L303): `ggsave( filename = "age region interaction volcano plot bottom vs peripheral.pdf", plot = p2,`
- [Line 339](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L339): `ggsave( filename = "age region interaction volcano plot bottom vs top.pdf", plot = p3,`
- [Line 374](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L374): `ggsave( filename = "age region interaction volcano plot central vs peripheral.pdf", plot = p4,`
- [Line 410](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L410): `ggsave( filename = "age region interaction volcano plot central vs top.pdf", plot = p5,`
- [Line 445](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L445): `ggsave( filename = "age region interaction volcano plot peripheral vs top.pdf", plot = p6,`
- [Line 462](../Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction/Script_8_multiple_factor_age%3Asex_and_age%3Aregion_interaction_DE.R#L462): `ggsave( filename = "Age_region_interaction_volcano_summary_2x3.pdf", plot = p_summary,`

## Script_9_region_cell_type_enrichment_with_scRNAseq/region protein marker in scRNAseq cell type/region_protein_marker_in_scRNAseq _cell_type.R

- [Line 9](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L9): `m_counts <- readMM("GSM3926541_mRat.raw.data.mtx.gz")  m_genes <- read.delim(`
- [Line 11](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L11): `m_genes <- read.delim( "GSM3926541_mRat.genes.tsv.gz", header = TRUE,`
- [Line 19](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L19): `m_barcodes <- read.delim( "GSM3926541_mRat.barcodes.tsv.gz", header = TRUE,`
- [Line 26](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L26): `f_counts <- readMM("GSM3926542_fRat.raw.data.mtx.gz")  f_genes <- read.delim(`
- [Line 28](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L28): `f_genes <- read.delim( "GSM3926542_fRat.genes.tsv.gz", header = TRUE,`
- [Line 35](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L35): `f_barcodes <- read.delim( "GSM3926542_fRat.barcodes.tsv.gz", header = TRUE,`
- [Line 225](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L225): `write.csv( markers, "Raredon_rat_FindAllMarkers.csv",`
- [Line 231](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L231): `write.csv( top10, "Raredon_rat_top10_markers_per_cluster.csv",`
- [Line 325](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L325): `meta <- read_excel("sample metadata.xlsx") expr <- read_excel("expression matrix.xlsx") expr<- as.data.frame(expr)`
- [Line 326](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L326): `expr <- read_excel("expression matrix.xlsx") expr<- as.data.frame(expr) rownames(expr) <- expr$Uniprot.ID`
- [Line 460](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L460): `export(region_specific_feature, "region_specific_protein_summary.xlsx")  `
- [Line 592](../Script_9_region_cell_type_enrichment_with_scRNAseq/region%20protein%20marker%20in%20scRNAseq%20cell%20type/region_protein_marker_in_scRNAseq%20_cell_type.R#L592): `export(region_specific_proteins, "region_specific_upregulated_proteins.xlsx")  `

## Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq cell type marker score in regions/Script_9_scRNA_cell_type_marker_in_regions.R

- [Line 11](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L11): `meta <- read_excel("sample metadata.xlsx") expr <- read_excel("expression matrix.xlsx") expr<- as.data.frame(expr)`
- [Line 12](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L12): `expr <- read_excel("expression matrix.xlsx") expr<- as.data.frame(expr) rownames(expr) <- expr$Uniprot.ID`
- [Line 17](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L17): `m_counts <- readMM("GSM3926541_mRat.raw.data.mtx.gz")  m_genes <- read.delim(`
- [Line 19](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L19): `m_genes <- read.delim( "GSM3926541_mRat.genes.tsv.gz", header = TRUE,`
- [Line 27](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L27): `m_barcodes <- read.delim( "GSM3926541_mRat.barcodes.tsv.gz", header = TRUE,`
- [Line 34](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L34): `f_counts <- readMM("GSM3926542_fRat.raw.data.mtx.gz")  f_genes <- read.delim(`
- [Line 36](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L36): `f_genes <- read.delim( "GSM3926542_fRat.genes.tsv.gz", header = TRUE,`
- [Line 43](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L43): `f_barcodes <- read.delim( "GSM3926542_fRat.barcodes.tsv.gz", header = TRUE,`
- [Line 234](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L234): `write.csv( markers, "Raredon_rat_FindAllMarkers.csv",`
- [Line 240](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L240): `write.csv( top10, "Raredon_rat_top10_markers_per_cluster.csv",`
- [Line 386](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L386): `write.csv( celltype_markers, "Raredon_rat_celltype_FindAllMarkers.csv",`
- [Line 436](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L436): `write.csv( celltype_marker_top30, "Raredon_rat_celltype_top30_proteomics_markers.csv",`
- [Line 866](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L866): `write.csv( emm_region, "celltype_region_emmeans.csv",`
- [Line 872](../Script_9_region_cell_type_enrichment_with_scRNAseq/scRNAseq%20cell%20type%20marker%20score%20in%20regions/Script_9_scRNA_cell_type_marker_in_regions.R#L872): `write.csv( celltype_region_omnibus, "celltype_region_LMM_FDR.csv",`

