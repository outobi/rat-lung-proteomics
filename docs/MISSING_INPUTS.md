# Missing literal input paths

Static check of quoted filenames in read/load/source calls, relative to each script directory. Variable-built paths are not covered. A file may be a generated upstream output or an external dependency rather than permanently unavailable. No filenames were automatically substituted.

| Script | Line | Referenced input |
|---|---:|---|
| [Script_1_protein_number_count/Script_1_protein_number_counts.R](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L14) | 14 | `report.pg_matrix.tsv` |
| [Script_1_protein_number_count/Script_1_protein_number_counts.R](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L23) | 23 | `20220113_uniprot_annotation_MW_IPA_matrisome.RData` |
| [Script_1_protein_number_count/Script_1_protein_number_counts.R](../Script_1_protein_number_count/Script_1_protein_number_counts.R#L46) | 46 | `Matrisome_annotation_output.xlsx` |
| [Script_3_data_preprocessing/Script_3_data_preprocessing.R](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L30) | 30 | `expression_matrix.rds` |
| [Script_3_data_preprocessing/Script_3_data_preprocessing.R](../Script_3_data_preprocessing/Script_3_data_preprocessing.R#L222) | 222 | `T:/Proteomics_Data_Analysis/SOP_Liang/Proteomics data analysis_R template/ggbiplot.r` |
| [Script_4_single_factor_DE/Script_4_Single_factor_DE.R](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L57) | 57 | `Sample_metadata.xlsx` |
| [Script_4_single_factor_DE/Script_4_Single_factor_DE.R](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L377) | 377 | `DE summary1.xlsx` |
| [Script_4_single_factor_DE/Script_4_Single_factor_DE.R](../Script_4_single_factor_DE/Script_4_Single_factor_DE.R#L454) | 454 | `DE summary1.xlsx` |
