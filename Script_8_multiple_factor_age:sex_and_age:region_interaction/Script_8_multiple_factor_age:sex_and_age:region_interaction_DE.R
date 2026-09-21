
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) ## set working dir as where the script locates

library(limma)
library(rio)
library(tidyr)

# DE analysis -------------------------------------------------------------


library(lme4) 
library(lmerTest) 
library(reshape2) 
library(emmeans) 
library(Rmisc) 
library(stringr) 
library(dplyr) 
library(lmerTest) 
library(multcomp) 


# Read the annotated imputed-intensity worksheet for output annotations.
# The first seven columns are assumed to contain protein annotations.

# Load the saved imputed matrices and sample information.
load("UW.Madison.rat.lung.proteomics_data tables_2026-09-15.RData")

data_imp <- read_excel("UW.Madison.rat.lung.proteomics_DIANN_protein intensity_2026-09-15.xlsx",
                       sheet = "data_imp"
)

data7 <- data_table$data_imp

workdata <- read_excel("Sample_metadata.xlsx")

sample_group_number <- read_excel("Real_values_number_in_sample_group.xlsx")
sample_group_number2 <- sample_group_number[sample_group_number$Uniprot %in% rownames(data7), ]

## age * sex interaction---------------------------
# no second round of filtering, all 3600 proteins included, BH correction
# Is the overall age effect, averaged across lung regions, different between females and males?



# factor
f <- info_table$Sample.Group
Region <- info_table$Sample.Region
Age <- info_table$Sample.Age
Sex <- info_table$Sample.Sex
ID <- as.factor(info_table$Sample.PID)



all_interactions_age_sex <- list()

for (i in 1:3600){
  
  workdata[5] <- t(data7[i,])
  
  res <- lmer(
    Protein ~ Age + Sex + Region +
      Age:Region + Sex:Region + Age:Sex +
      (1|Subject),
    data = workdata
  )
  
  # Estimated marginal means for Age x Sex,
  # averaged over the four regions
  emm_age_sex <- emmeans(res, ~ Age * Sex)
  
  # Formal Age x Sex interaction:
  # (Aged - Young in Female) - (Aged - Young in Male)
  interaction_age_sex <- contrast(
    emm_age_sex,
    interaction = "pairwise",
    adjust = "none"
  )
  
  all_interactions_age_sex[[i]] <-
    as.data.frame(summary(interaction_age_sex))
}

age_sex_interaction_df <-
  bind_rows(all_interactions_age_sex)

age_sex_interaction_df$FDR <-
  p.adjust(age_sex_interaction_df$p.value, method = "BH")

rownames(age_sex_interaction_df) <- rownames(data7)

age_sex_sig_summary <- age_sex_interaction_df %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::summarise(
    Female_greater = sum(estimate > 0, na.rm = TRUE),
    Male_greater   = sum(estimate < 0, na.rm = TRUE),
    Total          = dplyr::n()
  )


age_sex_sig_summary
#Female_greater Male_greater Total
#      23           32    55

age_sex_interaction_df <- cbind(file02, age_sex_interaction_df)
export(age_sex_interaction_df,"age sex interaction protein list.xlsx")



# age * sex interaction volcano plot---------------------------------------------------------------------

res1 <- age_sex_interaction_df
res1$LogPValue <- (-1)*log10(res1$FDR)
#res2 <- res1[res1$estimate>1|res1$estimate<(-1),]
res2 <- res1[res1$FDR<0.05,] 
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res2up <- res2[res2$estimate > 0.5849625, ]
res2down <- res2[res2$estimate < (-0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category),]    ## USED ADJUSTED P VALUE FOR CUTOFF

fc_cut <- 1.5 ## set ESTIMATE cut-off


p1 <- ggplot(res1,aes(estimate,LogPValue))
p1 <- p1 + geom_point(shape = 1,color = "gray") #+coord_cartesian(xlim = c(-3,3)) 
p1 <- p1 + geom_point(data = res2,aes(x=estimate,y=LogPValue,colour=estimate),shape = 19,size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour")#+ scale_colour_gradient(high="red",low = "green")
p1 <- p1 + geom_point(data = res3,aes(x=estimate,y=LogPValue),colour = "black", shape = 1,size = 2)
p1 <- p1 + geom_text_repel(data=res3,aes(x=estimate,y=LogPValue,label=Gene.symbol),size=4.5,arrow=arrow(length=unit(0.01,'npc')),point.padding=unit(0.6,'lines'),max.overlaps = 1600)
#p <- p+ geom_vline(xintercept = c(-log2(fc_cut),log2(fc_cut)),colour="red",size=0.5,linetype=2)

p1 <- p1 + geom_hline(yintercept = -log10(0.05),colour="red",size=0.5,linetype=2)
p1 <- p1 + geom_vline(xintercept = 0.5849625,colour="red",size=0.5,linetype=2)
p1 <- p1 + geom_vline(xintercept = -0.5849625,colour="red",size=0.5,linetype=2)

p1 <- p1 + theme_bw(base_size = 16) + labs(x="estimate", y = "-Log10(FDR)") + ggtitle("Age×Sex interaction volcano plot, up in red female, n = 23; up in green male, n = 32")
plot(p1)

ggsave(
  filename = "age sex interaction volcano plot.pdf",
  plot = p1,
  device = "pdf",
  width = 11,
  height = 8.5
)


##age * region interaction---------------

# Is the overall age effect, averaged across sex, different among different regions?


all_interactions_age_region <- list()

for (i in 1:3600){
  
  workdata[5] <- t(data7[i,])
  
  res <- lmer(
    Protein ~ Age + Sex + Region +
      Age:Region + Sex:Region + Age:Sex +
      (1|Subject),
    data = workdata
  )
  
  # Estimated marginal means for Age x Region,
  # averaged over Sex
  emm_age_region <- emmeans(res, ~ Age * Region)
  
  # Formal Age x Region interaction:
  # compare the age effect between pairs of regions
  interaction_age_region <- contrast(
    emm_age_region,
    interaction = "pairwise",
    adjust = "none"
  )
  
  all_interactions_age_region[[i]] <-
    as.data.frame(summary(interaction_age_region))
}

age_region_interaction_df <-
  bind_rows(all_interactions_age_region)


age_region_interaction_df <- age_region_interaction_df %>%
  group_by(Region_pairwise) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()


age_region_sig_summary <- age_region_interaction_df %>%
  dplyr::group_by(Region_pairwise) %>%
  dplyr::summarise(
    n_tested = dplyr::n(),
    n_FDR_005 = sum(FDR < 0.05, na.rm = TRUE),
    Region1_greater = sum(FDR < 0.05 & estimate > 0, na.rm = TRUE),
    Region2_greater = sum(FDR < 0.05 & estimate < 0, na.rm = TRUE)
  )

age_region_sig_summary

#Region_pairwise      n_tested n_FDR_005 Region1_greater Region2_greater
#<fct>                   <int>     <int>           <int>           <int>
#1 Bottom - Central         3600        60              30              30
#2 Bottom - Peripheral      3600        63              28              35
#3 Bottom - Top             3600        72              42              30
#4 Central - Peripheral     3600        75              35              40
#5 Central - Top            3600        68              40              28
#6 Peripheral - Top         3600        84              57              27

library(openxlsx)

# split into a list
age_region_list <- split(
  age_region_interaction_df,
  age_region_interaction_df$Region_pairwise
)

# optional: make sheet names shorter
names(age_region_list) <- gsub(
  " - ",
  "_vs_",
  names(age_region_list)
)

# add file02 if needed
age_region_list <- lapply(age_region_list, function(x){
  cbind(file02, x)
})

# export one Excel file with multiple sheets
openxlsx::write.xlsx(
  age_region_list,
  file = "age region interaction protein list.xlsx",
  rowNames = FALSE
)


# volcano plot-------------
# 1 bottom vs central
res1 <- age_region_list$Bottom_vs_Central
res1$LogPValue <- (-1)*log10(res1$FDR)
res2 <- res1[res1$FDR<0.05,] 
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res2up <- res2[res2$estimate > 0.5849625, ]
res2down <- res2[res2$estimate < (-0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category),]    ## USED ADJUSTED P VALUE FOR CUTOFF

fc_cut <- 1.5 ## set ESTIMATE cut-off


p1 <- ggplot(res1,aes(estimate,LogPValue))
p1 <- p1 + geom_point(shape = 1,color = "gray") #+coord_cartesian(xlim = c(-3,3)) 
p1 <- p1 + geom_point(data = res2,aes(x=estimate,y=LogPValue,colour=estimate),shape = 19,size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour")#+ scale_colour_gradient(high="red",low = "green")
p1 <- p1 + geom_point(data = res3,aes(x=estimate,y=LogPValue),colour = "black", shape = 1,size = 2)
p1 <- p1 + geom_text_repel(data=res3,aes(x=estimate,y=LogPValue,label=Gene.symbol),size=4.5,arrow=arrow(length=unit(0.01,'npc')),point.padding=unit(0.6,'lines'),max.overlaps = 1600)


p1 <- p1 + geom_hline(yintercept = -log10(0.05),colour="red",size=0.5,linetype=2)
p1 <- p1 + geom_vline(xintercept = 0.5849625,colour="red",size=0.5,linetype=2)
p1 <- p1 + geom_vline(xintercept = -0.5849625,colour="red",size=0.5,linetype=2)

p1 <- p1 + theme_bw(base_size = 16) + labs(x="estimate", y = "-Log10(FDR)") + ggtitle("Age×Region interaction volcano plot, up in red bottom, n = 30; up in green central, n = 30")
plot(p1)


ggsave(
  filename = "age region interaction volcano plot bottom vs central.pdf",
  plot = p1,
  device = "pdf",
  width = 11,
  height = 8.5
)



# 2 bottom vs peripheral
res1 <- age_region_list$Bottom_vs_Peripheral
res1$LogPValue <- (-1)*log10(res1$FDR)
res2 <- res1[res1$FDR<0.05,] 
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res2up <- res2[res2$estimate > 0.5849625, ]
res2down <- res2[res2$estimate < (-0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category),]    ## USED ADJUSTED P VALUE FOR CUTOFF

fc_cut <- 1.5 ## set ESTIMATE cut-off


p2 <- ggplot(res1,aes(estimate,LogPValue))
p2 <- p2 + geom_point(shape = 1,color = "gray") #+coord_cartesian(xlim = c(-3,3)) 
p2 <- p2 + geom_point(data = res2,aes(x=estimate,y=LogPValue,colour=estimate),shape = 19,size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour")#+ scale_colour_gradient(high="red",low = "green")
p2 <- p2 + geom_point(data = res3,aes(x=estimate,y=LogPValue),colour = "black", shape = 1,size = 2)
p2 <- p2 + geom_text_repel(data=res3,aes(x=estimate,y=LogPValue,label=Gene.symbol),size=4.5,arrow=arrow(length=unit(0.01,'npc')),point.padding=unit(0.6,'lines'),max.overlaps = 1600)


p2 <- p2 + geom_hline(yintercept = -log10(0.05),colour="red",size=0.5,linetype=2)
p2 <- p2 + geom_vline(xintercept = 0.5849625,colour="red",size=0.5,linetype=2)
p2 <- p2 + geom_vline(xintercept = -0.5849625,colour="red",size=0.5,linetype=2)

p2 <- p2 + theme_bw(base_size = 16) + labs(x="estimate", y = "-Log10(FDR)") + ggtitle("Age×Region interaction volcano plot, up in red bottom, n = 28; up in green peripheral, n = 35")
plot(p2)


ggsave(
  filename = "age region interaction volcano plot bottom vs peripheral.pdf",
  plot = p2,
  device = "pdf",
  width = 11,
  height = 8.5
)



# 3 bottom vs top
res1 <- age_region_list$Bottom_vs_Top
res1$LogPValue <- (-1)*log10(res1$FDR)
res2 <- res1[res1$FDR<0.05,] 
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res2up <- res2[res2$estimate > 0.5849625, ]
res2down <- res2[res2$estimate < (-0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category),]    ## USED ADJUSTED P VALUE FOR CUTOFF

fc_cut <- 1.5 ## set ESTIMATE cut-off


p3 <- ggplot(res1,aes(estimate,LogPValue))
p3 <- p3 + geom_point(shape = 1,color = "gray") #+coord_cartesian(xlim = c(-3,3)) 
p3 <- p3 + geom_point(data = res2,aes(x=estimate,y=LogPValue,colour=estimate),shape = 19,size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour")#+ scale_colour_gradient(high="red",low = "green")
p3 <- p3 + geom_point(data = res3,aes(x=estimate,y=LogPValue),colour = "black", shape = 1,size = 2)
p3 <- p3 + geom_text_repel(data=res3,aes(x=estimate,y=LogPValue,label=Gene.symbol),size=4.5,arrow=arrow(length=unit(0.01,'npc')),point.padding=unit(0.6,'lines'),max.overlaps = 1600)


p3 <- p3 + geom_hline(yintercept = -log10(0.05),colour="red",size=0.5,linetype=2)
p3 <- p3 + geom_vline(xintercept = 0.5849625,colour="red",size=0.5,linetype=2)
p3 <- p3 + geom_vline(xintercept = -0.5849625,colour="red",size=0.5,linetype=2)

p3 <- p3 + theme_bw(base_size = 16) + labs(x="estimate", y = "-Log10(FDR)") + ggtitle("Age×Region interaction volcano plot, up in red bottom, n = 42; up in green top, n = 30")
plot(p3)

ggsave(
  filename = "age region interaction volcano plot bottom vs top.pdf",
  plot = p3,
  device = "pdf",
  width = 11,
  height = 8.5
)



# 4 central vs peripheral
res1 <- age_region_list$Central_vs_Peripheral
res1$LogPValue <- (-1)*log10(res1$FDR)
res2 <- res1[res1$FDR<0.05,] 
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res2up <- res2[res2$estimate > 0.5849625, ]
res2down <- res2[res2$estimate < (-0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category),]    ## USED ADJUSTED P VALUE FOR CUTOFF

fc_cut <- 1.5 ## set ESTIMATE cut-off


p4 <- ggplot(res1,aes(estimate,LogPValue))
p4 <- p4 + geom_point(shape = 1,color = "gray") #+coord_cartesian(xlim = c(-3,3)) 
p4 <- p4 + geom_point(data = res2,aes(x=estimate,y=LogPValue,colour=estimate),shape = 19,size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour")#+ scale_colour_gradient(high="red",low = "green")
p4 <- p4 + geom_point(data = res3,aes(x=estimate,y=LogPValue),colour = "black", shape = 1,size = 2)
p4 <- p4 + geom_text_repel(data=res3,aes(x=estimate,y=LogPValue,label=Gene.symbol),size=4.5,arrow=arrow(length=unit(0.01,'npc')),point.padding=unit(0.6,'lines'),max.overlaps = 1600)

p4 <- p4 + geom_hline(yintercept = -log10(0.05),colour="red",size=0.5,linetype=2)
p4 <- p4 + geom_vline(xintercept = 0.5849625,colour="red",size=0.5,linetype=2)
p4 <- p4 + geom_vline(xintercept = -0.5849625,colour="red",size=0.5,linetype=2)

p4 <- p4 + theme_bw(base_size = 16) + labs(x="estimate", y = "-Log10(FDR)") + ggtitle("Age×Region interaction volcano plot, up in central, n = 35; up in green peripheral, n = 40")
plot(p4)

ggsave(
  filename = "age region interaction volcano plot central vs peripheral.pdf",
  plot = p4,
  device = "pdf",
  width = 11,
  height = 8.5
)




# 5 central vs top
res1 <- age_region_list$Central_vs_Top
res1$LogPValue <- (-1)*log10(res1$FDR)
res2 <- res1[res1$FDR<0.05,] 
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res2up <- res2[res2$estimate > 0.5849625, ]
res2down <- res2[res2$estimate < (-0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category),]    ## USED ADJUSTED P VALUE FOR CUTOFF

fc_cut <- 1.5 ## set ESTIMATE cut-off


p5 <- ggplot(res1,aes(estimate,LogPValue))
p5 <- p5 + geom_point(shape = 1,color = "gray") #+coord_cartesian(xlim = c(-3,3)) 
p5 <- p5 + geom_point(data = res2,aes(x=estimate,y=LogPValue,colour=estimate),shape = 19,size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour")#+ scale_colour_gradient(high="red",low = "green")
p5 <- p5 + geom_point(data = res3,aes(x=estimate,y=LogPValue),colour = "black", shape = 1,size = 2)
p5 <- p5 + geom_text_repel(data=res3,aes(x=estimate,y=LogPValue,label=Gene.symbol),size=4.5,arrow=arrow(length=unit(0.01,'npc')),point.padding=unit(0.6,'lines'),max.overlaps = 1600)

p5 <- p5 + geom_hline(yintercept = -log10(0.05),colour="red",size=0.5,linetype=2)
p5 <- p5 + geom_vline(xintercept = 0.5849625,colour="red",size=0.5,linetype=2)
p5 <- p5 + geom_vline(xintercept = -0.5849625,colour="red",size=0.5,linetype=2)

p5 <- p5 + theme_bw(base_size = 16) + labs(x="estimate", y = "-Log10(FDR)") + ggtitle("Age×Region interaction volcano plot, up in red central, n = 40; up in green top, n = 28")
plot(p5)

ggsave(
  filename = "age region interaction volcano plot central vs top.pdf",
  plot = p5,
  device = "pdf",
  width = 11,
  height = 8.5
)



# 6 peripheral vs top
res1 <- age_region_list$Peripheral_vs_Top
res1$LogPValue <- (-1)*log10(res1$FDR)
res2 <- res1[res1$FDR<0.05,] 
res2 <- res2[abs(res2$estimate) > 0.5849625, ]
res2up <- res2[res2$estimate > 0.5849625, ]
res2down <- res2[res2$estimate < (-0.5849625), ]
res3 <- res2[!is.na(res2$Matrisome.Sub.category),]    ## USED ADJUSTED P VALUE FOR CUTOFF

fc_cut <- 1.5 ## set ESTIMATE cut-off


p6 <- ggplot(res1,aes(estimate,LogPValue))
p6 <- p6 + geom_point(shape = 1,color = "gray") #+coord_cartesian(xlim = c(-3,3)) 
p6 <- p6 + geom_point(data = res2,aes(x=estimate,y=LogPValue,colour=estimate),shape = 19,size = 2) + scale_colour_gradient2(low = "green", mid = "lightyellow", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour")#+ scale_colour_gradient(high="red",low = "green")
p6 <- p6 + geom_point(data = res3,aes(x=estimate,y=LogPValue),colour = "black", shape = 1,size = 2)
p6 <- p6 + geom_text_repel(data=res3,aes(x=estimate,y=LogPValue,label=Gene.symbol),size=4.5,arrow=arrow(length=unit(0.01,'npc')),point.padding=unit(0.6,'lines'),max.overlaps = 1600)

p6 <- p6 + geom_hline(yintercept = -log10(0.05),colour="red",size=0.5,linetype=2)
p6 <- p6 + geom_vline(xintercept = 0.5849625,colour="red",size=0.5,linetype=2)
p6 <- p6 + geom_vline(xintercept = -0.5849625,colour="red",size=0.5,linetype=2)

p6 <- p6 + theme_bw(base_size = 16) + labs(x="estimate", y = "-Log10(FDR)") + ggtitle("Age×Region interaction volcano plot, up in red peripheral, n = 57; up in green top, n = 27")
plot(p6)

ggsave(
  filename = "age region interaction volcano plot peripheral vs top.pdf",
  plot = p6,
  device = "pdf",
  width = 11,
  height = 8.5
)

library(patchwork)

p_summary <-
  (p1 | p2 | p3) /
  (p4 | p5 | p6)

p_summary


ggsave(
  filename = "Age_region_interaction_volcano_summary_2x3.pdf",
  plot = p_summary,
  device = pdf,
  width = 39,
  height = 26
)


