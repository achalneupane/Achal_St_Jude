####################################################
## Read Phenotype data from 2_genetic_data_P_LP.r ##
####################################################
PHENO.ANY_SN <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/2_genetic_data_P_LP.rds")
# PHENO.ANY_SN <- readRDS("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/2_genetic_data_P_LP.rds")


#########################################################
## Read PRS score files and merge them to PHENO.ANY_SN ##
#########################################################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/prs_out")
# setwd("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/prs_out")

## Meningioma_Claus_prs.profile
meningioma_Claus <- read.table("Meningioma_Claus_prs.profile", header = T)
PHENO.ANY_SN$Meningiom_Claus_PRS <-  meningioma_Claus$SCORE [match(PHENO.ANY_SN$sjlid, meningioma_Claus$IID)]

## Meningioma_Dobbins_prs.profile
Meningioma_Dobbins <- read.table("Meningioma_Dobbins_prs.profile", header = T)
PHENO.ANY_SN$Meningioma_Dobbins_PRS <-  Meningioma_Dobbins$SCORE [match(PHENO.ANY_SN$sjlid, Meningioma_Dobbins$IID)]


## Meningioma_prs.profile
Meningioma <- read.table("Meningioma_prs.profile", header = T)
PHENO.ANY_SN$Meningioma_PRS <-  Meningioma_Dobbins$SCORE [match(PHENO.ANY_SN$sjlid, Meningioma$IID)]


## Pleiotropy_Bi_directional_Increasing_prs.profile
Pleiotropy_Bi_directional_Increasing <- read.table("Pleiotropy_Bi_directional_Increasing_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Increasing_PRS <-  Pleiotropy_Bi_directional_Increasing$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Increasing$IID)]


## Pleiotropy_Bi_directional_Increasing_Significant_prs.profile
Pleiotropy_Bi_directional_Increasing_Sig <- read.table("Pleiotropy_Bi_directional_Increasing_Significant_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Increasing_Sig_PRS <-  Pleiotropy_Bi_directional_Increasing_Sig$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Increasing_Sig$IID)]


## Pleiotropy_Bi_directional_Decreasing_prs.profile
Pleiotropy_Bi_directional_Decreasing <- read.table("Pleiotropy_Bi_directional_Decreasing_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Decreasing_PRS <-  Pleiotropy_Bi_directional_Decreasing$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Decreasing$IID)]


## Pleiotropy_Bi_directional_Decreasing_Significant_prs.profile
Pleiotropy_Bi_directional_Decreasing_Sig <- read.table("Pleiotropy_Bi_directional_Decreasing_Significant_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Decreasing_Sig_PRS <-  Pleiotropy_Bi_directional_Decreasing_Sig$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Decreasing_Sig$IID)]


## Pleiotropy_Meta_analysis_prs.profile
Pleiotropy_Meta_analysis <- read.table("Pleiotropy_Meta_analysis_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Meta_analysis_PRS <-  Pleiotropy_Meta_analysis$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Meta_analysis$IID)]

## Pleiotropy_One_cohort_prs.profile
Pleiotropy_One_cohort <- read.table("Pleiotropy_One_cohort_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_One_cohort_PRS <-  Pleiotropy_One_cohort$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_One_cohort$IID)]


## Pleiotropy_PRSWEB_prs.profile
Pleiotropy_PRSWEB <- read.table("Pleiotropy_PRSWEB_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_PRSWEB_PRS <-  Pleiotropy_PRSWEB$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_PRSWEB$IID)]


## Pleiotropy_One_directional_prs.profile
Pleiotropy_One_directional <- read.table("Pleiotropy_One_directional_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_One_directional_PRS <-  Pleiotropy_One_directional$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_One_directional$IID)]


## Pleiotropy_One_directional_Significant_prs.profile
Pleiotropy_One_directional_Significant <- read.table("Pleiotropy_One_directional_Significant_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_One_directional_Significant_PRS <-  Pleiotropy_One_directional_Significant$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_One_directional_Significant$IID)]


## Pleiotropy_Replication_prior_studies_prs.profile
Pleiotropy_Replication_prior_studies <- read.table("Pleiotropy_Replication_prior_studies_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Replication_prior_studies_PRS <-  Pleiotropy_Replication_prior_studies$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Replication_prior_studies$IID)]


## Sarcoma_Machiela_prs.profile
Sarcoma_Machiela <- read.table("Sarcoma_Machiela_prs.profile", header = T)
PHENO.ANY_SN$Sarcoma_Machiela_PRS <-  Sarcoma_Machiela$SCORE [match(PHENO.ANY_SN$sjlid, Sarcoma_Machiela$IID)]


## ALL_Vijayakrishnan_prs.profile
ALL_Vijayakrishnan <- read.table("ALL_Vijayakrishnan_prs.profile", header = T)
PHENO.ANY_SN$ALL_Vijayakrishnan_PRS <-  ALL_Vijayakrishnan$SCORE [match(PHENO.ANY_SN$sjlid, ALL_Vijayakrishnan$IID)]


## Allman_African_Breast_prs.profile
Allman_African_Breast <- read.table("Allman_African_Breast_prs.profile", header = T)
PHENO.ANY_SN$Allman_African_Breast_PRS <-  Allman_African_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Allman_African_Breast$IID)]

# Allman_Hispanic_Breast_prs.profile
Allman_Hispanic_Breast <- read.table("Allman_Hispanic_Breast_prs.profile", header = T)
PHENO.ANY_SN$Allman_Hispanic_Breast_PRS <-  Allman_Hispanic_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Allman_Hispanic_Breast$IID)]

## Khera_2018_Breast_prs.profile
Khera_2018_Breast <- read.table("Khera_2018_Breast_prs.profile", header = T)
PHENO.ANY_SN$Khera_2018_Breast_PRS <-  Khera_2018_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Khera_2018_Breast$IID)]


## Mavaddat_2015_ER_NEG_Breast_prs.profile
Mavaddat_2015_ER_NEG_Breast <- read.table("Mavaddat_2015_ER_NEG_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2015_ER_NEG_Breast_PRS <-  Mavaddat_2015_ER_NEG_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2015_ER_NEG_Breast$IID)]


## Mavaddat_2015_ER_POS_Breast_prs.profile
Mavaddat_2015_ER_POS_Breast <- read.table("Mavaddat_2015_ER_POS_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2015_ER_POS_Breast_PRS <-  Mavaddat_2015_ER_POS_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2015_ER_POS_Breast$IID)]


## Mavaddat_2015_ER_OVERALL_Breast_prs.profile
Mavaddat_2015_ER_OVERALL_Breast <- read.table("Mavaddat_2015_ER_OVERALL_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2015_ER_OVERALL_Breast_PRS <-  Mavaddat_2015_ER_OVERALL_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2015_ER_OVERALL_Breast$IID)]


## Mavaddat_2019_ER_NEG_Breast_prs.profile
Mavaddat_2019_ER_NEG_Breast <- read.table("Mavaddat_2019_ER_NEG_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_NEG_Breast_PRS <-  Mavaddat_2019_ER_NEG_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2019_ER_NEG_Breast$IID)]


## Mavaddat_2019_ER_POS_Breast_prs.profile
Mavaddat_2019_ER_POS_Breast <- read.table("Mavaddat_2019_ER_POS_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_POS_Breast_PRS <-  Mavaddat_2019_ER_POS_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2019_ER_POS_Breast$IID)]


## Mavaddat_2019_ER_OVERALL_Breast_prs.profile
Mavaddat_2019_ER_OVERALL_Breast <- read.table("Mavaddat_2019_ER_OVERALL_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_OVERALL_Breast_PRS <-  Mavaddat_2019_ER_OVERALL_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2019_ER_OVERALL_Breast$IID)]


## MichiganWeb_ER_NEG_Breast_prs.profile
MichiganWeb_ER_NEG_Breast <- read.table("MichiganWeb_ER_NEG_Breast_prs.profile", header = T)
PHENO.ANY_SN$MichiganWeb_ER_NEG_Breast_PRS <-  MichiganWeb_ER_NEG_Breast$SCORE [match(PHENO.ANY_SN$sjlid, MichiganWeb_ER_NEG_Breast$IID)]


## MichiganWeb_ER_POS_Breast_prs.profile
MichiganWeb_ER_POS_Breast <- read.table("MichiganWeb_ER_POS_Breast_prs.profile", header = T)
PHENO.ANY_SN$MichiganWeb_ER_POS_Breast_PRS <-  MichiganWeb_ER_POS_Breast$SCORE [match(PHENO.ANY_SN$sjlid, MichiganWeb_ER_POS_Breast$IID)]


## MichiganWeb_ER_OVERALL_Breast_prs.profile
MichiganWeb_ER_OVERALL_Breast <- read.table("MichiganWeb_ER_OVERALL_Breast_prs.profile", header = T)
PHENO.ANY_SN$MichiganWeb_ER_OVERALL_Breast_PRS <-  MichiganWeb_ER_OVERALL_Breast$SCORE [match(PHENO.ANY_SN$sjlid, MichiganWeb_ER_OVERALL_Breast$IID)]


## Wang_African_Breast_prs.profile
Wang_African_Breast <- read.table("Wang_African_Breast_prs.profile", header = T)
PHENO.ANY_SN$Wang_African_Breast_PRS <-  Wang_African_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Wang_African_Breast$IID)]

colnames(PHENO.ANY_SN)

saved.PHENO.ANY_SN <- PHENO.ANY_SN

# PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "MRN", "gender", "dob", "diagdt", "agedx", "diaggrp", "agelstcontact", "AnyRT", "anyrt_5", "brainrt_yn", "maxsegrtdose", "chestrt_yn", "maxchestrtdose", "anthra_jco_dose_any", "anthra_jco_dose_5","neckrt_yn", 
#                                "maxneckrtdose", "pelvisrt_yn", "maxpelvisrtdose","abdomenrt_yn", "maxabdrtdose", "aa_class_dose_any", "aa_class_dose_5", "epitxn_dose_any", "epitxn_dose_5",
#                                "cisplat_dose_any", "cisplateq_dose_5", "aa_hvymtl_dose_any", "aa_hvymtl_dose_5", "Zhaoming_carriers", "Qin_carriers", "Qin_carriers.HR.pathways", "Qin_carriers.FA.pathways",
#                                "Qin_carriers.MMR.pathways", "Qin_carriers.BER.pathways", "Qin_carriers.NER.pathways", "Qin_carriers.NHEJ.pathways", "PCA.ethnicity")]





# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores.RDATA")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores.RDATA")

## Removing sample that has no genetic data SJL5450006
PHENO.ANY_SN <- PHENO.ANY_SN[!grepl("SJL5450006", PHENO.ANY_SN$sjlid),]
## Change PRS to categories
PRS.to.categorize <- colnames(PHENO.ANY_SN)[grepl("_PRS$", colnames(PHENO.ANY_SN))]

# ## Breaks are not unique for "Meningiom_Claus_PRS"
# PRS.to.categorize <- PRS.to.categorize[!grepl("Meningiom_Claus_PRS", PRS.to.categorize)]

## Tertile categories
for(i in 1:length(PRS.to.categorize)){
TERT = unname(quantile(PHENO.ANY_SN[PRS.to.categorize[i]][PHENO.ANY_SN[PRS.to.categorize[i]] !=0], c(1/3, 2/3, 1), na.rm = T))
if(sum(duplicated(TERT)) > 0) next
print(TERT)
print (i)
# PHENO.ANY_SN$tmp.tert.category[PHENO.ANY_SN[PRS.to.categorize[i]] ==0| is.na(PHENO.ANY_SN[,PRS.to.categorize[i]])] <- "None"
# PHENO.ANY_SN$tmp.tert.category[is.na(PHENO.ANY_SN[,PRS.to.categorize[i]])] <- "None"
# PHENO.ANY_SN$tmp.tert.category[!grepl("None", PHENO.ANY_SN$tmp.tert.category)] <- as.character(cut(PHENO.ANY_SN[,PRS.to.categorize[i]][!grepl("None", PHENO.ANY_SN$tmp.tert.category)], breaks = c(0, TERT),
#                                            labels = c("1st", "2nd", "3rd"),
#                                            include.lowest = TRUE))

PHENO.ANY_SN$tmp.tert.category <- as.character(cut(PHENO.ANY_SN[,PRS.to.categorize[i]], breaks = c(0, TERT),
                                                                                                   labels = c("1st", "2nd", "3rd"),
                                                                                                   include.lowest = TRUE))

PHENO.ANY_SN$tmp.tert.category <- factor(PHENO.ANY_SN$tmp.tert.category, levels = c("1st", "2nd", "3rd"))
colnames(PHENO.ANY_SN)[colnames(PHENO.ANY_SN) == "tmp.tert.category"] <- paste0(PRS.to.categorize[i], ".tertile.category")
}


## Decile categories
## Breaks are not unique for Meningioma_Dobbins_PRS
for(i in 1:length(PRS.to.categorize)){
TERT = unname(quantile(PHENO.ANY_SN[PRS.to.categorize[i]][PHENO.ANY_SN[PRS.to.categorize[i]] !=0], c(1/10, 1/9, 1), na.rm = T))
if(sum(duplicated(TERT)) > 0) next
print(TERT)
print (i)
# PHENO.ANY_SN$tmp.tert.category[PHENO.ANY_SN[PRS.to.categorize[i]] ==0| is.na(PHENO.ANY_SN[,PRS.to.categorize[i]])] <- "None"
# PHENO.ANY_SN$tmp.tert.category[is.na(PHENO.ANY_SN[,PRS.to.categorize[i]])] <- "None"
# PHENO.ANY_SN$tmp.tert.category[!grepl("None", PHENO.ANY_SN$tmp.tert.category)] <- as.character(cut(PHENO.ANY_SN[,PRS.to.categorize[i]][!grepl("None", PHENO.ANY_SN$tmp.tert.category)], breaks = c(0, TERT),
#                                                                                                    labels = c("1st", "2nd", "3rd"),
#                                                                                                    include.lowest = TRUE))

PHENO.ANY_SN$tmp.tert.category <- as.character(cut(PHENO.ANY_SN[,PRS.to.categorize[i]], breaks = c(0, TERT),
                                                                                                   labels = c("1st", "2nd", "3rd"),
                                                                                                   include.lowest = TRUE))

PHENO.ANY_SN$tmp.tert.category <- factor(PHENO.ANY_SN$tmp.tert.category, levels = c("1st", "2nd", "3rd"))
colnames(PHENO.ANY_SN)[colnames(PHENO.ANY_SN) == "tmp.tert.category"] <- paste0(PRS.to.categorize[i], ".decile.category")
}


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_categories.RDATA")
