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

## Meningioma_prs.profile
Meningioma <- read.table("Meningioma_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Meningioma_PRS <-  Meningioma$SCORE [match(PHENO.ANY_SN$sjlid, Meningioma$IID)]


## Pleiotropy_PRSWEB_prs.profile
Pleiotropy_PRSWEB <- read.table("Pleiotropy_PRSWEB_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_PRSWEB_PRS <-  Pleiotropy_PRSWEB$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_PRSWEB$IID)]


## Sarcoma_Machiela_prs.profile
Sarcoma_Machiela <- read.table("Sarcoma_Machiela_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Sarcoma_Machiela_PRS <-  Sarcoma_Machiela$SCORE [match(PHENO.ANY_SN$sjlid, Sarcoma_Machiela$IID)]


## ALL_Vijayakrishnan_prs.profile
ALL_Vijayakrishnan <- read.table("ALL_Vijayakrishnan_prs.profile", header = T)
PHENO.ANY_SN$ALL_Vijayakrishnan_PRS <-  ALL_Vijayakrishnan$SCORE [match(PHENO.ANY_SN$sjlid, ALL_Vijayakrishnan$IID)]


## Mavaddat_2019_ER_NEG_Breast_prs.profile
Mavaddat_2019_ER_NEG_Breast <- read.table("Mavaddat_2019_ER_NEG_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_NEG_Breast_PRS <-  Mavaddat_2019_ER_NEG_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2019_ER_NEG_Breast$IID)]


## Mavaddat_2019_ER_POS_Breast_prs.profile
Mavaddat_2019_ER_POS_Breast <- read.table("Mavaddat_2019_ER_POS_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_POS_Breast_PRS <-  Mavaddat_2019_ER_POS_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2019_ER_POS_Breast$IID)]


## Mavaddat_2019_ER_OVERALL_Breast_prs.profile
Mavaddat_2019_ER_OVERALL_Breast <- read.table("Mavaddat_2019_ER_OVERALL_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_OVERALL_Breast_PRS <-  Mavaddat_2019_ER_OVERALL_Breast$SCORE [match(PHENO.ANY_SN$sjlid, Mavaddat_2019_ER_OVERALL_Breast$IID)]


## Thyroid cancer
THYROID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/thyroid_preQC_VCF/sjlife_thyroid.profile", header = T)
PHENO.ANY_SN$Thyroid_PRS <-  THYROID$SCORE [match(PHENO.ANY_SN$sjlid, THYROID$IID)]

## NMSCs
BASALcell <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/prs_out/Basal_cell_carcinoma_PRSWeb_prs.profile", header = T)
PHENO.ANY_SN$BASALcell_PRS <-  BASALcell$SCORE [match(PHENO.ANY_SN$sjlid, BASALcell$IID)]


SQUAMOUScell <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/prs_out/Squamous_cell_carcinoma_PRSWeb_prs.profile", header = T)
PHENO.ANY_SN$SQUAMOUScell_PRS <-  SQUAMOUScell$SCORE [match(PHENO.ANY_SN$sjlid, SQUAMOUScell$IID)]


colnames(PHENO.ANY_SN)

saved.PHENO.ANY_SN <- PHENO.ANY_SN

# PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "MRN", "gender", "dob", "diagdt", "agedx", "diaggrp", "agelstcontact", "AnyRT", "anyrt_5", "brainrt_yn", "maxsegrtdose", "chestrt_yn", "maxchestrtdose", "anthra_jco_dose_any", "anthra_jco_dose_5","neckrt_yn", 
#                                "maxneckrtdose", "pelvisrt_yn", "maxpelvisrtdose","abdomenrt_yn", "maxabdrtdose", "aa_class_dose_any", "aa_class_dose_5", "epitxn_dose_any", "epitxn_dose_5",
#                                "cisplat_dose_any", "cisplateq_dose_5", "aa_hvymtl_dose_any", "aa_hvymtl_dose_5", "Zhaoming_carriers", "Qin_carriers", "Qin_carriers.HR.pathways", "Qin_carriers.FA.pathways",
#                                "Qin_carriers.MMR.pathways", "Qin_carriers.BER.pathways", "Qin_carriers.NER.pathways", "Qin_carriers.NHEJ.pathways", "PCA.ethnicity")]





# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_v11.RDATA")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_v11.RDATA")

## Removing sample that has no genetic data SJL5450006
PHENO.ANY_SN <- PHENO.ANY_SN[!grepl("SJL5450006", PHENO.ANY_SN$sjlid),]
## Change PRS to categories
PRS.to.categorize <- colnames(PHENO.ANY_SN)[grepl("_PRS$", colnames(PHENO.ANY_SN))]


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


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_categories_v11.RDATA")
