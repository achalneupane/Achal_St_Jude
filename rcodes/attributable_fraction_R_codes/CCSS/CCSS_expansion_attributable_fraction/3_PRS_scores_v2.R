## Read CCSS Phenotype data
PHENO.ANY_SN <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data/ExportedCCSS_data.txt", header = T, sep = "\t")

#########################################################
## Read PRS score files and merge them to PHENO.ANY_SN ##
#########################################################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out")
# setwd("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/prs_out")


## Meningioma_prs.profile
Meningioma <- read.table("Meningioma_prs.profile", header = T)
sum(Meningioma$IID %in% PHENO.ANY_SN$ccssid)
# 2936
Meningioma$IID[!Meningioma$IID %in% PHENO.ANY_SN$ccssid]
# 0 # mimatch ID


## Keep pheno that have CCSS EXP WGS data only
PHENO.ANY_SN <- PHENO.ANY_SN[PHENO.ANY_SN$ccssid %in% Meningioma$IID,]
dim(PHENO.ANY_SN)
# [1] 2936   26

PHENO.ANY_SN$Meningioma_PRS <-  Meningioma$SCORE [match(PHENO.ANY_SN$ccssid, Meningioma$IID)]


## Pleiotropy_PRSWEB_prs.profile
Pleiotropy_PRSWEB <- read.table("Pleiotropy_PRSWEB_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_PRSWEB_PRS <-  Pleiotropy_PRSWEB$SCORE [match(PHENO.ANY_SN$ccssid, Pleiotropy_PRSWEB$IID)]


## Sarcoma_Machiela_prs.profile
Sarcoma_Machiela <- read.table("Sarcoma_Machiela_prs.profile", header = T)
PHENO.ANY_SN$Sarcoma_Machiela_PRS <-  Sarcoma_Machiela$SCORE [match(PHENO.ANY_SN$ccssid, Sarcoma_Machiela$IID)]


## ALL_Vijayakrishnan_prs.profile
ALL_Vijayakrishnan <- read.table("ALL_Vijayakrishnan_prs.profile", header = T)
PHENO.ANY_SN$ALL_Vijayakrishnan_PRS <-  ALL_Vijayakrishnan$SCORE [match(PHENO.ANY_SN$ccssid, ALL_Vijayakrishnan$IID)]


## Mavaddat_2019_ER_NEG_Breast_prs.profile
Mavaddat_2019_ER_NEG_Breast <- read.table("Mavaddat_2019_ER_NEG_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_NEG_Breast_PRS <-  Mavaddat_2019_ER_NEG_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_NEG_Breast$IID)]


## Mavaddat_2019_ER_POS_Breast_prs.profile
Mavaddat_2019_ER_POS_Breast <- read.table("Mavaddat_2019_ER_POS_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_POS_Breast_PRS <-  Mavaddat_2019_ER_POS_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_POS_Breast$IID)]


## Mavaddat_2019_ER_OVERALL_Breast_prs.profile
Mavaddat_2019_ER_OVERALL_Breast <- read.table("Mavaddat_2019_ER_OVERALL_Breast_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_OVERALL_Breast_PRS <-  Mavaddat_2019_ER_OVERALL_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_OVERALL_Breast$IID)]


## Thyroid cancer
THYROID <- read.table("THYROID_PGS_prs.profile", header = T)
PHENO.ANY_SN$Thyroid_PRS <-  THYROID$SCORE [match(PHENO.ANY_SN$ccssid, THYROID$IID)]

## NMSCs
BASALcell <- read.table("Basal_cell_carcinoma_PRSWeb_prs.profile", header = T)
PHENO.ANY_SN$BASALcell_PRS <-  BASALcell$SCORE [match(PHENO.ANY_SN$ccssid, BASALcell$IID)]


SQUAMOUScell <- read.table("Squamous_cell_carcinoma_PRSWeb_prs.profile", header = T)
PHENO.ANY_SN$SQUAMOUScell_PRS <-  SQUAMOUScell$SCORE [match(PHENO.ANY_SN$ccssid, SQUAMOUScell$IID)]


colnames(PHENO.ANY_SN)

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


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_categories_CCSS_expansion.RDATA")
