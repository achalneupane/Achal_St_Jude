#################################
## 1....................SJLIFE ##
#################################
## Diagnosis variable mapping based on SJLIFE <-> CCSS DXs mapping; Yadav's email: on 05/02/2023


library(haven)
## Load PHenotype
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)



wgsdiag <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/wgsdiag.sas7bdat")
head(wgsdiag)

PHENO.ANY_SN$diag <- wgsdiag$diag[match(PHENO.ANY_SN$MRN, wgsdiag$MRN)]

dim(PHENO.ANY_SN)
# View(PHENO.ANY_SN)

## Sex
table(PHENO.ANY_SN$gender)
Male = sum(PHENO.ANY_SN$gender == "Male")
Female = sum(PHENO.ANY_SN$gender != "Male")

# ## PCA ethnicity
# table(PHENO.ANY_SN$PCA.ethnicity)
# Black = sum(PHENO.ANY_SN$PCA.ethnicity == "AFR")
# White = sum(PHENO.ANY_SN$PCA.ethnicity == "EUR")
# # Asian = sum(PHENO.ANY_SN$PCA.ethnicity == "EAS")
# Other = 4401 - sum(Black, White)

## Ethnicity (Self-reported)
table(PHENO.ANY_SN$ethnic)
non_hispanic = sum(grepl("Non Hispanic", PHENO.ANY_SN$ethnic))
hispanic = sum(!grepl("Non Hispanic", PHENO.ANY_SN$ethnic)) - 37 # (Uknown + caribbean)

## Race (Self-reported)
Black = sum(PHENO.ANY_SN$race == "Black")
White = sum(PHENO.ANY_SN$race == "White")
Other = 4401 - sum(Black, White)

table(PHENO.ANY_SN$ethnic)
# Caribbean                              Cuban                    Mexican/Chicano 
# 5                                  2                                 27 
# Non Spanish speaking, Non Hispanic        NOS Spanish,Hispanic,Latino                       Puerto Rican 
# 4225                                 60                                 10 
# South or Central Amercian                          Southwest                            Unknown 
# 39                                  1                                 32 

as.data.frame(table(PHENO.ANY_SN$diaggrp))
# as.data.frame(table(PHENO.ANY_SN$diag))

###############
## DIAGNOSIS ##
###############

# ## Acute lymphoblastic leukemia
# ALL <- sum(grepl("lymphoblastic", PHENO.ANY_SN$diaggrp, ignore.case = T))
# AML <- sum(grepl("Acute myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
# other.leukemia <- sum(grepl("Other leukemia|Chronic myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
# 
# ## CNS
# CNS.tumors <- PHENO.ANY_SN[PHENO.ANY_SN$diaggrp == "Central nervous system (CNS)",]
# CNS.counts <- nrow(CNS.tumors)
# astrocytoma_glioma <- sum(grepl("astrocytoma|glioma", CNS.tumors$diag, ignore.case = T))
# medulloblastoma_PNET <- sum(grepl("Medulloblastoma|PNET", CNS.tumors$diag, ignore.case = T))
# ependymoma <- sum(grepl("Ependymoma", CNS.tumors$diag, ignore.case = T))
# other_CNS.tumor <- CNS.counts - (astrocytoma_glioma + medulloblastoma_PNET + ependymoma)
# 
# ## Lymphoma
# hodgkin.lymphoma <- sum(grepl("^Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# non.hodgkin.lymphoma <- sum(grepl("non-Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T)) 
# 
# ## Sarcoma
# sarcoma <- sum(grepl("sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# ewing.sarcoma <- sum(grepl("ewing", PHENO.ANY_SN$diaggrp, ignore.case = T))
# osteosarcoma <- sum(grepl("osteosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# rhabdomyosarcoma <- sum(grepl("^Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# nonrhabdomyosarcoma <- sum(grepl("Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# 
# 
# ## Embryonal
# wilms <- sum(grepl("wilms", PHENO.ANY_SN$diaggrp, ignore.case = T))
# neuroblastoma <- sum(grepl("neuroblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# germ.cell.tumor <- sum(grepl("Germ cell tumor", PHENO.ANY_SN$diaggrp, ignore.case = T))
# retinoblastoma <- sum(grepl("Retinoblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# hepatoblastoma <- sum(grepl("Liver malignancies", PHENO.ANY_SN$diaggrp, ignore.case = T))
# 
# ## Other
# melanoma <- sum(grepl("melanoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# carcinoma <- sum(grepl("carcinoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# # Others

## Diagnosis as in CCSS (with common disease types)
bone.cancer <- sum(grepl("Osteosarcoma|Ewing sarcoma family of tumors", PHENO.ANY_SN$diaggrp, ignore.case = T))
CNS <- sum(grepl("Central nervous system", PHENO.ANY_SN$diaggrp, ignore.case = T))
leukemia <- sum(grepl("leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
# leukemia  <- sum(grepl("Acute lymphoblastic leukemia|Acute myeloid leukemia|Other leukemia|MDS/Acute myeloid", PHENO.ANY_SN$diaggrp, ignore.case = T))
neuroblastoma  <- sum(grepl("Neuroblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
wilms  <- sum(grepl("Wilms tumor", PHENO.ANY_SN$diaggrp, ignore.case = T))
HD  <- sum(grepl("^Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
NHL  <- sum(grepl("Non-Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
soft.tissue.sarcoma  <- sum(grepl("Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
Rhabdomyosarcoma  <- sum(grepl("Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
##################
## Radiotherapy ##
##################
# brainRT <- sum(PHENO.ANY_SN$brainorheadrt_yn == "Y", na.rm = T)
# neckRT <- sum(PHENO.ANY_SN$neckrt_yn == "Y", na.rm = T)
# chestRT <- sum(PHENO.ANY_SN$chestrt_yn == "Y", na.rm = T)
# abdomenRT <- sum(PHENO.ANY_SN$abdomenrt_yn == "Y", na.rm = T)
# pelvisRT <- sum(PHENO.ANY_SN$pelvisrt_yn == "Y", na.rm = T)

##
brainRT <- sum(PHENO.ANY_SN$maxsegrtdose > 200, na.rm = T)
neckRT <- sum(PHENO.ANY_SN$maxneckrtdose > 200, na.rm = T)
chestRT <- sum(PHENO.ANY_SN$maxchestrtdose > 200, na.rm = T)
abdomenRT <- sum(PHENO.ANY_SN$maxabdrtdose > 200, na.rm = T)
pelvisRT <- sum(PHENO.ANY_SN$maxpelvisrtdose > 200, na.rm = T)

##################
## Chemotherapy ##
##################
# alkylating <- sum(PHENO.ANY_SN$aa_class_dose_5_yn == "Y", na.rm = T)
alkylating <- sum(PHENO.ANY_SN$aa_class_dose_5 > 0, na.rm = T)
# anthracyclines <- sum(PHENO.ANY_SN$anthra_jco_dose_5_yn == "Y", na.rm = T)
anthracyclines <- sum(PHENO.ANY_SN$anthra_jco_dose_5 > 0, na.rm = T)
# epipodophyllotoxins <- sum(PHENO.ANY_SN$epitxn_dose_5_yn == "Y", na.rm = T)
epipodophyllotoxins <- sum(PHENO.ANY_SN$epitxn_dose_5 > 0, na.rm = T)

## age at diagnosis
median.agedx <- round(median(PHENO.ANY_SN$agedx), 1)
agedx.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agedx, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agedx.IQR <- gsub(" ", "-", agedx.IQR)
agedx <- paste0(median.agedx, " (", agedx.IQR, ")")

## Age at follow up
median.age.followup <- round(median(PHENO.ANY_SN$agelstcontact), 1)
age.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
age.followup.IQR <- gsub(" ", "-", age.followup.IQR)
age.at.followup <- paste0(median.age.followup, " (", age.followup.IQR, ")")

# ## Length of follow up
# PHENO.ANY_SN$length.followup <- PHENO.ANY_SN$agelstcontact -PHENO.ANY_SN$agedx
# median.length.followup <- round(median(PHENO.ANY_SN$length.followup), 1)
# length.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$length.followup, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
# length.followup.IQR <- gsub(" ", "-", length.followup.IQR)
# lenght.followup <- paste0(median.length.followup, " (", length.followup.IQR, ")")


## rbind all
# df <- as.data.frame(t(cbind.data.frame(Male, Female, White, Black, Other, ALL, AML, other.leukemia, astrocytoma_glioma, medulloblastoma_PNET, ependymoma, other_CNS.tumor, hodgkin.lymphoma, non.hodgkin.lymphoma, ewing.sarcoma, osteosarcoma, rhabdomyosarcoma, nonrhabdomyosarcoma, wilms, neuroblastoma, germ.cell.tumor, retinoblastoma, hepatoblastoma, melanoma, carcinoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup, lenght.followup)))


# df.sjlife <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup, lenght.followup)))
df.sjlife <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, Rhabdomyosarcoma, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup)))
df.sjlife$percent <- c(round((as.numeric(df.sjlife$V1[!grepl("agedx|age.at", rownames(df.sjlife))])/4401)*100, 1), NA, NA)


###############################
## 2....................CCSS ##
###############################
## Load PHenotype
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

dim(PHENO.ANY_SN)
# View(PHENO.ANY_SN)

## Sex
table(PHENO.ANY_SN$gender)
Male = sum(PHENO.ANY_SN$gender == "Male")
Female = sum(PHENO.ANY_SN$gender != "Male")

# ## PCA ethnicity
# table(PHENO.ANY_SN$PCA.ethnicity)
# Black = sum(PHENO.ANY_SN$PCA.ethnicity == "AFR")
# White = sum(PHENO.ANY_SN$PCA.ethnicity == "EUR")
# # Asian = sum(PHENO.ANY_SN$PCA.ethnicity == "EAS")
# Other = 4401 - sum(Black, White)

# ## Ethnicity (Self-reported)
# table(PHENO.ANY_SN$ethnic)
# non_hispanic = sum(grepl("Non Hispanic", PHENO.ANY_SN$ethnic))
# hispanic = sum(!grepl("Non Hispanic", PHENO.ANY_SN$ethnic)) - 37 # (Uknown + caribbean)
# 
# ## Race (Self-reported)
# Black = sum(PHENO.ANY_SN$race == "Black")
# White = sum(PHENO.ANY_SN$race == "White")
# Other = 4401 - sum(Black, White)


as.data.frame(table(PHENO.ANY_SN$diagnose))
# as.data.frame(table(PHENO.ANY_SN$diag))

###############
## DIAGNOSIS ##
###############
bone.cancer <- sum(PHENO.ANY_SN$diagnose == "Bone cancer", na.rm = T)
CNS <- sum(PHENO.ANY_SN$diagnose == "CNS", na.rm = T)
HD <- sum(PHENO.ANY_SN$diagnose == "HD", na.rm = T)
wilms <- sum(PHENO.ANY_SN$diagnose == "Kidney (Wilms)", na.rm = T)
leukemia <- sum(PHENO.ANY_SN$diagnose == "Leukemia", na.rm = T)
neuroblastoma <- sum(PHENO.ANY_SN$diagnose == "Neuroblastoma", na.rm = T)
NHL <- sum(PHENO.ANY_SN$diagnose == "NHL", na.rm = T)
NHL <- sum(PHENO.ANY_SN$diagnose == "NHL", na.rm = T)
soft.tissue.sarcoma <- sum(PHENO.ANY_SN$diagnose == "Soft tissue sarcoma", na.rm = T)
Rhabdomyosarcoma  <- sum(grepl("Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
##################
## Radiotherapy ##
##################
PHENO.ANY_SN$maxsegrtdose <- as.numeric(PHENO.ANY_SN$maxsegrtdose)
PHENO.ANY_SN$neckmaxrtdose <- as.numeric(PHENO.ANY_SN$neckmaxrtdose)
PHENO.ANY_SN$chestmaxrtdose <- as.numeric(PHENO.ANY_SN$chestmaxrtdose)
PHENO.ANY_SN$abdmaxrtdose <- as.numeric(PHENO.ANY_SN$abdmaxrtdose)
PHENO.ANY_SN$pelvismaxrtdose <- as.numeric(PHENO.ANY_SN$pelvismaxrtdose)

brainRT <- sum(PHENO.ANY_SN$maxsegrtdose > 0, na.rm = T)
round((brainRT/7943)*100,1)

neckRT <- sum(PHENO.ANY_SN$neckmaxrtdose > 0, na.rm = T)
round((neckRT/7943)*100,1)

chestRT <- sum(PHENO.ANY_SN$chestmaxrtdose > 0, na.rm = T)
round((chestRT/7943)*100,1)

abdomenRT <- sum(PHENO.ANY_SN$abdmaxrtdose > 0, na.rm = T)
round((abdomenRT/7943)*100,1)

pelvisRT <- sum(PHENO.ANY_SN$pelvismaxrtdose > 0, na.rm = T)
round((pelvisRT/7943)*100,1)

##################
## Chemotherapy ##
##################
alkylating <- sum(PHENO.ANY_SN$alk_CED5 > 0, na.rm = T)
round((alkylating/7943)*100,1)

anthracyclines <- sum(PHENO.ANY_SN$anth_DED5 > 0, na.rm = T)
round((anthracyclines/7943)*100,1)

epipodophyllotoxins <- sum(PHENO.ANY_SN$epipdose5 > 0, na.rm = T)
round((epipodophyllotoxins/7943)*100,1)
## age at diagnosis
median.agedx <- round(median(PHENO.ANY_SN$agedx), 1)
agedx.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agedx, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agedx.IQR <- gsub(" ", "-", agedx.IQR)
agedx <- paste0(median.agedx, " (", agedx.IQR, ")")

## Age at follow up
median.age.followup <- round(median(PHENO.ANY_SN$agelstcontact), 1)
age.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
age.followup.IQR <- gsub(" ", "-", age.followup.IQR)
age.at.followup <- paste0(median.age.followup, " (", age.followup.IQR, ")")

## Length of follow up
# PHENO.ANY_SN$length.followup <- PHENO.ANY_SN$agelstcontact -PHENO.ANY_SN$agedx
# median.length.followup <- round(median(PHENO.ANY_SN$length.followup), 1)
# length.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$length.followup, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
# length.followup.IQR <- gsub(" ", "-", length.followup.IQR)
# lenght.followup <- paste0(median.length.followup, " (", length.followup.IQR, ")")


## Common primary dx in CCSS and SJLIFE
# CCSS: SJLIFE
# Bone cancer: Osteosarcoma|Ewing sarcoma family of tumors|
#   CNS: Central nervous system (CNS)
# Leukemia: Acute lymphoblastic leukemia|Acute myeloid leukemia|Other leukemia|MDS/Acute myeloid leukemia|Chronic myeloid leukemia
# Neuroblastoma: Neuroblastoma
# Kidney (Wilms): Wilms tumor
# HD: Hodgkin lymphoma
# NHL: Non-Hodgkin lymphoma
# Soft tissue sarcoma: Soft tissue sarcoma


# df.ccss <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup, lenght.followup)))
df.ccss <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, Rhabdomyosarcoma, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup)))
df.ccss$percent <- c(round((as.numeric(df.ccss$V1[!grepl("agedx|age.at", rownames(df.ccss))])/7943)*100, 1), NA, NA)

colnames(df.sjlife) <- c("SJLIFE (n)", "sjlife%")
colnames(df.ccss) <- c("CCSS (n)", "CCSS%")

df <- cbind(df.sjlife, df.ccss)

