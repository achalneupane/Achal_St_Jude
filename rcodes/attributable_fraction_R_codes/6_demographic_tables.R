#################################
## 1....................SJLIFE ##
#################################
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


as.data.frame(table(PHENO.ANY_SN$diaggrp))
# as.data.frame(table(PHENO.ANY_SN$diag))

###############
## DIAGNOSIS ##
###############

## Acute lymphoblastic leukemia
ALL <- sum(grepl("lymphoblastic", PHENO.ANY_SN$diaggrp, ignore.case = T))
AML <- sum(grepl("Acute myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
other.leukemia <- sum(grepl("Other leukemia|Chronic myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))

## CNS
CNS.tumors <- PHENO.ANY_SN[PHENO.ANY_SN$diaggrp == "Central nervous system (CNS)",]
CNS.counts <- nrow(CNS.tumors)
astrocytoma_glioma <- sum(grepl("astrocytoma|glioma", CNS.tumors$diag, ignore.case = T))
medulloblastoma_PNET <- sum(grepl("Medulloblastoma|PNET", CNS.tumors$diag, ignore.case = T))
ependymoma <- sum(grepl("Ependymoma", CNS.tumors$diag, ignore.case = T))
other_CNS.tumor <- CNS.counts - (astrocytoma_glioma + medulloblastoma_PNET + ependymoma)

## Lymphoma
hodgkin.lymphoma <- sum(grepl("^Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
non.hodgkin.lymphoma <- sum(grepl("non-Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T)) 

## Sarcoma
sarcoma <- sum(grepl("sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
ewing.sarcoma <- sum(grepl("ewing", PHENO.ANY_SN$diaggrp, ignore.case = T))
osteosarcoma <- sum(grepl("osteosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
rhabdomyosarcoma <- sum(grepl("^Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
nonrhabdomyosarcoma <- sum(grepl("Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))


## Embryonal
wilms <- sum(grepl("wilms", PHENO.ANY_SN$diaggrp, ignore.case = T))
neuroblastoma <- sum(grepl("neuroblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
germ.cell.tumor <- sum(grepl("Germ cell tumor", PHENO.ANY_SN$diaggrp, ignore.case = T))
retinoblastoma <- sum(grepl("Retinoblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
hepatoblastoma <- sum(grepl("Liver malignancies", PHENO.ANY_SN$diaggrp, ignore.case = T))

## Other
melanoma <- sum(grepl("melanoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
carcinoma <- sum(grepl("carcinoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# Others

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
alkylating <- sum(PHENO.ANY_SN$aa_class_dose_5_yn == "Y", na.rm = T)
anthracyclines <- sum(PHENO.ANY_SN$anthra_jco_dose_5_yn == "Y", na.rm = T)
epipodophyllotoxins <- sum(PHENO.ANY_SN$epitxn_dose_5_yn == "Y", na.rm = T)

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
PHENO.ANY_SN$length.followup <- PHENO.ANY_SN$agelstcontact -PHENO.ANY_SN$agedx
median.length.followup <- round(median(PHENO.ANY_SN$length.followup), 1)
length.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$length.followup, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
length.followup.IQR <- gsub(" ", "-", length.followup.IQR)
lenght.followup <- paste0(median.length.followup, " (", length.followup.IQR, ")")


## rbind all
df <- as.data.frame(t(cbind.data.frame(Male, Female, White, Black, Other, ALL, AML, other.leukemia, astrocytoma_glioma, medulloblastoma_PNET, ependymoma, other_CNS.tumor, hodgkin.lymphoma, non.hodgkin.lymphoma, ewing.sarcoma, osteosarcoma, rhabdomyosarcoma, nonrhabdomyosarcoma, wilms, neuroblastoma, germ.cell.tumor, retinoblastoma, hepatoblastoma, melanoma, carcinoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup, lenght.followup)))
df$percent <- c(round((as.numeric(df$V1[1:33])/4401)*100, 1), NA, NA, NA)


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

## Acute lymphoblastic leukemia
ALL <- sum(grepl("lymphoblastic", PHENO.ANY_SN$diaggrpm, ignore.case = T))
AML <- sum(grepl("Acute myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
other.leukemia <- sum(grepl("Other leukemia|Chronic myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))

## CNS
CNS.tumors <- PHENO.ANY_SN[PHENO.ANY_SN$diaggrp == "Central nervous system (CNS)",]
CNS.tumors <- as.data.frame(sort(table(CNS.tumors$diag), decreasing = T))
CNS.counts <- nrow(CNS.tumors)
astrocytoma_glioma <- sum(grepl("astrocytoma|glioma", CNS.tumors$diag, ignore.case = T))
medulloblastoma_PNET <- sum(grepl("Medulloblastoma|PNET", CNS.tumors$diag, ignore.case = T))
ependymoma <- sum(grepl("Ependymoma", CNS.tumors$diag, ignore.case = T))
other_CNS.tumor <- CNS.counts - (astrocytoma_glioma + medulloblastoma_PNET + ependymoma)

## Lymphoma
hodgkin.lymphoma <- sum(grepl("^Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
non.hodgkin.lymphoma <- sum(grepl("non-Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T)) 

## Sarcoma
sarcoma <- sum(grepl("sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
ewing.sarcoma <- sum(grepl("ewing", PHENO.ANY_SN$diaggrp, ignore.case = T))
osteosarcoma <- sum(grepl("osteosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
rhabdomyosarcoma <- sum(grepl("^Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
nonrhabdomyosarcoma <- sum(grepl("Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))


## Embryonal
wilms <- sum(grepl("wilms", PHENO.ANY_SN$diaggrp, ignore.case = T))
neuroblastoma <- sum(grepl("neuroblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
germ.cell.tumor <- sum(grepl("Germ cell tumor", PHENO.ANY_SN$diaggrp, ignore.case = T))
retinoblastoma <- sum(grepl("Retinoblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
hepatoblastoma <- sum(grepl("Liver malignancies", PHENO.ANY_SN$diaggrp, ignore.case = T))

## Other
melanoma <- sum(grepl("melanoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
carcinoma <- sum(grepl("carcinoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# Others

##################
## Radiotherapy ##
##################
brainRT <- sum(PHENO.ANY_SN$brainorheadrt_yn == "Y", na.rm = T)
neckRT <- sum(PHENO.ANY_SN$neckrt_yn == "Y", na.rm = T)
chestRT <- sum(PHENO.ANY_SN$chestrt_yn == "Y", na.rm = T)
abdomenRT <- sum(PHENO.ANY_SN$abdomenrt_yn == "Y", na.rm = T)
pelvisRT <- sum(PHENO.ANY_SN$pelvisrt_yn == "Y", na.rm = T)

##################
## Chemotherapy ##
##################
alkylating <- sum(PHENO.ANY_SN$aa_class_dose_5_yn == "Y", na.rm = T)
anthracyclines <- sum(PHENO.ANY_SN$anthra_jco_dose_5_yn == "Y", na.rm = T)
epipodophyllotoxins <- sum(PHENO.ANY_SN$epitxn_dose_5_yn == "Y", na.rm = T)

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
PHENO.ANY_SN$length.followup <- PHENO.ANY_SN$agelstcontact -PHENO.ANY_SN$agedx
median.length.followup <- round(median(PHENO.ANY_SN$length.followup), 1)
length.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$length.followup, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
length.followup.IQR <- gsub(" ", "-", length.followup.IQR)
lenght.followup <- paste0(median.length.followup, " (", length.followup.IQR, ")")

