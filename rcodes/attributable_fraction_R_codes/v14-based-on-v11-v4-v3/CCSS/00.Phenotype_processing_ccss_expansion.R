setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2")
rm()

CCSS_data <- read.delim("ExportedCCSS_data.txt", header = T, sep = "\t", stringsAsFactors = F)


BMI.PA.SMK.DRK <- read.delim("ExportedCCSS_BMI_PA_Smk_drink.txt", header = T, sep = "\t", stringsAsFactors = F)

BMI.PA.SMK.DRK <- BMI.PA.SMK.DRK[c("ccssid", "a_base", "a_fu1", "a_fu2", "a_fu3", "a_fu2007", "a_fu5", "a_fu6",
                          "cbmi_0",  "cbmi_2", "cbmi_2007", "cbmi_5",
                          "t_eqmodfu2", "cdc_fu2", "t_eqmodfu5", "cdc_fu5", "t_eqmodfu6", "cdc_fu6", 
                          "smkcatb", "smkcatf2", "smkcatf07", "smkcatf5",
                          "riskyb", "riskyf7", "riskyfu5")]

colnames(BMI.PA.SMK.DRK) <- c("ccssid", "base.age", "fu1.age", "fu2.age", "fu3.age", "fu7.age", "fu5.age", "fu6.age",
                          "base.bmi",  "fu2.bmi", "fu7.bmi", "fu5.bmi",
                          "fu2.MET", "fu2.CDC", "fu5.MET", "fu5.CDC", "fu6.MET", "fu6.CDC", 
                          "base.smk", "fu2.smk", "fu7.smk", "fu5.smk",
                          "base.riskydrk", "fu7.riskydrk", "fu5.riskydrk")


# Make columns uniform for all follow-ups by adding missing columns
BMI.PA.SMK.DRK$base.MET <- NA
BMI.PA.SMK.DRK$base.CDC <- NA

BMI.PA.SMK.DRK$fu1.bmi <- NA
BMI.PA.SMK.DRK$fu1.MET <- NA
BMI.PA.SMK.DRK$fu1.CDC <- NA
BMI.PA.SMK.DRK$fu1.smk <- NA
BMI.PA.SMK.DRK$fu1.riskydrk <- NA

BMI.PA.SMK.DRK$fu2.riskydrk <- NA

BMI.PA.SMK.DRK$fu3.bmi <- NA
BMI.PA.SMK.DRK$fu3.bmi <- NA
BMI.PA.SMK.DRK$fu3.MET <- NA
BMI.PA.SMK.DRK$fu3.CDC <- NA
BMI.PA.SMK.DRK$fu3.smk <- NA
BMI.PA.SMK.DRK$fu3.riskydrk <- NA

BMI.PA.SMK.DRK$fu7.MET <- NA
BMI.PA.SMK.DRK$fu7.CDC <- NA

BMI.PA.SMK.DRK$fu6.bmi <- NA
BMI.PA.SMK.DRK$fu6.smk <- NA
BMI.PA.SMK.DRK$fu6.riskydrk <- NA

BMI.PA.SMK.DRK[BMI.PA.SMK.DRK == "."] <- NA

## Reshape BMI.PA.SMK.DRK to long format
# dw <- BMI.PA.SMK.DRK[1:2,]

## OR
# suppressPackageStartupMessages({
#   library(tidyr)
# })
#  dw |>
#   pivot_longer(
#     cols = -ccssid,
#     names_to = c("var", ".value"),
#     names_pattern = "(.*)\\.(.*)"
#   )

## make prefixes to suffixes
names(BMI.PA.SMK.DRK) <- strsplit(names(BMI.PA.SMK.DRK), '\\.') |> lapply(rev) |> sapply(paste, collapse='.')

cc <- reshape(BMI.PA.SMK.DRK, direction='l', idvar='ccssid', varying=sort(names(BMI.PA.SMK.DRK)[-1]))




cc$age <- as.numeric(cc$age)

## BMI
bmi_iid_dob_18 = subset(cc, age >= 18)
bmi_iid_dob_18 <- bmi_iid_dob_18[!is.na(bmi_iid_dob_18$bmi), ]
bmi_iid_dob_18_sorted = bmi_iid_dob_18[order(bmi_iid_dob_18$ccssid, bmi_iid_dob_18$age, decreasing = FALSE),]
bmi_iid_dob_18_uniq = bmi_iid_dob_18_sorted[!duplicated(bmi_iid_dob_18_sorted$ccssid),]

## MET
MET_iid_dob_18 = subset(cc, age >= 18)
MET_iid_dob_18 <- MET_iid_dob_18[!is.na(MET_iid_dob_18$MET), ]
MET_iid_dob_18_sorted = MET_iid_dob_18[order(MET_iid_dob_18$ccssid, MET_iid_dob_18$age, decreasing = FALSE),]
MET_iid_dob_18_uniq = MET_iid_dob_18_sorted[!duplicated(MET_iid_dob_18_sorted$ccssid),]

## smk
smk_iid_dob_18 = subset(cc, age >= 18)
smk_iid_dob_18 <- smk_iid_dob_18[!is.na(smk_iid_dob_18$smk), ]
smk_iid_dob_18_sorted = smk_iid_dob_18[order(smk_iid_dob_18$ccssid, smk_iid_dob_18$age, decreasing = FALSE),]
smk_iid_dob_18_uniq = smk_iid_dob_18_sorted[!duplicated(smk_iid_dob_18_sorted$ccssid),]

## drinking
drk_iid_dob_18 = subset(cc, age >= 18)
drk_iid_dob_18 <- drk_iid_dob_18[!is.na(drk_iid_dob_18$riskydrk), ]
drk_iid_dob_18_sorted = drk_iid_dob_18[order(drk_iid_dob_18$ccssid, drk_iid_dob_18$age, decreasing = FALSE),]
drk_iid_dob_18_uniq = drk_iid_dob_18_sorted[!duplicated(drk_iid_dob_18_sorted$ccssid),]


##############
## CCSS_exp ##
##############
ccss_exp.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/merged.dat.fam")
sum(CCSS_data$ccssid %in% ccss_exp.samples$V2)
# 3078
sum(ccss_exp.samples$V2 %in% CCSS_data$ccssid)
# 2936

#############################
## Add lifestyle variables ##
#############################
ccss_exp <- CCSS_data[CCSS_data$ccssid %in% ccss_exp.samples$V2,]


# Obesity
ccss_exp$BMI <- as.numeric(bmi_iid_dob_18_uniq$bmi[match(ccss_exp$ccssid, bmi_iid_dob_18_uniq$ccssid)])
ccss_exp$Not_obese_yn_agesurvey <- as.numeric(bmi_iid_dob_18_uniq$age[match(ccss_exp$ccssid, bmi_iid_dob_18_uniq$ccssid)])
ccss_exp$Not_obese_yn <- factor(ifelse(as.numeric(bmi_iid_dob_18_uniq$bmi[match(ccss_exp$ccssid, bmi_iid_dob_18_uniq$ccssid)]) < 30, 1, 0))
ccss_exp$Not_obese_yn <- factor(ccss_exp$Not_obese_yn, level = c(1, 0, "Unknown")) 
ccss_exp$Not_obese_yn[is.na(ccss_exp$Not_obese_yn)] <- "Unknown";

# Physical activity
ccss_exp$CDC <- MET_iid_dob_18_uniq$CDC[match(ccss_exp$ccssid, MET_iid_dob_18_uniq$ccssid)]
ccss_exp$PhysicalActivity_yn_agesurvey <- as.numeric(MET_iid_dob_18_uniq$age[match(ccss_exp$ccssid, MET_iid_dob_18_uniq$ccssid)])
ccss_exp$PhysicalActivity_yn <- factor(MET_iid_dob_18_uniq$CDC[match(ccss_exp$ccssid, MET_iid_dob_18_uniq$ccssid)])
ccss_exp$PhysicalActivity_yn <- ifelse (ccss_exp$PhysicalActivity_yn == "Yes", 1, 0)
ccss_exp$PhysicalActivity_yn[is.na(ccss_exp$PhysicalActivity_yn)] <- "Unknown"
ccss_exp$PhysicalActivity_yn <- factor(ccss_exp$PhysicalActivity_yn, level = c(1, 0, "Unknown")) 

# Smoker
ccss_exp$SMK <- smk_iid_dob_18_uniq$smk[match(ccss_exp$ccssid, smk_iid_dob_18_uniq$ccssid)]
ccss_exp$smoker_former_or_never_yn_agesurvey <- as.numeric(smk_iid_dob_18_uniq$age[match(ccss_exp$ccssid, smk_iid_dob_18_uniq$ccssid)])
ccss_exp$smoker_former_or_never_yn <- factor(smk_iid_dob_18_uniq$smk[match(ccss_exp$ccssid, smk_iid_dob_18_uniq$ccssid)])
ccss_exp$smoker_former_or_never_yn <- factor(ifelse(ccss_exp$smoker_former_or_never_yn != 3, 1, 0))
ccss_exp$smoker_former_or_never_yn <- factor(ccss_exp$smoker_former_or_never_yn, level = c(1, 0, "Unknown")) 
ccss_exp$smoker_former_or_never_yn[is.na(ccss_exp$smoker_former_or_never_yn)] <- "Unknown"

# drinker
ccss_exp$DRK <- drk_iid_dob_18_uniq$riskydrk[match(ccss_exp$ccssid, smk_iid_dob_18_uniq$ccssid)]
ccss_exp$NOT_RiskyHeavyDrink_yn_agesurvey <- as.numeric(drk_iid_dob_18_uniq$age[match(ccss_exp$ccssid, smk_iid_dob_18_uniq$ccssid)])
ccss_exp$NOT_RiskyHeavyDrink_yn <- factor(ifelse(factor(drk_iid_dob_18_uniq$riskydrk[match(ccss_exp$ccssid, smk_iid_dob_18_uniq$ccssid)]) == "No", 1, 0))
ccss_exp$NOT_RiskyHeavyDrink_yn <- factor(ccss_exp$NOT_RiskyHeavyDrink_yn, level = c(1, 0, "Unknown")) 
ccss_exp$NOT_RiskyHeavyDrink_yn[is.na(ccss_exp$NOT_RiskyHeavyDrink_yn)] <- "Unknown"

#############################
## Harmonize age variables ##
#############################
## Use this variable ccss_exp_ANY_SN for ANY_SN
ccss_exp$AGE.ANY_SN <- as.numeric(ccss_exp$a_candx)
ccss_exp$agedx <- as.numeric(ccss_exp$a_dx)
ccss_exp$agelstcontact <- as.numeric(ccss_exp$a_end)


## Age at diagnosis
# PHENO.ANY_SN$agedx <- floor(PHENO.ANY_SN$agedx)
ccss_exp$AGE_AT_DIAGNOSIS[ccss_exp$agedx >= 0 & ccss_exp$agedx < 5 ] <- "0-4"
ccss_exp$AGE_AT_DIAGNOSIS[ccss_exp$agedx >= 5 & ccss_exp$agedx < 10 ] <- "5-9"
ccss_exp$AGE_AT_DIAGNOSIS[ccss_exp$agedx >= 10 & ccss_exp$agedx < 15 ] <- "10-14"
ccss_exp$AGE_AT_DIAGNOSIS[ccss_exp$agedx >= 15 ] <- ">=15"
ccss_exp$AGE_AT_DIAGNOSIS <- factor(ccss_exp$AGE_AT_DIAGNOSIS, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference


## Age at last contact
ccss_exp$AGE_AT_LAST_CONTACT[ccss_exp$agelstcontact >= 0 & ccss_exp$agelstcontact < 25 ] <- "0-24"
ccss_exp$AGE_AT_LAST_CONTACT[ccss_exp$agelstcontact >= 25 & ccss_exp$agelstcontact < 35 ] <- "25-34"
ccss_exp$AGE_AT_LAST_CONTACT[ccss_exp$agelstcontact >= 35 & ccss_exp$agelstcontact < 45 ] <- "35-44"
ccss_exp$AGE_AT_LAST_CONTACT[ccss_exp$agelstcontact >= 45 ] <- ">=45"
ccss_exp$AGE_AT_LAST_CONTACT <- factor(ccss_exp$AGE_AT_LAST_CONTACT, levels = c("0-24", "25-34", "35-44", ">=45")) # first level will be treated as reference

ccss_exp$gradedt <- as.numeric(ccss_exp$a_candx)


subneo <- ccss_exp
PHENO.ANY_SN <- ccss_exp[c('ccssid', 'SEX', 'agedx', 'd_candx', 'diagnose', 'AGE_AT_DIAGNOSIS', 'agelstcontact', 
                           'chestrtgrp', 'neckrtgrp', 'abdomenrtgrp', 'abdomenrtgrp', 'brainrtgrp', 'pelvisrtgrp', 
                           'smoker_former_or_never_yn_agesurvey', 'smoker_former_or_never_yn', 'NOT_RiskyHeavyDrink_yn_agesurvey', 'NOT_RiskyHeavyDrink_yn',
                           'Not_obese_yn_agesurvey', 'Not_obese_yn', 'PhysicalActivity_yn_agesurvey', 'PhysicalActivity_yn',
                           'chestmaxrtdose', 'neckmaxrtdose', 'pelvismaxrtdose', 'abdmaxrtdose', 'maxsegrtdose', 'anth_DED5', 'alk_CED5', 'epipdose5', 'pt_cisED5')]
PHENO.ANY_SN <- PHENO.ANY_SN[!duplicated(PHENO.ANY_SN$ccssid),]



add_cubic_spline <- function(ccss_exp){
  ## Age at last contact (cubic spline)
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/cubic_spline.r")
  
breaks = seq(5, 95, 22.5)

cp = quantile(ccss_exp$agelstcontact, breaks/100, na.rm = T)

cs = cubic_spline(ccss_exp$agelstcontact, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
ccss_exp <- cbind.data.frame(ccss_exp, cs)
# Merge cs to your original data.frame and adjust in the logistic regression
# add cubic spline
return(ccss_exp)
}




####################################################
## Phenotype: radiation and chemotherapy tertiles ##
####################################################
## Create tertiles for ccss org and ccss exp separately

add_therapy_tertiles <- function(ccss_exp){
## Anthracyclines (Y/N and Tertiles)
ccss_exp$anth_DED5 <- as.numeric(ccss_exp$anth_DED5)
ccss_exp$anthra_jco_dose_5_yn <- factor(ifelse(ccss_exp$anth_DED5 == 0, "N", "Y"))

TERT = unname(quantile(ccss_exp$anth_DED5[ccss_exp$anth_DED5 !=0], c(1/3, 2/3, 1), na.rm = T))
ccss_exp$anthra_jco_dose_5.category <- cut(ccss_exp$anth_DED5, breaks = c(0, 0.001, TERT),
                                               labels = c("None", "1st", "2nd", "3rd"),
                                               include.lowest = TRUE)
levels(ccss_exp$anthra_jco_dose_5.category) <- c(levels(ccss_exp$anthra_jco_dose_5.category), "Unknown")
ccss_exp$anthra_jco_dose_5.category [is.na(ccss_exp$anthra_jco_dose_5.category)] <- "Unknown"


## Alkylating agents (Y/N and Tertiles)
ccss_exp$alk_CED5 <- as.numeric(ccss_exp$alk_CED5)
ccss_exp$aa_class_dose_5_yn <- factor(ifelse(ccss_exp$alk_CED5 == 0, "N", "Y"))

TERT = unname(quantile(ccss_exp$alk_CED5[ccss_exp$alk_CED5 !=0], c(1/3, 2/3, 1), na.rm = T))
ccss_exp$aa_class_dose_5.category <- cut(ccss_exp$alk_CED5, breaks = c(0, 0.001, TERT),
                                           labels = c("None", "1st", "2nd", "3rd"),
                                           include.lowest = TRUE)
levels(ccss_exp$aa_class_dose_5.category) <- c(levels(ccss_exp$aa_class_dose_5.category), "Unknown")
ccss_exp$aa_class_dose_5.category [is.na(ccss_exp$aa_class_dose_5.category)] <- "Unknown"

## Epipodophyllotoxin agents (Y/N and Tertiles)
ccss_exp$epipdose5 <- as.numeric(ccss_exp$epipdose5)
ccss_exp$epitxn_dose_5_yn <- factor(ifelse(ccss_exp$epipdose5 == 0, "N", "Y"))

TERT = unname(quantile(ccss_exp$epipdose5[ccss_exp$epipdose5 !=0], c(1/3, 2/3, 1), na.rm = T))
ccss_exp$epitxn_dose_5.category <- cut(ccss_exp$epipdose5, breaks = c(0, 0.001, TERT),
                                         labels = c("None", "1st", "2nd", "3rd"),
                                         include.lowest = TRUE)
levels(ccss_exp$epitxn_dose_5.category) <- c(levels(ccss_exp$epitxn_dose_5.category), "Unknown")
ccss_exp$epitxn_dose_5.category [is.na(ccss_exp$epitxn_dose_5.category)] <- "Unknown"

## Cis-platinum (Y/N and Tertiles)
ccss_exp$pt_cisED5 <- as.numeric(ccss_exp$pt_cisED5)
ccss_exp$cisplateq_dose_5_yn <- factor(ifelse(ccss_exp$pt_cisED5 == 0, "N", "Y"))

TERT = unname(quantile(ccss_exp$pt_cisED5[ccss_exp$pt_cisED5 !=0], c(1/3, 2/3, 1), na.rm = T))
ccss_exp$cisplateq_dose_5.category <- cut(ccss_exp$pt_cisED5, breaks = c(0, 0.001, TERT),
                                         labels = c("None", "1st", "2nd", "3rd"),
                                         include.lowest = TRUE)
levels(ccss_exp$cisplateq_dose_5.category) <- c(levels(ccss_exp$cisplateq_dose_5.category), "Unknown")
ccss_exp$cisplateq_dose_5.category [is.na(ccss_exp$cisplateq_dose_5.category)] <- "Unknown"

# Making variable names consistent with the SJLIFE variables
ccss_exp$gender <- factor(ccss_exp$SEX, levels = c("Male", "Female"))
ccss_exp$maxchestrtdose.category <- factor(ccss_exp$chestrtgrp, levels = c("None", "0-20", ">=20", "Unknown"))
ccss_exp$maxneckrtdose.category <- factor(ccss_exp$neckrtgrp, levels = c("None", "0-11", "11-20", "20-30", ">=30", "Unknown"))
ccss_exp$maxabdrtdose.category <- factor(ccss_exp$abdomenrtgrp, levels = c("None", "0-30", ">=30", "Unknown"))
ccss_exp$maxsegrtdose.category <- factor(ccss_exp$brainrtgrp, levels = c("None", "0-18", "18-30", ">=30", "Unknown"))
ccss_exp$maxpelvisrtdose.category <- factor(ccss_exp$pelvisrtgrp, levels = c("None", "0-20", ">=20", "Unknown"))

# table(PHENO.ANY_SN$cisplateq_dose_5.category)
# table(ccss_exp$cisplateq_dose_5.category)
return(ccss_exp)
}



# # Now KEEP the unique PHENO samples
# sum(duplicated(PHENO.ANY_SN$ccssid))
# # 142
# PHENO.ANY_SN <- PHENO.ANY_SN[!duplicated(PHENO.ANY_SN$ccssid),]




################
## PRS scores ##
################
add_PRS_to_PHENO <-function(PHENO.ANY_SN){
PHENO.ANY_SN <- PHENO.ANY_SN[c('ccssid', 'gender', 'agedx', 'd_candx', 'diagnose', 'agelstcontact', 'AGE_AT_DIAGNOSIS', 
                                 "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", 
                                 'maxchestrtdose.category', 'maxneckrtdose.category', 'maxabdrtdose.category', 'maxsegrtdose.category', 'maxpelvisrtdose.category',
                                 'Not_obese_yn_agesurvey', 'Not_obese_yn', 'PhysicalActivity_yn_agesurvey', 'PhysicalActivity_yn', 
                                 'smoker_former_or_never_yn_agesurvey', 'smoker_former_or_never_yn', 'NOT_RiskyHeavyDrink_yn_agesurvey', 'NOT_RiskyHeavyDrink_yn', 
                                 'anthra_jco_dose_5.category', 'aa_class_dose_5.category', 'epitxn_dose_5.category', 'cisplateq_dose_5.category',
                               'chestmaxrtdose', 'neckmaxrtdose', 'pelvismaxrtdose', 'abdmaxrtdose', 'maxsegrtdose', 'anth_DED5', 'alk_CED5', 'epipdose5', 'pt_cisED5')]
## Meningioma_from_variants_also_in_CCSS_org_prs.profile-----------------------------------------------
Meningioma <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Meningioma_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Meningioma_PRS <-  Meningioma$SCORE [match(PHENO.ANY_SN$ccssid, Meningioma$IID)]

## Pleiotropy_PRSWEB_from_variants_also_in_CCSS_org_prs.profile--------- 
Pleiotropy_PRSWEB <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Pleiotropy_PRSWEB_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_PRSWEB_PRS <-  Pleiotropy_PRSWEB$SCORE [match(PHENO.ANY_SN$ccssid, Pleiotropy_PRSWEB$IID)]


## Sarcoma_Machiela_from_variants_also_in_CCSS_org_prs.profile-----------------------------------------
Sarcoma_Machiela <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Sarcoma_Machiela_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Sarcoma_Machiela_PRS <-  Sarcoma_Machiela$SCORE [match(PHENO.ANY_SN$ccssid, Sarcoma_Machiela$IID)]


## ALL_Vijayakrishnan_prs.profile---------------------------------------
ALL_Vijayakrishnan <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/ALL_Vijayakrishnan_prs.profile", header = T)
PHENO.ANY_SN$ALL_Vijayakrishnan_PRS <-  ALL_Vijayakrishnan$SCORE [match(PHENO.ANY_SN$ccssid, ALL_Vijayakrishnan$IID)]


## Breast cancer--------------------------------------------------------
## Mavaddat_2019_ER_NEG_Breast_from_variants_also_in_CCSS_org_prs.profile
Mavaddat_2019_ER_NEG_Breast <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Mavaddat_2019_ER_NEG_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_NEG_Breast_PRS <-  Mavaddat_2019_ER_NEG_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_NEG_Breast$IID)]


## Mavaddat_2019_ER_POS_Breast_from_variants_also_in_CCSS_org_prs.profile
Mavaddat_2019_ER_POS_Breast <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Mavaddat_2019_ER_POS_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_POS_Breast_PRS <-  Mavaddat_2019_ER_POS_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_POS_Breast$IID)]


## Mavaddat_2019_ER_OVERALL_Breast_from_variants_also_in_CCSS_org_prs.profile
Mavaddat_2019_ER_OVERALL_Breast <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Mavaddat_2019_ER_OVERALL_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
PHENO.ANY_SN$Mavaddat_2019_ER_OVERALL_Breast_PRS <-  Mavaddat_2019_ER_OVERALL_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_OVERALL_Breast$IID)]


## Thyroid cancer------------------------
THYROID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/THYROID_PGS_prs.profile", header = T)
PHENO.ANY_SN$Thyroid_PRS <-  THYROID$SCORE [match(PHENO.ANY_SN$ccssid, THYROID$IID)]

## NMSCs---------------------------------
BASALcell <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Basal_cell_carcinoma_PRSWeb_prs.profile", header = T)
PHENO.ANY_SN$BASALcell_PRS <-  BASALcell$SCORE [match(PHENO.ANY_SN$ccssid, BASALcell$IID)]

SQUAMOUScell <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Squamous_cell_carcinoma_PRSWeb_prs.profile", header = T)
PHENO.ANY_SN$SQUAMOUScell_PRS <-  SQUAMOUScell$SCORE [match(PHENO.ANY_SN$ccssid, SQUAMOUScell$IID)]


## Change PRS to categories
PRS.to.categorize <- colnames(PHENO.ANY_SN)[grepl("_PRS$", colnames(PHENO.ANY_SN))]


## Tertile categories
for(i in 1:length(PRS.to.categorize)){
  TERT = unname(quantile(PHENO.ANY_SN[PRS.to.categorize[i]][PHENO.ANY_SN[PRS.to.categorize[i]] !=0], c(1/3, 2/3, 1), na.rm = T))
  if(sum(duplicated(TERT)) > 0) next
  print(TERT)
  print (i)
  PHENO.ANY_SN$tmp.tert.category <- as.character(cut(PHENO.ANY_SN[,PRS.to.categorize[i]], breaks = c(0, TERT),
                                                     labels = c("1st", "2nd", "3rd"),
                                                     include.lowest = TRUE))
  
  PHENO.ANY_SN$tmp.tert.category <- factor(PHENO.ANY_SN$tmp.tert.category, levels = c("1st", "2nd", "3rd"))
  colnames(PHENO.ANY_SN)[colnames(PHENO.ANY_SN) == "tmp.tert.category"] <- paste0(PRS.to.categorize[i], ".tertile.category")
}


return(PHENO.ANY_SN)
}

PHENO.ANY_SN <- add_cubic_spline(PHENO.ANY_SN)
PHENO.ANY_SN <- add_therapy_tertiles(PHENO.ANY_SN)
PHENO.ANY_SN <- add_PRS_to_PHENO(PHENO.ANY_SN)

rm(list=ls()[!grepl(c("PHENO.ANY_SN|subneo"), ls())])


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_exp_Genetic_data_P_LP_v14.Rdata")

