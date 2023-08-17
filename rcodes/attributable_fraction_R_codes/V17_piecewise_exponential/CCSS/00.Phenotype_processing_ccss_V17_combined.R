rm(list=ls())

##############
## CCSS_exp ##
##############
ccss_exp.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/merged.dat.fam")

##############
## CCSS_org ##
##############
ccss_org.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/plink_data/merged_batch2.fam")
ccss_org.samples$V1 <- sub("^(\\d+)_.*", "\\1", ccss_org.samples$V1)

ccss_samples <- c(ccss_exp.samples$V1, ccss_org.samples$V1)
overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt", header = T)

sum(ccss_samples %in% overlaps$ccssid)
ccss_samples <- ccss_samples[!(ccss_samples %in% overlaps$ccssid)]

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2")

CCSS_data <- read.delim("ExportedCCSS_data.txt", header = T, sep = "\t", stringsAsFactors = F)
CCSS_data <- CCSS_data[CCSS_data$ccssid %in% ccss_samples,]

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

## Keep only those present in ccss WGS
BMI.PA.SMK.DRK <- BMI.PA.SMK.DRK[BMI.PA.SMK.DRK$ccssid %in% ccss_samples,]

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
cc <- subset(cc, age >= 18)
dim(cc)

table(is.na(cc$MET))
table(is.na(cc$smk))
table(is.na(cc$riskydrk))
table(is.na(cc$bmi))


######################
## recode lifestyle ##
######################
## BMI
bmi_iid_dob_18 = subset(cc, age >= 18)
bmi_iid_dob_18 <- bmi_iid_dob_18[!is.na(bmi_iid_dob_18$bmi), ]
bmi_iid_dob_18_sorted = bmi_iid_dob_18[order(bmi_iid_dob_18$ccssid, bmi_iid_dob_18$age, decreasing = FALSE),]
bmi_iid_dob_18_uniq = bmi_iid_dob_18_sorted[!duplicated(bmi_iid_dob_18_sorted$ccssid),]
bmi_iid_dob_18_uniq = bmi_iid_dob_18_uniq[,c("ccssid", "time", "age", "bmi")]

## MET
MET_iid_dob_18 = subset(cc, age >= 18)
MET_iid_dob_18 <- MET_iid_dob_18[!is.na(MET_iid_dob_18$MET), ]
MET_iid_dob_18_sorted = MET_iid_dob_18[order(MET_iid_dob_18$ccssid, MET_iid_dob_18$age, decreasing = FALSE),]
MET_iid_dob_18_uniq = MET_iid_dob_18_sorted[!duplicated(MET_iid_dob_18_sorted$ccssid),]
MET_iid_dob_18_uniq = MET_iid_dob_18_uniq[,c("ccssid", "time", "age", "CDC")]

## smk
smk_iid_dob_18 = subset(cc, age >= 18)
smk_iid_dob_18 <- smk_iid_dob_18[!is.na(smk_iid_dob_18$smk), ]
smk_iid_dob_18_sorted = smk_iid_dob_18[order(smk_iid_dob_18$ccssid, smk_iid_dob_18$age, decreasing = FALSE),]
smk_iid_dob_18_uniq = smk_iid_dob_18_sorted[!duplicated(smk_iid_dob_18_sorted$ccssid),]
smk_iid_dob_18_uniq = smk_iid_dob_18_uniq[,c("ccssid", "time", "age", "smk")]

## drinking
drk_iid_dob_18 = subset(cc, age >= 18)
drk_iid_dob_18 <- drk_iid_dob_18[!is.na(drk_iid_dob_18$riskydrk), ]
drk_iid_dob_18_sorted = drk_iid_dob_18[order(drk_iid_dob_18$ccssid, drk_iid_dob_18$age, decreasing = FALSE),]
drk_iid_dob_18_uniq = drk_iid_dob_18_sorted[!duplicated(drk_iid_dob_18_sorted$ccssid),]
drk_iid_dob_18_uniq = drk_iid_dob_18_uniq[,c("ccssid", "time", "age", "riskydrk")]

merged_df <- Reduce(function(x, y) left_join(x, y, by = 'ccssid'), list(bmi_iid_dob_18_uniq, smk_iid_dob_18_uniq, drk_iid_dob_18_uniq, MET_iid_dob_18_uniq))

test <- merged_df[grepl("age", colnames(merged_df))]
colnames(test) <- c("age_bmi", "age_smk", "age_drk", "age_PA")

count_same <- function(row) {
  sum(row == row[1], na.rm = TRUE)
}

# Apply the function row-wise to the dataframe
test$same_count <- apply(test, 1, count_same)
table(test$same_count)

table(merged_df$ccssid %in% ccss_samples)

test <- test[c("age_bmi", "age_smk", "age_drk")]
test$same_count <- apply(test, 1, count_same)
table(test$same_count)

#############################
## Add lifestyle variables ##
#############################



# Obesity
CCSS_data$BMI <- as.numeric(bmi_iid_dob_18_uniq$bmi[match(CCSS_data$ccssid, bmi_iid_dob_18_uniq$ccssid)])
CCSS_data$Obese_yn_agesurvey <- as.numeric(bmi_iid_dob_18_uniq$age[match(CCSS_data$ccssid, bmi_iid_dob_18_uniq$ccssid)])
CCSS_data$Obese_yn <- factor(ifelse(CCSS_data$BMI < 30, "No", "Yes"))
CCSS_data$Obese_yn <- factor(CCSS_data$Obese_yn, level = c("No", "Yes", "Unknown")) 
CCSS_data$Obese_yn[is.na(CCSS_data$Obese_yn)] <- "Unknown";

# Physical activity
CCSS_data$CDC <- MET_iid_dob_18_uniq$CDC[match(CCSS_data$ccssid, MET_iid_dob_18_uniq$ccssid)]
CCSS_data$PhysicalActivity_yn_agesurvey <- as.numeric(MET_iid_dob_18_uniq$age[match(CCSS_data$ccssid, MET_iid_dob_18_uniq$ccssid)])
CCSS_data$PhysicalActivity_yn <- factor(CCSS_data$CDC)
CCSS_data$PhysicalActivity_yn <- ifelse (CCSS_data$PhysicalActivity_yn == "Yes", "Yes", "No")
CCSS_data$PhysicalActivity_yn[is.na(CCSS_data$PhysicalActivity_yn)] <- "Unknown"
CCSS_data$PhysicalActivity_yn <- factor(CCSS_data$PhysicalActivity_yn, level = c("Yes", "No", "Unknown")) 

# Smoker
CCSS_data$SMK <- smk_iid_dob_18_uniq$smk[match(CCSS_data$ccssid, smk_iid_dob_18_uniq$ccssid)]
CCSS_data$Current_smoker_yn_agesurvey <- as.numeric(smk_iid_dob_18_uniq$age[match(CCSS_data$ccssid, smk_iid_dob_18_uniq$ccssid)])
CCSS_data$Current_smoker_yn <- factor(CCSS_data$SMK)
CCSS_data$Current_smoker_yn <- factor(ifelse(CCSS_data$Current_smoker_yn != 3, "No", "Yes"))
CCSS_data$Current_smoker_yn <- factor(CCSS_data$Current_smoker_yn, level = c("No", "Yes", "Unknown")) 
CCSS_data$Current_smoker_yn[is.na(CCSS_data$Current_smoker_yn)] <- "Unknown"

# drinker
CCSS_data$DRK <- drk_iid_dob_18_uniq$riskydrk[match(CCSS_data$ccssid, drk_iid_dob_18_uniq$ccssid)]
CCSS_data$RiskyHeavyDrink_yn_agesurvey <- as.numeric(drk_iid_dob_18_uniq$age[match(CCSS_data$ccssid, drk_iid_dob_18_uniq$ccssid)])
CCSS_data$RiskyHeavyDrink_yn <- factor(ifelse(factor(CCSS_data$DRK) == "No", "No", "Yes"))
CCSS_data$RiskyHeavyDrink_yn <- factor(CCSS_data$RiskyHeavyDrink_yn, level = c("No", "Yes", "Unknown")) 
CCSS_data$RiskyHeavyDrink_yn[is.na(CCSS_data$RiskyHeavyDrink_yn)] <- "Unknown"


# # remove those with all 4 lifestyle missing
# CCSS_data <- CCSS_data[!(is.na(CCSS_data$Obese_yn_agesurvey) & is.na(CCSS_data$Current_smoker_yn_agesurvey) &
#                            is.na(CCSS_data$RiskyHeavyDrink_yn_agesurvey) & is.na(CCSS_data$PhysicalActivity_yn_agesurvey)), ] ##$$


test <- CCSS_data[grepl("agesurvey", colnames(CCSS_data))]
# > colnames(test)
# [1] "Obese_yn_agesurvey"            "PhysicalActivity_yn_agesurvey" "Current_smoker_yn_agesurvey"   "RiskyHeavyDrink_yn_agesurvey" 

count_same <- function(row) {
  sum(row == row[1], na.rm = TRUE)
}

# Apply the function row-wise to the dataframe
test$same_count <- apply(test, 1, count_same)
table(test$same_count)

test <- CCSS_data[!duplicated(CCSS_data$ccssid),]
test <- test[c("Obese_yn_agesurvey", "Current_smoker_yn_agesurvey", "RiskyHeavyDrink_yn_agesurvey")]

test$same_count <- apply(test, 1, count_same)
table(test$same_count)


#############################
## Harmonize age variables ##
#############################
## Use this variable CCSS_data_ANY_SN for ANY_SN
CCSS_data$AGE.ANY_SN <- as.numeric(CCSS_data$a_candx)
CCSS_data$agedx <- as.numeric(CCSS_data$a_dx)
CCSS_data$agelstcontact <- as.numeric(CCSS_data$a_end)


## Age at diagnosis
# PHENO.ANY_SN$agedx <- floor(PHENO.ANY_SN$agedx)
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 0 & CCSS_data$agedx < 5 ] <- "0-4"
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 5 & CCSS_data$agedx < 10 ] <- "5-9"
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 10 & CCSS_data$agedx < 15 ] <- "10-14"
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 15 ] <- ">=15"
CCSS_data$AGE_AT_DIAGNOSIS <- factor(CCSS_data$AGE_AT_DIAGNOSIS, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference


## Age at last contact
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 0 & CCSS_data$agelstcontact < 25 ] <- "0-24"
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 25 & CCSS_data$agelstcontact < 35 ] <- "25-34"
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 35 & CCSS_data$agelstcontact < 45 ] <- "35-44"
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 45 ] <- ">=45"
CCSS_data$AGE_AT_LAST_CONTACT <- factor(CCSS_data$AGE_AT_LAST_CONTACT, levels = c("0-24", "25-34", "35-44", ">=45")) # first level will be treated as reference

CCSS_data$gradedt <- as.numeric(CCSS_data$a_candx)


subneo <- CCSS_data
PHENO.ANY_SN <- CCSS_data[c('ccssid', 'SEX', 'agedx', 'AGE_AT_DIAGNOSIS', 'agelstcontact', 
                           'chestrtgrp', 'neckrtgrp', 'abdomenrtgrp', 'abdomenrtgrp', 'brainrtgrp', 'pelvisrtgrp', 
                           'chestmaxrtdose', 'neckmaxrtdose', 'pelvismaxrtdose', 'abdmaxrtdose', 'maxsegrtdose', 'anth_DED5', 'alk_CED5', 'epipdose5', 'pt_cisED5', 
                           "Obese_yn_agesurvey", "Obese_yn", "PhysicalActivity_yn_agesurvey", "PhysicalActivity_yn", "Current_smoker_yn_agesurvey", "Current_smoker_yn", "RiskyHeavyDrink_yn_agesurvey", "RiskyHeavyDrink_yn")]
PHENO.ANY_SN <- PHENO.ANY_SN[!duplicated(PHENO.ANY_SN$ccssid),]



add_cubic_spline <- function(CCSS_data){
  ## Age at last contact (cubic spline)
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/cubic_spline.r")
  
breaks = seq(5, 95, 22.5)

cp = quantile(CCSS_data$agelstcontact, breaks/100, na.rm = T)

cs = cubic_spline(CCSS_data$agelstcontact, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
CCSS_data <- cbind.data.frame(CCSS_data, cs)
# Merge cs to your original data.frame and adjust in the logistic regression
# add cubic spline
return(CCSS_data)
}




####################################################
## Phenotype: radiation and chemotherapy tertiles ##
####################################################
## Create tertiles for ccss org and ccss exp separately

add_therapy_tertiles <- function(CCSS_data){
## Anthracyclines (Y/N and Tertiles)
CCSS_data$anth_DED5 <- as.numeric(CCSS_data$anth_DED5)
CCSS_data$anthra_jco_dose_5_yn <- factor(ifelse(CCSS_data$anth_DED5 == 0, "N", "Y"))

TERT = unname(quantile(CCSS_data$anth_DED5[CCSS_data$anth_DED5 !=0], c(1/3, 2/3, 1), na.rm = T))
CCSS_data$anthra_jco_dose_5.category <- cut(CCSS_data$anth_DED5, breaks = c(0, 0.001, TERT),
                                               labels = c("None", "1st", "2nd", "3rd"),
                                               include.lowest = TRUE)
levels(CCSS_data$anthra_jco_dose_5.category) <- c(levels(CCSS_data$anthra_jco_dose_5.category), "Unknown")
CCSS_data$anthra_jco_dose_5.category [is.na(CCSS_data$anthra_jco_dose_5.category)] <- "Unknown"


## Alkylating agents (Y/N and Tertiles)
CCSS_data$alk_CED5 <- as.numeric(CCSS_data$alk_CED5)
CCSS_data$aa_class_dose_5_yn <- factor(ifelse(CCSS_data$alk_CED5 == 0, "N", "Y"))

TERT = unname(quantile(CCSS_data$alk_CED5[CCSS_data$alk_CED5 !=0], c(1/3, 2/3, 1), na.rm = T))
CCSS_data$aa_class_dose_5.category <- cut(CCSS_data$alk_CED5, breaks = c(0, 0.001, TERT),
                                           labels = c("None", "1st", "2nd", "3rd"),
                                           include.lowest = TRUE)
levels(CCSS_data$aa_class_dose_5.category) <- c(levels(CCSS_data$aa_class_dose_5.category), "Unknown")
CCSS_data$aa_class_dose_5.category [is.na(CCSS_data$aa_class_dose_5.category)] <- "Unknown"

## Epipodophyllotoxin agents (Y/N and Tertiles)
CCSS_data$epipdose5 <- as.numeric(CCSS_data$epipdose5)
CCSS_data$epitxn_dose_5_yn <- factor(ifelse(CCSS_data$epipdose5 == 0, "N", "Y"))

TERT = unname(quantile(CCSS_data$epipdose5[CCSS_data$epipdose5 !=0], c(1/3, 2/3, 1), na.rm = T))
CCSS_data$epitxn_dose_5.category <- cut(CCSS_data$epipdose5, breaks = c(0, 0.001, TERT),
                                         labels = c("None", "1st", "2nd", "3rd"),
                                         include.lowest = TRUE)
levels(CCSS_data$epitxn_dose_5.category) <- c(levels(CCSS_data$epitxn_dose_5.category), "Unknown")
CCSS_data$epitxn_dose_5.category [is.na(CCSS_data$epitxn_dose_5.category)] <- "Unknown"

## Cis-platinum (Y/N and Tertiles)
CCSS_data$pt_cisED5 <- as.numeric(CCSS_data$pt_cisED5)
CCSS_data$cisplateq_dose_5_yn <- factor(ifelse(CCSS_data$pt_cisED5 == 0, "N", "Y"))

TERT = unname(quantile(CCSS_data$pt_cisED5[CCSS_data$pt_cisED5 !=0], c(1/3, 2/3, 1), na.rm = T))
CCSS_data$cisplateq_dose_5.category <- cut(CCSS_data$pt_cisED5, breaks = c(0, 0.001, TERT),
                                         labels = c("None", "1st", "2nd", "3rd"),
                                         include.lowest = TRUE)
levels(CCSS_data$cisplateq_dose_5.category) <- c(levels(CCSS_data$cisplateq_dose_5.category), "Unknown")
CCSS_data$cisplateq_dose_5.category [is.na(CCSS_data$cisplateq_dose_5.category)] <- "Unknown"

# Making variable names consistent with the SJLIFE variables
CCSS_data$gender <- factor(CCSS_data$SEX, levels = c("Male", "Female"))
CCSS_data$maxchestrtdose.category <- factor(CCSS_data$chestrtgrp, levels = c("None", "0-20", ">=20", "Unknown"))
CCSS_data$maxneckrtdose.category <- factor(CCSS_data$neckrtgrp, levels = c("None", "0-11", "11-20", "20-30", ">=30", "Unknown"))
CCSS_data$maxabdrtdose.category <- factor(CCSS_data$abdomenrtgrp, levels = c("None", "0-30", ">=30", "Unknown"))
CCSS_data$maxsegrtdose.category <- factor(CCSS_data$brainrtgrp, levels = c("None", "0-18", "18-30", ">=30", "Unknown"))
CCSS_data$maxpelvisrtdose.category <- factor(CCSS_data$pelvisrtgrp, levels = c("None", "0-20", ">=20", "Unknown"))

# table(PHENO.ANY_SN$cisplateq_dose_5.category)
# table(CCSS_data$cisplateq_dose_5.category)
return(CCSS_data)
}



# # Now KEEP the unique PHENO samples
# sum(duplicated(PHENO.ANY_SN$ccssid))
# # 142
# PHENO.ANY_SN <- PHENO.ANY_SN[!duplicated(PHENO.ANY_SN$ccssid),]




################
## PRS scores ##
################
add_PRS_to_PHENO <-function(PHENO.ANY_SN){
PHENO.ANY_SN <- PHENO.ANY_SN[c('ccssid', 'gender', 'agelstcontact', 'AGE_AT_DIAGNOSIS', 
                                 "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", 
                                 'maxchestrtdose.category', 'maxneckrtdose.category', 'maxabdrtdose.category', 'maxsegrtdose.category', 'maxpelvisrtdose.category',
                                 'anthra_jco_dose_5.category', 'aa_class_dose_5.category', 'epitxn_dose_5.category', 'cisplateq_dose_5.category',
                               'chestmaxrtdose', 'neckmaxrtdose', 'pelvismaxrtdose', 'abdmaxrtdose', 'maxsegrtdose', 'anth_DED5', 'alk_CED5', 'epipdose5', 'pt_cisED5',
                               "Obese_yn_agesurvey", "Obese_yn", "PhysicalActivity_yn_agesurvey", "PhysicalActivity_yn", "Current_smoker_yn_agesurvey", "Current_smoker_yn", "RiskyHeavyDrink_yn_agesurvey", "RiskyHeavyDrink_yn")]
## Meningioma_from_variants_also_in_CCSS_org_prs.profile-----------------------------------------------
Meningioma.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Meningioma_from_variants_also_in_CCSS_org_prs.profile", header = T)
Meningioma.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Meningioma_prs.profile", header = T)
Meningioma.org$IID <- sub("^(\\d+)_.*", "\\1", Meningioma.org$IID)
Meningioma <- rbind.data.frame(Meningioma.exp, Meningioma.org)
PHENO.ANY_SN$Meningioma_PRS <-  Meningioma$SCORE [match(PHENO.ANY_SN$ccssid, Meningioma$IID)]


## Pleiotropy_PRSWEB_from_variants_also_in_CCSS_org_prs.profile--------- 
Pleiotropy_PRSWEB.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Pleiotropy_PRSWEB_from_variants_also_in_CCSS_org_prs.profile", header = T)
Pleiotropy_PRSWEB.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Meningioma_prs.profile", header = T)
Pleiotropy_PRSWEB.org$IID <- sub("^(\\d+)_.*", "\\1", Pleiotropy_PRSWEB.org$IID)
Pleiotropy_PRSWEB <- rbind.data.frame(Pleiotropy_PRSWEB.exp, Pleiotropy_PRSWEB.org)
PHENO.ANY_SN$Pleiotropy_PRSWEB_PRS <-  Pleiotropy_PRSWEB$SCORE [match(PHENO.ANY_SN$ccssid, Pleiotropy_PRSWEB$IID)]

## Sarcoma_Machiela_from_variants_also_in_CCSS_org_prs.profile-----------------------------------------
Sarcoma_Machiela.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Sarcoma_Machiela_from_variants_also_in_CCSS_org_prs.profile", header = T)
Sarcoma_Machiela.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Sarcoma_Machiela_prs.profile", header = T)
Sarcoma_Machiela.org$IID <- sub("^(\\d+)_.*", "\\1", Sarcoma_Machiela.org$IID)
Sarcoma_Machiela <- rbind.data.frame(Sarcoma_Machiela.exp, Sarcoma_Machiela.org)
PHENO.ANY_SN$Sarcoma_Machiela_PRS <-  Sarcoma_Machiela$SCORE [match(PHENO.ANY_SN$ccssid, Sarcoma_Machiela$IID)]


## ALL_Vijayakrishnan_prs.profile---------------------------------------
ALL_Vijayakrishnan.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/ALL_Vijayakrishnan_prs.profile", header = T)
ALL_Vijayakrishnan.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/ALL_Vijayakrishnan_prs.profile", header = T)
ALL_Vijayakrishnan.org$IID <- sub("^(\\d+)_.*", "\\1", ALL_Vijayakrishnan.org$IID)
ALL_Vijayakrishnan <- rbind.data.frame(ALL_Vijayakrishnan.exp, ALL_Vijayakrishnan.org)
PHENO.ANY_SN$ALL_Vijayakrishnan_PRS <-  ALL_Vijayakrishnan$SCORE [match(PHENO.ANY_SN$ccssid, ALL_Vijayakrishnan$IID)]


## Breast cancer--------------------------------------------------------
## Mavaddat_2019_ER_NEG_Breast_from_variants_also_in_CCSS_org_prs.profile
Mavaddat_2019_ER_NEG_Breast.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Mavaddat_2019_ER_NEG_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
Mavaddat_2019_ER_NEG_Breast.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Mavaddat_2019_ER_NEG_Breast_prs.profile", header = T)
Mavaddat_2019_ER_NEG_Breast.org$IID <- sub("^(\\d+)_.*", "\\1", Mavaddat_2019_ER_NEG_Breast.org$IID)
Mavaddat_2019_ER_NEG_Breast <- rbind.data.frame(Mavaddat_2019_ER_NEG_Breast.exp, Mavaddat_2019_ER_NEG_Breast.org)
PHENO.ANY_SN$Mavaddat_2019_ER_NEG_Breast_PRS <-  Mavaddat_2019_ER_NEG_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_NEG_Breast$IID)]


## Mavaddat_2019_ER_POS_Breast_from_variants_also_in_CCSS_org_prs.profile
Mavaddat_2019_ER_POS_Breast.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Mavaddat_2019_ER_POS_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
Mavaddat_2019_ER_POS_Breast.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Mavaddat_2019_ER_POS_Breast_prs.profile", header = T)
Mavaddat_2019_ER_POS_Breast.org$IID <- sub("^(\\d+)_.*", "\\1", Mavaddat_2019_ER_POS_Breast.org$IID)
Mavaddat_2019_ER_POS_Breast <- rbind.data.frame(Mavaddat_2019_ER_POS_Breast.exp, Mavaddat_2019_ER_POS_Breast.org)
PHENO.ANY_SN$Mavaddat_2019_ER_POS_Breast_PRS <-  Mavaddat_2019_ER_POS_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_POS_Breast$IID)]


## Mavaddat_2019_ER_OVERALL_Breast_from_variants_also_in_CCSS_org_prs.profile
Mavaddat_2019_ER_OVERALL_Breast.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Mavaddat_2019_ER_OVERALL_Breast_from_variants_also_in_CCSS_org_prs.profile", header = T)
Mavaddat_2019_ER_OVERALL_Breast.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Mavaddat_2019_ER_OVERALL_Breast_prs.profile", header = T)
Mavaddat_2019_ER_OVERALL_Breast.org$IID <- sub("^(\\d+)_.*", "\\1", Mavaddat_2019_ER_OVERALL_Breast.org$IID)
Mavaddat_2019_ER_OVERALL_Breast <- rbind.data.frame(Mavaddat_2019_ER_OVERALL_Breast.exp, Mavaddat_2019_ER_OVERALL_Breast.org)
PHENO.ANY_SN$Mavaddat_2019_ER_OVERALL_Breast_PRS <-  Mavaddat_2019_ER_OVERALL_Breast$SCORE [match(PHENO.ANY_SN$ccssid, Mavaddat_2019_ER_OVERALL_Breast$IID)]


## Thyroid cancer------------------------
THYROID.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/THYROID_PGS_prs.profile", header = T)
THYROID.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/THYROID_PGS_prs.profile", header = T)
THYROID.org$IID <- sub("^(\\d+)_.*", "\\1", THYROID.org$IID)
THYROID <- rbind.data.frame(THYROID.exp, THYROID.org)
PHENO.ANY_SN$Thyroid_PRS <-  THYROID$SCORE [match(PHENO.ANY_SN$ccssid, THYROID$IID)]

## NMSCs---------------------------------
BASALcell.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Basal_cell_carcinoma_PRSWeb_prs.profile", header = T)
BASALcell.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Basal_cell_carcinoma_PRSWeb_prs.profile", header = T)
BASALcell.org$IID <- sub("^(\\d+)_.*", "\\1", BASALcell.org$IID)
BASALcell <- rbind.data.frame(BASALcell.exp, BASALcell.org)
PHENO.ANY_SN$BASALcell_PRS <-  BASALcell$SCORE [match(PHENO.ANY_SN$ccssid, BASALcell$IID)]

SQUAMOUScell.exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/prs_out/Squamous_cell_carcinoma_PRSWeb_prs.profile", header = T)
SQUAMOUScell.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/prs_out/Squamous_cell_carcinoma_PRSWeb_prs.profile", header = T)
SQUAMOUScell.org$IID <- sub("^(\\d+)_.*", "\\1", SQUAMOUScell.org$IID)
SQUAMOUScell <- rbind.data.frame(SQUAMOUScell.exp, SQUAMOUScell.org)
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


PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# None 0-18 18-30 >=30 Unknown -->>> None >0-<18 >=18-<30 >=30 Unknown
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "0-18"] <- ">0-<18"
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "18-30"] <- ">=18-<30"
PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<18", ">=18-<30", ">=30", "Unknown"))
table(PHENO.ANY_SN$maxsegrtdose.category)

PHENO.ANY_SN$maxneckrtdose.category <- as.character(PHENO.ANY_SN$maxneckrtdose.category)
# None 0-11 11-20 20-30 >=30 Unknown -->>> None >0-<11 >=11-<20 >=20-<30 >=30 Unknown
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "0-11"] <- ">0-<11"
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "11-20"] <- ">=11-<20"
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "20-30"] <- ">=20-<30"
PHENO.ANY_SN$maxneckrtdose.category <- factor(PHENO.ANY_SN$maxneckrtdose.category, levels = c("None", ">0-<11", ">=11-<20", ">=20-<30", ">=30", "Unknown"))
table(PHENO.ANY_SN$maxneckrtdose.category)

PHENO.ANY_SN$maxabdrtdose.category <- as.character(PHENO.ANY_SN$maxabdrtdose.category)
#  None 0-30 >=30 Unknown -->>> None >0-<30 >=30 Unknown
PHENO.ANY_SN$maxabdrtdose.category[PHENO.ANY_SN$maxabdrtdose.category == "0-30"] <- ">0-<30"
PHENO.ANY_SN$maxabdrtdose.category <- factor(PHENO.ANY_SN$maxabdrtdose.category, levels = c("None", ">0-<30", ">=30", "Unknown"))
table(PHENO.ANY_SN$maxabdrtdose.category)

PHENO.ANY_SN$maxchestrtdose.category <- as.character(PHENO.ANY_SN$maxchestrtdose.category)
#  None 0-20 >=20 Unknown -->>> None >0-<20 >=20 Unknown
PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == "0-20"] <- ">0-<20"
PHENO.ANY_SN$maxchestrtdose.category <- factor(PHENO.ANY_SN$maxchestrtdose.category, levels = c("None", ">0-<20", ">=20", "Unknown"))
table(PHENO.ANY_SN$maxchestrtdose.category)

PHENO.ANY_SN$maxpelvisrtdose.category <- as.character(PHENO.ANY_SN$maxpelvisrtdose.category)
#  None 0-20 >=20 Unknown -->>> None >0-<20 >=20 Unknown
PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == "0-20"] <- ">0-<20"
PHENO.ANY_SN$maxpelvisrtdose.category <- factor(PHENO.ANY_SN$maxpelvisrtdose.category, levels = c("None", ">0-<20", ">=20", "Unknown"))
table(PHENO.ANY_SN$maxpelvisrtdose.category)


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

