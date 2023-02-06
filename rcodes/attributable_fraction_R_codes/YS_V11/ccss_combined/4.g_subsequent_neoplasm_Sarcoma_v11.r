## Load CCSS org
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_org_Genetic_data_P_LP_V11.Rdata")

library(haven)
library(benchmarkme)
library(dplyr)
library(plyr)
library(data.table)
library (birk)
library(gtools)
library(stringr)
# library(tidyverse)
library(lubridate)
# benchmarkme::get_ram()

table(is.na(PHENO.ANY_SN$smoker_former_or_never_yn))
#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1351

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$ccssid))
# 19


#############
## Sarcoma ##
#############
SARCOMA <- subneo[grepl("sarcoma", subneo$groupdx3, ignore.case = T),]
SARCOMA <- setDT(SARCOMA)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(SARCOMA)
# 89
# Removing samples with SNs within 5 years of childhood cancer
SARCOMA <- SARCOMA[!SARCOMA$ccssid %in% subneo.within5$ccssid,]
nrow(SARCOMA)
# 86
PHENO.ANY_SN$SARCOMA <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% SARCOMA$ccssid, 0, 1))
table(PHENO.ANY_SN$SARCOMA)
# 0    1 
# 4921   86 


PHENO.ANY_SN$AGE.ANY_SN <- SARCOMA$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, SARCOMA$ccssid)]
PHENO.ANY_SN$ANY_SN_TYPE <- SARCOMA$groupdx3[match(PHENO.ANY_SN$ccssid, SARCOMA$ccssid)]
#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
PHENO.ANY_SN$CACO <- PHENO.ANY_SN$SARCOMA

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "smoker_former_or_never_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "PhysicalActivity_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "NOT_RiskyHeavyDrink_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "Not_obese_yn"] <- "Unknown"

# SARCOMA.2 <- PHENO.ANY_SN[PHENO.ANY_SN$SARCOMA == 1,]
#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


## Load CCSS exp
CCSS.org.ANY_SN <- PHENO.ANY_SN # CCSS org




#########################
## Subsequent Neoplasm ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_exp_Genetic_data_P_LP_V11.Rdata")

#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 269

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$ccssid))
# 8

#############
## Sarcoma ##
#############
SARCOMA <- subneo[grepl("sarcoma", subneo$groupdx3, ignore.case = T),]
SARCOMA <- setDT(SARCOMA)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(SARCOMA)
# 20
# Removing samples with SNs within 5 years of childhood cancer
SARCOMA <- SARCOMA[!SARCOMA$ccssid %in% subneo.within5$ccssid,]
nrow(SARCOMA)
# 18
PHENO.ANY_SN$SARCOMA <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% SARCOMA$ccssid, 0, 1))

table(PHENO.ANY_SN$SARCOMA)
# 0    1 
# 2918   18 

PHENO.ANY_SN$AGE.ANY_SN <- SARCOMA$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, SARCOMA$ccssid)]
PHENO.ANY_SN$ANY_SN_TYPE <- SARCOMA$groupdx3[match(PHENO.ANY_SN$ccssid, SARCOMA$ccssid)]
#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
PHENO.ANY_SN$CACO <- PHENO.ANY_SN$SARCOMA

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "smoker_former_or_never_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "PhysicalActivity_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "NOT_RiskyHeavyDrink_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "Not_obese_yn"] <- "Unknown"


#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])

CCSS.exp.ANY_SN <- PHENO.ANY_SN ## CCSS exp


sum(colnames(CCSS.org.ANY_SN) == colnames(CCSS.exp.ANY_SN))
# 54
# Since columns are same, we can simply rbind the dataframes
PHENO.ANY_SN <- rbind.data.frame(CCSS.org.ANY_SN, CCSS.exp.ANY_SN)
rm(list=setdiff(ls(), c("PHENO.ANY_SN")))

table(PHENO.ANY_SN$CACO)
# 0    1 
# 7839  104

## Distribution of overlapping "Unknown" lifestyle factors
source("https://raw.githubusercontent.com/achalneupane/Achal_St_Jude/main/rcodes/attributable_fraction_R_codes/YS_V11_without_diet/get_missing_combination.R")
get_missing_combinations(PHENO.ANY_SN[c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "Not_obese_yn")])
###########################################
## Check data in each category/cross tab ##
###########################################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- CROSS_CASES.df[c("SARCOMA", "smoker_former_or_never_yn", "PhysicalActivity_yn",
                                   "NOT_RiskyHeavyDrink_yn", "Not_obese_yn",  "aa_class_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, SARCOMA = "SARCOMA",  
                               smoker_former_or_never_yn = "smoker_former_or_never_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               NOT_RiskyHeavyDrink_yn = "NOT_RiskyHeavyDrink_yn", Not_obese_yn = "Not_obese_yn", aa_class_dose_5.category= "aa_class_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(SARCOMA, list( smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, Not_obese_yn, aa_class_dose_5.category))))

cc.SARCOMA <- as.data.frame(t(CROSS_CASES.df %>%
                                cross_cases(SARCOMA, list( smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, Not_obese_yn, aa_class_dose_5.category))))

rownames(cc.SARCOMA) <- NULL 


##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SARCOMA ~ Sarcoma_Machiela_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + aa_class_dose_5.category +
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + Not_obese_yn +
                EAS + AFR,
              family = binomial,
              data = dat_all)


summary(fit_all)


##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

###############
## Treatment ##
###############

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all

dat_tx$aa_class_dose_5.category = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.026
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_plp.prs$Sarcoma_Machiela_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# -0.045
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$smoker_former_or_never_yn =
dat_lifestyle$PhysicalActivity_yn =
dat_lifestyle$NOT_RiskyHeavyDrink_yn =
dat_lifestyle$Not_obese_yn = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# 0.494
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$aa_class_dose_5.category = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Sarcoma_Machiela_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn =
dat_tx.plp.prs.lifestyle$PhysicalActivity_yn =
dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn =
dat_tx.plp.prs.lifestyle$Not_obese_yn = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
# 0.491
SARCOMA.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
SARCOMA.res

counts <- cbind.data.frame(cc.SN, cc.SMN, cc.NMSC, cc.BREAST, cc.THYROID, cc.MENINGIOMA, cc.THYROID)
# View(counts)

all.res <- cbind.data.frame(SN=SN.res, SMN=SMN.res, NMSC=NMSC.res, BREAST=BREAST.res, THYROID=THYROID.res, MENINGIOMA=MENINGIOMA.res, SARCOMA=SARCOMA.res)
# View(all.res)
View(t(all.res))