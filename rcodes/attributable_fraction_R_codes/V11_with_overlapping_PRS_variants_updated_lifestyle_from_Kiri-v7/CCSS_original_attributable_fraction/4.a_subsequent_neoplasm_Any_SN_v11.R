load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_org_Genetic_data_P_LP_v11-6.Rdata")

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
## Any SNs ##
#############
## GET SN 18 or older
subneo <- subneo[subneo$AGE.ANY_SN >= 18,]


PHENO.ANY_SN[c("PhysicalActivity_yn_agesurvey", "smoker_former_or_never_yn_agesurvey", "NOT_RiskyHeavyDrink_yn_agesurvey",
                "Not_obese_yn_agesurvey")] <- sapply(PHENO.ANY_SN[c("PhysicalActivity_yn_agesurvey", "smoker_former_or_never_yn_agesurvey", "NOT_RiskyHeavyDrink_yn_agesurvey",
                                                                                                  "Not_obese_yn_agesurvey")], floor)

PHENO.ANY_SN$SURVEY_MIN <- apply(PHENO.ANY_SN[c("PhysicalActivity_yn_agesurvey", "smoker_former_or_never_yn_agesurvey", "NOT_RiskyHeavyDrink_yn_agesurvey",
                                                  "Not_obese_yn_agesurvey")], 1, min)

# Anyone before the first survey, remove them
subneo$survey_First <- PHENO.ANY_SN$SURVEY_MIN[match(subneo$ccssid,PHENO.ANY_SN$ccssid)]
table(subneo$survey_First > subneo$AGE.ANY_SN)

## Anyone not on the first survey age
PHENO.ANY_SN$PhysicalActivity_yn[which(PHENO.ANY_SN$PhysicalActivity_yn_agesurvey != PHENO.ANY_SN$SURVEY_MIN)] <- NA
PHENO.ANY_SN$Current_smoker_yn[which(PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey != PHENO.ANY_SN$SURVEY_MIN)] <- NA
PHENO.ANY_SN$RiskyHeavyDrink_yn[which(PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey != PHENO.ANY_SN$SURVEY_MIN)] <- NA
PHENO.ANY_SN$Obese_yn[which(PHENO.ANY_SN$Not_obese_yn_agesurvey != PHENO.ANY_SN$SURVEY_MIN)] <- NA



subneo <- subneo[!subneo$survey_First > subneo$AGE.ANY_SN,]
#############

# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]


# Removing samples with SN within the 5 years of childhood cancer
ANY_SNs <- ANY_SNs[!ANY_SNs$ccssid %in% subneo.within5$ccssid,]
dim(ANY_SNs)
# 1343

PHENO.ANY_SN$ANY_SN <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% ANY_SNs$ccssid, 0, 1))
table(PHENO.ANY_SN$ANY_SN)
# 0    1 
# 3664 1343 
PHENO.ANY_SN$AGE.ANY_SN <- ANY_SNs$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, ANY_SNs$ccssid)]
#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
PHENO.ANY_SN$CACO <- PHENO.ANY_SN$ANY_SN

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "Current_smoker_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "PhysicalActivity_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "RiskyHeavyDrink_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "Obese_yn"] <- "Unknown"

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


######################################
## Attributable fraction of Any SNs ##
######################################

dat_all = PHENO.ANY_SN
fit_all = glm(formula = ANY_SN ~ Pleiotropy_PRSWEB_PRS.tertile.category + 
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category + 
                Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
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

dat_tx$maxsegrtdose.category [!grepl("Unknown", dat_tx$maxsegrtdose.category)] =
dat_tx$maxabdrtdose.category [!grepl("Unknown", dat_tx$maxabdrtdose.category)] =
dat_tx$maxchestrtdose.category [!grepl("Unknown", dat_tx$maxchestrtdose.category)] =
dat_tx$epitxn_dose_5.category [!grepl("Unknown", dat_tx$epitxn_dose_5.category)] = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.245
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_plp.prs$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# 0.052
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$Current_smoker_yn [!grepl("Unknown", dat_lifestyle$Current_smoker_yn)] =
dat_lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_lifestyle$PhysicalActivity_yn)] =
dat_lifestyle$RiskyHeavyDrink_yn [!grepl("Unknown", dat_lifestyle$RiskyHeavyDrink_yn)] =
dat_lifestyle$Obese_yn [!grepl("Unknown", dat_lifestyle$Obese_yn)] = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# -0.085

#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxsegrtdose.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$maxsegrtdose.category)] =
  dat_tx.plp.prs.lifestyle$maxabdrtdose.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$maxabdrtdose.category)] =
  dat_tx.plp.prs.lifestyle$maxchestrtdose.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$maxchestrtdose.category)] =
  dat_tx.plp.prs.lifestyle$epitxn_dose_5.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$epitxn_dose_5.category)] = "None"

## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$Current_smoker_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Current_smoker_yn)] =
  dat_tx.plp.prs.lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$PhysicalActivity_yn)] =
  dat_tx.plp.prs.lifestyle$RiskyHeavyDrink_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$RiskyHeavyDrink_yn)] =
  dat_tx.plp.prs.lifestyle$Obese_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Obese_yn)] = "1"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
# 0.227
