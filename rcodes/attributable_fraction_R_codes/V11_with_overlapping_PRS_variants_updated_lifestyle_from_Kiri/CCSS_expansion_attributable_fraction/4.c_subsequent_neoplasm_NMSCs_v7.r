load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_exp_Genetic_data_P_LP.Rdata")

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
# 269

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$ccssid))
# 8




###########
## NMSCs ##
###########
# This will include basal cell, squamous cell and melanoma
NMSCs <- subneo[grepl("skin", subneo$groupdx3, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(NMSCs)
# 11

# Removing samples with SNs within 5 years of childhood cancer
NMSCs <- NMSCs[!NMSCs$ccssid %in% subneo.within5$ccssid,]
nrow(NMSCs)
# 11
PHENO.ANY_SN$NMSC <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% NMSCs$ccssid, 0, 1))



















#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/final.5.Q_header2_SJLIFE_only", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AMR", "SAS", "AFR")])


PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]
PHENO.ANY_SN.AFR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'AFR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]


###########################################
## Check data in each category/cross tab ##
###########################################
library(expss)

# Getting counts for non-missing data only; 26 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$AMR),]

CROSS_CASES.df <- CROSS_CASES.df[c("NMSC", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn",
                                   "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, NMSC = "NMSC", smoker_never_yn = "smoker_never_yn", 
                               smoker_former_or_never_yn = "smoker_former_or_never_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               NOT_RiskyHeavyDrink_yn = "NOT_RiskyHeavyDrink_yn", HEALTHY_Diet_yn = "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(NMSC, list(smoker_never_yn , smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Not_obese_yn))))



###################################
## MODEL TEST for Qin's variants ##
###################################

###########################
## 1. Qin baseline model ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(NMSC ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(NMSC ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = NMSC ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
                BASALcell_PRS.tertile.category + SQUAMOUScell_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn  + HEALTHY_Diet_yn + Not_obese_yn +
                EAS + AMR + SAS + AFR,
              family = binomial,
              data = dat_all)

summary(fit_all)

# ## With HEI Score
# dat_all = PHENO.ANY_SN
# fit_all = glm(formula = NMSC ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
#                 BASALcell_PRS.tertile.category + SQUAMOUScell_PRS.tertile.category +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
#                 smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + HEI2015_TOTAL_SCORE.tertile.category + Not_obese_yn +
#                 EAS + AMR + SAS + AFR,
#               family = binomial,
#               data = dat_all)
# 
# summary(fit_all)


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
dat_tx$maxpelvisrtdose.category [!grepl("Unknown", dat_tx$maxpelvisrtdose.category)] = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.189
# 0.186 (Without diet)
# 0.187 [With HEI 2015 (instead of Diet) and other 4 lifestyle factors]
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"
dat_plp.prs$SQUAMOUScell_PRS.tertile.category = dat_plp.prs$BASALcell_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# 0.22
# 0.22 (Without diet)
# 0.229 [With HEI 2015 (instead of Diet) and other 4 lifestyle factors]
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$smoker_former_or_never_yn [!grepl("Unknown", dat_lifestyle$smoker_former_or_never_yn)] =
dat_lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_lifestyle$PhysicalActivity_yn)] =
dat_lifestyle$NOT_RiskyHeavyDrink_yn [!grepl("Unknown", dat_lifestyle$NOT_RiskyHeavyDrink_yn)] =
dat_lifestyle$HEALTHY_Diet_yn [!grepl("Unknown", dat_lifestyle$HEALTHY_Diet_yn)] =
dat_lifestyle$Not_obese_yn [!grepl("Unknown", dat_lifestyle$Not_obese_yn)] = "1"

# ## HEI
# dat_lifestyle$smoker_former_or_never_yn [!grepl("Unknown", dat_lifestyle$smoker_former_or_never_yn)] =
# dat_lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_lifestyle$PhysicalActivity_yn)] =
# dat_lifestyle$NOT_RiskyHeavyDrink_yn [!grepl("Unknown", dat_lifestyle$NOT_RiskyHeavyDrink_yn)] =
# dat_lifestyle$Not_obese_yn [!grepl("Unknown", dat_lifestyle$Not_obese_yn)] = "1"
# dat_lifestyle$HEI2015_TOTAL_SCORE.tertile.category [!grepl("Unknown", dat_lifestyle$HEI2015_TOTAL_SCORE.tertile.category)] = "3rd"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# -0.311
# -0.16
# -0.233 [With HEI 2015 (instead of Diet) and other 4 lifestyle factors]


#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxsegrtdose.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$maxsegrtdose.category)] =
  dat_tx.plp.prs.lifestyle$maxabdrtdose.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$maxabdrtdose.category)] =
  dat_tx.plp.prs.lifestyle$maxpelvisrtdose.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$maxpelvisrtdose.category)] = "None"


## Nullify Genetics
dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$SQUAMOUScell_PRS.tertile.category = dat_tx.plp.prs.lifestyle$BASALcell_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn)] =
  dat_tx.plp.prs.lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$PhysicalActivity_yn)] =
  dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn)] =
  dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn)] =
  dat_tx.plp.prs.lifestyle$Not_obese_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Not_obese_yn)] = "1"


# ## HEI
# dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn)] =
# dat_tx.plp.prs.lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$PhysicalActivity_yn)] =
# dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn)] =
# dat_tx.plp.prs.lifestyle$Not_obese_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Not_obese_yn)] = "1"
# dat_tx.plp.prs.lifestyle$HEI2015_TOTAL_SCORE.tertile.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$HEI2015_TOTAL_SCORE.tertile.category)] = "3rd"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
# 0.218 With diet
# 0.282 Without diet
