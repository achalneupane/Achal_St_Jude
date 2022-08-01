#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores.RDATA")
#########################
## Subsequent Neoplasm ##
#########################

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


load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores.RDATA")


subneo <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/subneo.sas7bdat")
head(subneo)
table(subneo$diaggrp)
dim(subneo)
# 1731 9

subneo <- subneo[subneo$sjlid %in% PHENO.ANY_SN$sjlid ,]
dim(subneo)
# 1717
# add diagnosis date 
subneo$diagdt <-  PHENO.ANY_SN$diagdt [match(subneo$sjlid , PHENO.ANY_SN$sjlid)]
subneo$agedx <-  PHENO.ANY_SN$agedx [match(subneo$sjlid , PHENO.ANY_SN$sjlid)]
# add DOB
subneo$DOB <- PHENO.ANY_SN$dob[match(subneo$sjlid, PHENO.ANY_SN$sjlid)]

subneo$AGE.ANY_SN <- time_length(interval(as.Date(subneo$DOB), as.Date(subneo$gradedt)), "years")

## These two dates should be the (almost) same
subneo$AGE.ANY_SN.after.childhood.cancer <- time_length(interval(as.Date(subneo$diagdt), as.Date(subneo$gradedt)), "years")
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$sjlid))
# 612

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 22
#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


# Removing samples with with SN within the 5 years of childhood cancer
ANY_SNs <- ANY_SNs[!ANY_SNs$sjlid %in% subneo.within5$sjlid,]
dim(ANY_SNs)
# 605

PHENO.ANY_SN$ANY_SN <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))

#############################
## Add Lifestyle variables ##
#############################








#########################
## Extract Ethnicities ##
#########################
PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]
PHENO.ANY_SN.AFR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'AFR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]


#################
## MODEL TESTS ##
#################

###########################
## 1. Qin baseline model ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(ANY_SN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(ANY_SN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


######################################
## Attributable fraction of Any SNs ##
######################################

## Predicted prevalence of ANY_SN
dat_all = PHENO.ANY_SN
# Fit the model - this needs to be the final model with all the variables of interest
fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_carriers.HR.pathways + Qin_carriers.FA.pathways + Qin_carriers.MMR.pathways + 
                Qin_carriers.BER.pathways + Qin_carriers.NER.pathways + Qin_carriers.NHEJ.pathways + AGE_AT_LAST_CONTACT.cs1 +
                AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category + Pleiotropy_Bi_directional_Increasing_Sig_PRS + 
                Pleiotropy_Bi_directional_Decreasing_Sig_PRS, family = binomial,
              data = dat_all)

# Get predicted values
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all
dat_tx$maxsegrtdose.category = dat_tx$maxabdrtdose.category = dat_tx$maxchestrtdose.category = dat_tx$epitxn_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation
# First get the "predicted" number of SNs
# Based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# paste0(round(af_by_tx *100, 2), "%")

## P/LP
dat_plp = dat_all
dat_plp$Zhaoming_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp = (N_all - N_no_plp) / N_all
print(af_by_plp)

# # H.C.Clin.LoF.Non.Ref.Counts
# dat_plp = dat_all
# dat_plp$H.C.Clin.LoF.Non.Ref.Counts = 0
# 
# dat_all$pred_no_H.C.Clin.LoF.Non.Ref.Counts = predict(fit_all, newdata = dat_plp, type = "response")
# N_no_pred_no_H.C.Clin.LoF.Non.Ref.Counts = sum(dat_all$pred_no_H.C.Clin.LoF.Non.Ref.Counts, na.rm = TRUE)
# af_by_H.C.Clin.LoF.Non.Ref.Counts = (N_all - N_no_pred_no_H.C.Clin.LoF.Non.Ref.Counts) / N_all
# print(af_by_H.C.Clin.LoF.Non.Ref.Counts)


# Pleiotropy_Bi_directional_Increasing_Sig_PRS
dat_plp = dat_all
dat_plp$Pleiotropy_Bi_directional_Increasing_Sig_PRS = 0

dat_all$pred_no_Pleiotropy_Bi_directional_Increasing_Sig_PRS = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_no_Pleiotropy_Bi_directional_Increasing_Sig_PRS = sum(dat_all$pred_no_Pleiotropy_Bi_directional_Increasing_Sig_PRS, na.rm = TRUE)
af_by_Pleiotropy_Bi_directional_Increasing_Sig_PRS = (N_all - N_no_pred_no_Pleiotropy_Bi_directional_Increasing_Sig_PRS) / N_all
print(af_by_Pleiotropy_Bi_directional_Increasing_Sig_PRS)

# Pleiotropy_Bi_directional_Decreasing_Sig_PRS
dat_plp = dat_all
dat_plp$Pleiotropy_Bi_directional_Decreasing_Sig_PRS = 0

dat_all$pred_no_Pleiotropy_Bi_directional_Decreasing_Sig_PRS = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_no_Pleiotropy_Bi_directional_Decreasing_Sig_PRS = sum(dat_all$pred_no_Pleiotropy_Bi_directional_Decreasing_Sig_PRS, na.rm = TRUE)
af_by_Pleiotropy_Bi_directional_Decreasing_Sig_PRS = (N_all - N_no_pred_no_Pleiotropy_Bi_directional_Decreasing_Sig_PRS) / N_all
print(af_by_Pleiotropy_Bi_directional_Decreasing_Sig_PRS)
