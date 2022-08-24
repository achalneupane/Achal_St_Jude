#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA")
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
##############
## Any SMNs ##
##############
# This will include any SNs excluding NMSCs
SMNs <- subneo[!grepl("basal cell|squamous cell", subneo$diag, ignore.case = T),]
SMNs <- setDT(SMNs)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(SMNs)
# 476
table(SMNs$diaggrp)

# Removing samples with SNs within 5 years of childhood cancer
SMNs <- SMNs[!SMNs$sjlid %in% subneo.within5$sjlid,]
nrow(SMNs)
# 454
PHENO.ANY_SN$SMN <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% SMNs$sjlid, 0, 1))

#############################
## Add Lifestyle variables ##
#############################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA")
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% SMNs$sjlid, 0, 1))


## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- SMNs$gradedt[match(ALL.LIFESTYLE$SJLIFEID, SMNs$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- SMNs$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, SMNs$sjlid)]

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "Not_Obese_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")] <- NA
# ALL.LIFESTYLE[c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "Not_Obese_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")]


CASES.ALL.LIFESTYLE <- ALL.LIFESTYLE[ALL.LIFESTYLE$CACO == 1,c("smoker_former_or_never_yn", "PhysicalActivity_yn", "Not_Obese_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")]

CASES.ALL.LIFESTYLE$Missing.VAlues <- rowSums(is.na(CASES.ALL.LIFESTYLE))
table(CASES.ALL.LIFESTYLE$Missing.VAlues)








#########################
## Extract Ethnicities ##
#########################
PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]
PHENO.ANY_SN.AFR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'AFR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]


#################
## MODEL TESTS ##
#################

#######################################
## 1. Qin baseline model for ANY SMNs##
#######################################
## SJLIFE (ALL) 
mod1 <- glm(SMN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(SMN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)

######################################
## Attributable fraction of Any SNs ##
######################################

dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category + favorable_lifestyle.category, family = binomial,
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

dat_tx$maxsegrtdose.category = dat_tx$maxabdrtdose.category =  dat_tx$maxchestrtdose.category = dat_tx$epitxn_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.513

##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"; dat_plp.prs$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# 0.2

###############
## Lifestyle ##
###############
dat_lifestyle = dat_all
dat_lifestyle$favorable_lifestyle.category <- "favorable"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# -0.342




