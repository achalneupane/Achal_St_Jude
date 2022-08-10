#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_categories.RDATA")
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
###################
## Breast cancer ##
###################
BREASTcancer <- subneo[grepl("breast", subneo$diag, ignore.case = T),]
BREASTcancer <- setDT(BREASTcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(BREASTcancer)
# 78
table(BREASTcancer$diaggrp)
# Removing samples with SNs within 5 years of childhood cancer
BREASTcancer <- BREASTcancer[!BREASTcancer$sjlid %in% subneo.within5$sjlid,]
nrow(BREASTcancer)
# 76
PHENO.ANY_SN$BREASTcancer <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% BREASTcancer$sjlid, 0, 1))

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
mod1 <- glm(BREASTcancer ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(BREASTcancer ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


############################
## Attributable Fractions ##
############################

# PRS 2015
dat_all = PHENO.ANY_SN
fit_all = glm(formula = BREASTcancer ~ Zhaoming_carriers + Qin_carriers + 
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts + 
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Mavaddat_2015_ER_POS_Breast_PRS.tertile.category +
                Mavaddat_2015_ER_OVERALL_Breast_PRS.tertile.category +
                Mavaddat_2015_ER_NEG_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2+ 
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial,
                data = dat_all)

summary(fit_all)

# PRS 2019
dat_all = PHENO.ANY_SN
fit_all = glm(formula = BREASTcancer ~ Zhaoming_carriers + Qin_carriers + 
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts + 
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Mavaddat_2019_ER_POS_Breast_PRS.tertile.category +
                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                Mavaddat_2019_ER_NEG_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2+ 
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial,
              data = dat_all)

summary(fit_all)

# PRSWEB
dat_all = PHENO.ANY_SN
fit_all = glm(formula = BREASTcancer ~ Zhaoming_carriers + Qin_carriers + 
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts + 
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                MichiganWeb_ER_NEG_Breast_PRS.tertile.category +
                MichiganWeb_ER_OVERALL_Breast_PRS.tertile.category +
                MichiganWeb_ER_POS_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2+ 
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial,
              data = dat_all)

summary(fit_all)

# Get predicted values
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all
dat_tx$maxchestrtdose.category = dat_tx$anthra_jco_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation
# First get the "predicted" number of SNs
# Based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.5470882
# 0.5516615
# 0.5523445

## maxchestrtdose.category
dat_tx = dat_all
dat_tx$maxchestrtdose.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.4945335
# 0.4954146
# 0.5004813

## anthra_jco_dose_5.category
dat_tx = dat_all
dat_tx$anthra_jco_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.100922
# 0.1084345
# 0.09832719

## P/LP Zhaoming
dat_plp = dat_all
dat_plp$Zhaoming_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp_Zhaoming = (N_all - N_no_plp) / N_all
print(af_by_plp_Zhaoming)
# 0.05583012
# 0.05341221
# 0.05524048

## P/LP Qin
dat_plp = dat_all
dat_plp$Qin_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp_Qin = (N_all - N_no_plp) / N_all
print(af_by_plp_Qin)
# 0.0008419519
# 0.0008096127
# -0.001769915

## H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts
dat_plp = dat_all
dat_plp$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts  = 0

dat_all$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = sum(dat_all$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts, na.rm = TRUE)
af_by_N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = (N_all - N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts) / N_all
print(af_by_N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts)
# -24.57465
# -23.21813
# -21.60368

## All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts
dat_plp = dat_all
dat_plp$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts  = 0

dat_all$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = sum(dat_all$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts, na.rm = TRUE)
af_by_N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = (N_all - N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts) / N_all
print(af_by_N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts)
# -4.260216
# -2.897749
# -6.729559

#########
## PRS ##
#########
# Mavaddat_2015_Breast
dat_prs = dat_all
dat_prs$Mavaddat_2015_ER_POS_Breast_PRS.tertile.category = dat_prs$Mavaddat_2015_ER_NEG_Breast_PRS.tertile.category = dat_prs$Mavaddat_2015_ER_OVERALL_Breast_PRS.tertile.category = "1st"

dat_all$pred_no_Mavaddat_2015.tertile.category = predict(fit_all, newdata = dat_prs, type = "response")
N_no_pred_no_Mavaddat_2015.tertile.category = sum(dat_all$pred_no_Mavaddat_2015.tertile.category, na.rm = TRUE)
af_by_N_no_pred_no_Mavaddat_2015.tertile.category = (N_all - N_no_pred_no_Mavaddat_2015.tertile.category) / N_all
print(af_by_N_no_pred_no_Mavaddat_2015.tertile.category)
# -2.472324


# Mavaddat_2019_Breast
dat_prs = dat_all
dat_prs$Mavaddat_2019_ER_POS_Breast_PRS.tertile.category = dat_prs$Mavaddat_2019_ER_NEG_Breast_PRS.tertile.category = dat_prs$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st"

dat_all$pred_no_Mavaddat_2019.tertile.category = predict(fit_all, newdata = dat_prs, type = "response")
N_no_pred_no_Mavaddat_2019.tertile.category = sum(dat_all$pred_no_Mavaddat_2019.tertile.category, na.rm = TRUE)
af_by_N_no_pred_no_Mavaddat_2019.tertile.category = (N_all - N_no_pred_no_Mavaddat_2019.tertile.category) / N_all
print(af_by_N_no_pred_no_Mavaddat_2019.tertile.category)
# -1.016479

# PRSWEB
dat_prs = dat_all
dat_prs$MichiganWeb_ER_NEG_Breast_PRS.tertile.category = dat_prs$MichiganWeb_ER_OVERALL_Breast_PRS.tertile.category = dat_prs$MichiganWeb_ER_POS_Breast_PRS.tertile.category = "1st"

dat_all$pred_no_MichiganWeb.tertile.category = predict(fit_all, newdata = dat_prs, type = "response")
N_no_pred_no_MichiganWeb.tertile.category = sum(dat_all$pred_no_MichiganWeb.tertile.category, na.rm = TRUE)
af_by_N_no_pred_no_MichiganWeb.tertile.category = (N_all - N_no_pred_no_MichiganWeb.tertile.category) / N_all
print(af_by_N_no_pred_no_MichiganWeb.tertile.category)
# -4.377463
