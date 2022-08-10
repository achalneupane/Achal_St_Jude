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

############################
## Attributable Fractions ##
############################
dat_all = PHENO.ANY_SN

# fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
#                 H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
#                 All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
#                 Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
#                 Pleiotropy_Meta_analysis_PRS.tertile.category +
#                 Pleiotropy_PRSWEB_PRS.tertile.category +
#                 Pleiotropy_One_directional_PRS.tertile.category +
#                 AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
#                 maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
#               data = dat_all)

# dat_all = PHENO.ANY_SN
# fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
#                 H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
#                 All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
#                 Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
#                 Pleiotropy_Meta_analysis_PRS.tertile.category +
#                 Pleiotropy_One_directional_PRS.tertile.category +
#                 AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
#                 maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
#               data = dat_all)
# 
# dat_all = PHENO.ANY_SN
# fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
#                 H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
#                 All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
#                 Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
#                 Pleiotropy_One_directional_PRS.tertile.category +
#                 AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
#                 maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
#               data = dat_all)

###############################
## Fitting One PRS at a time ##
###############################
dat_all = PHENO.ANY_SN
# Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
#   Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
#   Pleiotropy_Meta_analysis_PRS.tertile.category +
#   Pleiotropy_PRSWEB_PRS.tertile.category +
#   Pleiotropy_One_directional_PRS.tertile.category +

# Pleiotropy_Bi_directional_Increasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_Meta_analysis_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Meta_analysis_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_PRSWEB_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)


# Pleiotropy_One_directional_PRS.tertile.category 
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_One_directional_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)


# Pleiotropy_Replication_prior_studies_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Replication_prior_studies_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)



# Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)


# Pleiotropy_One_directional_Significant_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_One_directional_Significant_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)


## Repeat without H.C and P/LP in the Genome

dat_all = PHENO.ANY_SN
# Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
#   Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
#   Pleiotropy_Meta_analysis_PRS.tertile.category +
#   Pleiotropy_PRSWEB_PRS.tertile.category +
#   Pleiotropy_One_directional_PRS.tertile.category +

# Pleiotropy_Bi_directional_Increasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_Meta_analysis_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Meta_analysis_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_PRSWEB_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)


# Pleiotropy_One_directional_PRS.tertile.category 
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_One_directional_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)


# Pleiotropy_Replication_prior_studies_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Replication_prior_studies_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)



# Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

# Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)


# Pleiotropy_One_directional_Significant_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_One_directional_Significant_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)

#############################################################################



dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers + 
                H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts + 
                All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)





summary(fit_all)

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
# 0.4947382

## maxsegrtdose.category
dat_tx = dat_all
dat_tx$maxsegrtdose.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.2138112

## maxabdrtdose.category
dat_tx = dat_all
dat_tx$maxabdrtdose.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.09324085


## maxchestrtdose.category
dat_tx = dat_all
dat_tx$maxchestrtdose.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.1786089

## epitxn_dose_5.category
dat_tx = dat_all
dat_tx$epitxn_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.07783655

## P/LP Zhaoming
dat_plp = dat_all
dat_plp$Zhaoming_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp_Zhaoming = (N_all - N_no_plp) / N_all
print(af_by_plp_Zhaoming)
# 0.02267418

## P/LP Qin
dat_plp = dat_all
dat_plp$Qin_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp_Qin = (N_all - N_no_plp) / N_all
print(af_by_plp_Qin)
# 0.01281991

## H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts
dat_plp = dat_all
dat_plp$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts  = 0

dat_all$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = sum(dat_all$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts, na.rm = TRUE)
af_by_N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = (N_all - N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts) / N_all
print(af_by_N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts)
# -3.298338

## All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts
dat_plp = dat_all
dat_plp$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts  = 0

dat_all$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = sum(dat_all$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts, na.rm = TRUE)
af_by_N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = (N_all - N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts) / N_all
print(af_by_N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts)
# -1.44166

#########
## PRS ##
#########
dat_prs = dat_all
dat_prs$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

dat_all$pred_no_Pleiotropy_PRSWEB_PRS.tertile.category = predict(fit_all, newdata = dat_prs, type = "response")
N_no_Pleiotropy_PRSWEB_PRS.tertile.category = sum(dat_all$pred_no_Pleiotropy_PRSWEB_PRS.tertile.category, na.rm = TRUE)
af_by_N_no_Pleiotropy_PRSWEB_PRS.tertile.category = (N_all - N_no_Pleiotropy_PRSWEB_PRS.tertile.category) / N_all
print(af_by_N_no_Pleiotropy_PRSWEB_PRS.tertile.category)
# -1.168303

