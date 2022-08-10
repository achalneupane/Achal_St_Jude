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
df.1 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])


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
df.2 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

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
df.3 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

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
df.4 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

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
df.5 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])


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
df.6 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])


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
df.7 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

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
df.8 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

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
df.9 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

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
df.1.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])
# Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.2.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_Meta_analysis_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Meta_analysis_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.3.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_PRSWEB_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.4.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_One_directional_PRS.tertile.category 
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_One_directional_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.5.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_Replication_prior_studies_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Replication_prior_studies_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.6.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])


# Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.7.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.8.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_One_directional_Significant_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_One_directional_Significant_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.9.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

paste0("df", 1:9,".WO")
df.ALL <- rbind.data.frame(df.1,df.2,df.3,df.4,df.5,df.6,df.7,df.8,df.9)
df.ALL.0 <- df.ALL[!grepl("H.C.Clin|All.P.LP",rownames(df.ALL)),]
df.ALL.0
colnames(df.ALL.0) <- c("OR_Pleotropy","P_Pleotropy" )

df.ALL.1 <- df.ALL[grepl("H.C.Clin",rownames(df.ALL)),]
cut <- rep(1:(nrow(df.ALL.1)/1), each = 1)
df.ALL.1 <- df.ALL.1[sapply(split(1:nrow(df.ALL.1), cut), c, NA), ]
colnames(df.ALL.1) <- c("OR_H.C_P_LP","P_H.C_P_LP" )

df.ALL.2 <- df.ALL[grepl("All.P.LP",rownames(df.ALL)),]
cut <- rep(1:(nrow(df.ALL.2)/1), each = 1)
df.ALL.2 <- df.ALL.2[sapply(split(1:nrow(df.ALL.2), cut), c, NA), ]
colnames(df.ALL.2) <- c("OR_ALL_P_LP","P_ALL_P_LP" )

df.ALL.WO <- rbind.data.frame(df.1.WO,df.2.WO,df.3.WO,df.4.WO,df.5.WO,df.6.WO,df.7.WO,df.8.WO,df.9.WO)
colnames(df.ALL.WO) <- c("OR_Pleotropy_WO_P_LP", "P_Pleotropy_WO_P_LP")

df.final.prs.check <- cbind.data.frame(df.ALL.0, df.ALL.1, df.ALL.2, df.ALL.WO)

#######################################
## Changing the P/LP without MetaSVM ##
#######################################
dat_all = PHENO.ANY_SN
# Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
#   Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
#   Pleiotropy_Meta_analysis_PRS.tertile.category +
#   Pleiotropy_PRSWEB_PRS.tertile.category +
#   Pleiotropy_One_directional_PRS.tertile.category +

# Pleiotropy_Bi_directional_Increasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.1 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])


# Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.2 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

# Pleiotropy_Meta_analysis_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Meta_analysis_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.3 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

# Pleiotropy_PRSWEB_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.4 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

# Pleiotropy_One_directional_PRS.tertile.category 
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_One_directional_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.5 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])


# Pleiotropy_Replication_prior_studies_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Replication_prior_studies_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.6 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])


# Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.7 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

# Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.8 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

# Pleiotropy_One_directional_Significant_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts +
                All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts +
                Pleiotropy_One_directional_Significant_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.9 <- cbind.data.frame(OR=exp(fit_all$coefficients[4:7]),P=(summary(fit_all))$coefficients[4:7,4])

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
df.1.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])
# Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.2.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_Meta_analysis_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Meta_analysis_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.3.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_PRSWEB_PRS.tertile.category  
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.4.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_One_directional_PRS.tertile.category 
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_One_directional_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.5.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_Replication_prior_studies_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Replication_prior_studies_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.6.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])


# Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Increasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.7.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_Bi_directional_Decreasing_Sig_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.8.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

# Pleiotropy_One_directional_Significant_PRS.tertile.category
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SMN ~ Zhaoming_carriers + Qin_carriers +
                Pleiotropy_One_directional_Significant_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
              data = dat_all)
summary(fit_all)
df.9.WO <- cbind.data.frame(OR=exp(fit_all$coefficients[4:5]),P=(summary(fit_all))$coefficients[4:5,4])

paste0("df", 1:9,".WO")
df.ALL <- rbind.data.frame(df.1,df.2,df.3,df.4,df.5,df.6,df.7,df.8,df.9)
df.ALL.0 <- df.ALL[!grepl("H.C.Clin|All.P.LP",rownames(df.ALL)),]
df.ALL.0
colnames(df.ALL.0) <- c("OR_Pleotropy","P_Pleotropy" )

df.ALL.1 <- df.ALL[grepl("H.C.Clin",rownames(df.ALL)),]
cut <- rep(1:(nrow(df.ALL.1)/1), each = 1)
df.ALL.1 <- df.ALL.1[sapply(split(1:nrow(df.ALL.1), cut), c, NA), ]
colnames(df.ALL.1) <- c("OR_H.C_P_LP","P_H.C_P_LP" )

df.ALL.2 <- df.ALL[grepl("All.P.LP",rownames(df.ALL)),]
cut <- rep(1:(nrow(df.ALL.2)/1), each = 1)
df.ALL.2 <- df.ALL.2[sapply(split(1:nrow(df.ALL.2), cut), c, NA), ]
colnames(df.ALL.2) <- c("OR_ALL_P_LP","P_ALL_P_LP" )

df.ALL.WO <- rbind.data.frame(df.1.WO,df.2.WO,df.3.WO,df.4.WO,df.5.WO,df.6.WO,df.7.WO,df.8.WO,df.9.WO)
colnames(df.ALL.WO) <- c("OR_Pleotropy_WO_P_LP", "P_Pleotropy_WO_P_LP")

df.final.prs.check <- cbind.data.frame(df.ALL.0, df.ALL.1, df.ALL.2, df.ALL.WO)


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

