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

# #############################
# ## Add Lifestyle variables ##
# #############################
# ## For each samples, get habits immediately after 18 years of age in agesurvey
# 
# # adlthabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.sas7bdat")
# adlthabits <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.txt", header = T, sep ="\t")
# head(adlthabits)
# # remove duplicated rows
# adlthabits <- distinct(adlthabits)
# ## Get DOB
# adlthabits$DOB <- PHENO.ANY_SN$dob [match(adlthabits$SJLIFEID, PHENO.ANY_SN$sjlid)]
# adlthabits <- adlthabits[!is.na(adlthabits$DOB),]
# # change the format of dates YYYY-MM-DD
# adlthabits$datecomp <- paste(sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 3), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 1), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 2), sep = "-")
# adlthabits$agesurvey <- time_length(interval(as.Date(adlthabits$DOB), as.Date(adlthabits$datecomp)), "years")
# # adlthabits$agesurvey <- floor(adlthabits$agesurvey_diff_of_dob_datecomp)
# 
# samples.sjlife <- unique(adlthabits$SJLIFEID)
# 
# lifestyle <- {}
# for (i in 1:length(samples.sjlife)){
#   print(paste0("Doing ", i))
#   dat <- adlthabits[adlthabits$SJLIFEID == samples.sjlife[i],]
#   if (max(dat$agesurvey) >= 18){
#     print("YES")
#     dat2 <- dat[dat$agesurvey >= 18,]
#     lifestyle.tmp <- dat2[which(dat2$agesurvey == min(dat2$agesurvey)),] # Keep the earliest age after 18 years
#   }
#   lifestyle <- rbind.data.frame(lifestyle, lifestyle.tmp)
# }
# 
# sum(duplicated(lifestyle$SJLIFEID))
# lifestyle$SJLIFEID[duplicated(lifestyle$SJLIFEID)]
# ## Remove duplicate row
# lifestyle <- lifestyle[!duplicated(lifestyle$SJLIFEID),]
# 
# ## Add all samples
# # lifestyle <- cbind.data.frame(wgspop[,1:2], lifestyle[match(wgspop$MRN, lifestyle$mrn), ])
# # lifestyle <- lifestyle[-c(3,4)]
# # tt <- lifestyle
# ## Recode categorical variables
# lifestyle$relation[lifestyle$relation == 1] <- "Self"
# lifestyle$relation[lifestyle$relation == 2] <- "Parent"
# lifestyle$relation[lifestyle$relation == 3] <- "Other"
# 
# ## Recode smoker
# lifestyle$smoker[lifestyle$smoker == 1] <- "Past"
# lifestyle$smoker[lifestyle$smoker == 2] <- "Current"
# lifestyle$smoker[lifestyle$smoker == 3] <- "Never"
# lifestyle$smoker_current_yn <- factor(ifelse(lifestyle$smoker != "Current", "N", "Y"))
# lifestyle$smoker_ever_yn <- factor(ifelse(grepl("Never", lifestyle$smoker), "N", "Y"))
# 
# ## Recode 1/2 or 0/1 to Y/N
# lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink|ltpaw|wtlt|vpa10|yoga", 
#                 colnames(lifestyle))][lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink|ltpaw|wtlt|vpa10|yoga", 
#                                                       colnames(lifestyle))] == 1 ] <- "Y"
# 
# lifestyle[grepl("nopa|ltpa|ltpaw|wtlt|vpa10|yoga", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa|ltpaw|wtlt|vpa10|yoga", colnames(lifestyle))] == 2 ] <- "N"
# 
# lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 0 ] <- "N"

#######################
## Adolescent habits ##
#######################
# 
# adolhabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adolhabits.sas7bdat")
# head(adolhabits)

# ###############
# ## Adult BMI ##
# ###############
# 
# adultbmi <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.sas7bdat")
# head(adultbmi)
# 
# lifestyle$agesurvey_floor <- floor(lifestyle$agesurvey)
# 
# ## Add BMI and Nutrition to Lifestyle. Here, I am extracting the the BMI and Nutrition for the same age for each sample
# lifestyle$BMI_KEY <- paste(lifestyle$SJLIFEID, lifestyle$agesurvey_floor, sep = ":")
# 
# length(unique(adultbmi$sjlid))
# # 3640
# adultbmi$BMI_KEY <- paste(adultbmi$sjlid, adultbmi$AGE, sep = ":")
# sum(lifestyle$BMI_KEY %in% adultbmi$BMI_KEY)
# # 2964
# ## samples that did not match by corresponding age
# cc <- lifestyle[!lifestyle$BMI_KEY %in% adultbmi$BMI_KEY,]
# lifestyle <- cbind.data.frame(lifestyle, adultbmi[match(lifestyle$BMI_KEY, adultbmi$BMI_KEY), c("BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")])


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

# ## SJLIFE (EUR)
# mod1.EUR <- glm(ANY_SN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
# summary(mod1.EUR)


######################################
## Attributable fraction of Any SNs ##
######################################

# Pleiotropy_Bi_directional_Increasing_PRS.decile.category +
# Pleiotropy_Bi_directional_Decreasing_PRS.decile.category +
# Pleiotropy_Meta_analysis_PRS.decile.category +
# Pleiotropy_PRSWEB_PRS.decile.category +
# Pleiotropy_One_directional_PRS.decile.category +

# Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
# Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
# Pleiotropy_Meta_analysis_PRS.tertile.category +
# Pleiotropy_PRSWEB_PRS.tertile.category +
# Pleiotropy_One_directional_PRS.tertile.category +



# H.C.Clin.LoF.MetaSVM.Non.Ref.Counts
# H.C.Clin.LoF.Non.Ref.Counts
# H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts



# We also considered bidirectional pleiotropic associations, wherein the same
# allele for a given variant was associated with an increased risk for some
# cancers but a decreased risk for others.

# For any pair of cancers associated with the same variant, the type of
# association falls in one of three categories: (1) SNPs identified in the
# one-directional analysis, where all associations are in the same direction;
# (2) SNPs identified in the bidirectional analysis, where both cancers in the
# pair are associated in the same direction (both risk increasing or both risk
# decreasing), even though at least one other cancer is associated in the
# opposite direction; and (3) SNPs identified in the bidirectional analysis,
# where the pair of cancers are associated in opposite directions (one risk
# increasing and one risk decreasing).


## Predicted prevalence of ANY_SN
dat_all = PHENO.ANY_SN
# Fit the model - this needs to be the final model with all the variables of interest
# fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_carriers + H.C.Clin.LoF.MetaSVM.Non.Ref.Counts + H.C.Clin.LoF.Non.Ref.Counts + 
#                 H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts + 
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 Pleiotropy_Bi_directional_Increasing_PRS.decile.category +
#                 Pleiotropy_Bi_directional_Decreasing_PRS.decile.category +
#                 Pleiotropy_Meta_analysis_PRS.decile.category +
#                 Pleiotropy_PRSWEB_PRS.decile.category +
#                 Pleiotropy_One_directional_PRS.decile.category +
#                 AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
#                 maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
#               data = dat_all)


# fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_carriers + 
#                 H.C.Clin.LoF.MetaSVM.Non.Ref.Counts + 
#                 H.C.Clin.LoF.Non.Ref.Counts + 
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


# fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_carriers + 
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



# fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_carriers + 
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


# fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_carriers + 
#                 H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts + 
#                 All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 Pleiotropy_Bi_directional_Increasing_PRS.tertile.category +
#                 Pleiotropy_Bi_directional_Decreasing_PRS.tertile.category +
#                 Pleiotropy_One_directional_PRS.tertile.category +
#                 AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
#                 maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
#               data = dat_all)


fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_carriers + 
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

## maxsegrtdose.category
dat_tx = dat_all
dat_tx$maxsegrtdose.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.2024548

## maxabdrtdose.category
dat_tx = dat_all
dat_tx$maxabdrtdose.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.07492755


## maxchestrtdose.category
dat_tx = dat_all
dat_tx$maxchestrtdose.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.1635751

## epitxn_dose_5.category
dat_tx = dat_all
dat_tx$epitxn_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
print(af_by_tx)
# 0.04908683

## P/LP Zhaoming
dat_plp = dat_all
dat_plp$Zhaoming_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp_Zhaoming = (N_all - N_no_plp) / N_all
print(af_by_plp_Zhaoming)
# 0.01999075

## P/LP Qin
dat_plp = dat_all
dat_plp$Qin_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp_Qin = (N_all - N_no_plp) / N_all
print(af_by_plp_Qin)
# -0.0006263049

## H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts
dat_plp = dat_all
dat_plp$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts  = 0

dat_all$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = sum(dat_all$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts, na.rm = TRUE)
af_by_N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = (N_all - N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts) / N_all
print(af_by_N_no_pred_H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts)
# -2.216773

## All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts
dat_plp = dat_all
dat_plp$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts  = 0

dat_all$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = predict(fit_all, newdata = dat_plp, type = "response")
N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = sum(dat_all$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts, na.rm = TRUE)
af_by_N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = (N_all - N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts) / N_all
print(af_by_N_no_pred_All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts)
# -2.424273

#########
## PRS ##
#########
dat_prs = dat_all
dat_prs$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

dat_all$pred_no_Pleiotropy_PRSWEB_PRS.tertile.category = predict(fit_all, newdata = dat_prs, type = "response")
N_no_Pleiotropy_PRSWEB_PRS.tertile.category = sum(dat_all$pred_no_Pleiotropy_PRSWEB_PRS.tertile.category, na.rm = TRUE)
af_by_N_no_Pleiotropy_PRSWEB_PRS.tertile.category = (N_all - N_no_Pleiotropy_PRSWEB_PRS.tertile.category) / N_all
print(af_by_N_no_Pleiotropy_PRSWEB_PRS.tertile.category)

