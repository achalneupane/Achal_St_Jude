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
###################
## Breast cancer ##
###################
BREASTcancer <- subneo[grepl("breast", subneo$diag, ignore.case = T),]
BREASTcancer <- setDT(BREASTcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(BREASTcancer)
# 78
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



########################################
## MODEL TEST for Zhaoming's variants ##
########################################
## SJLIFE (ALL)
# mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + aa_class_dose_any_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN)
mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
# mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + aa_class_dose_any_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


## Prevalence
prevalence.counts <- sum(Zhaoming_vars$Zhaoming_Non.Ref.Counts > 0)
prop.test(prevalence.counts, 4507)


###################################
## MODEL TEST for Qin's variants ##
###################################

###########################
## 1. Qin baseline model ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(BREASTcancer ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)
estimates <- as.data.frame(summary(mod1)$coefficients)[c(1,2)]
estimates$Estimate <- as.numeric(estimates$Estimate)
estimates$RR <- round(exp(estimates$Estimate),2)
## CI
# exp(5.319e-01-(1.96*1.009e-01))
# exp(5.319e-01+(1.96*1.009e-01))
estimates$S.error <- as.numeric(as.character(estimates$`Std. Error`))
CI <-  paste0("(", paste0(round(exp(estimates$Estimate - (1.96*estimates$S.error)),1), " to ", round(exp(estimates$Estimate + (1.96*estimates$S.error)),1)),")")

estimates$RR <- paste(estimates$RR, CI, sep = " ")
estimates
write.table(estimates, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Additional_files/estimates_Qin.txt", sep = "\t", col.names = T, row.names = T, quote = F)

## SJLIFE (EUR)
mod1.EUR <- glm(BREASTcancer ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


###########################
## 2. Qin carriers (all) ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(BREASTcancer ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(BREASTcancer ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)

####################
## 3. HR Pathways ##
####################
## SJLIFE (ALL) 
## Checking with Qi's data
# sum(ANY_SNs$sjlid %in% qi.df.BREAST.filtered$sjlid)
BREASTcancer <- subneo[grepl("breast", subneo$diag, ignore.case = T),]
BREASTcancer <- setDT(BREASTcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
BREASTcancer <- BREASTcancer[BREASTcancer$sjlid %in% qi.df.BREAST.filtered$sjlid,]
dim(BREASTcancer)
# 57
PHENO.ANY_SN$BREASTcancer <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% BREASTcancer$sjlid, 0, 1))

## SJLIFE (ALL) 
mod1 <- glm(BREASTcancer ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)
estimates <- as.data.frame(summary(mod1)$coefficients)[c(1,2)]
estimates$Estimate <- as.numeric(estimates$Estimate)
estimates$OR <- round(exp(estimates$Estimate),2)
## CI
# exp(5.319e-01-(1.96*1.009e-01))
# exp(5.319e-01+(1.96*1.009e-01))
estimates$S.error <- as.numeric(as.character(estimates$`Std. Error`))
CI <-  paste0("(", paste0(round(exp(estimates$Estimate - (1.96*estimates$S.error)),1), " to ", round(exp(estimates$Estimate + (1.96*estimates$S.error)),1)),")")
estimates$OR <- paste(estimates$OR, CI, sep = " ")
estimates

## SJLIFE (EUR)
mod1.EUR <- glm(BREASTcancer ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + maxchestrtdose.category + anthra_jco_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)

