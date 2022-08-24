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

# write.table(PHENO.ANY_SN$sjlid, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/4401_attributable_fraction_ids.txt", col.names = F, row.names = F, quote = F)
#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% ANY_SNs$sjlid, 0, 1))

## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- ANY_SNs$gradedt[match(ALL.LIFESTYLE$SJLIFEID, ANY_SNs$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- ANY_SNs$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, ANY_SNs$sjlid)]

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$agesurvey > ALL.LIFESTYLE$AGE.ANY_SN),
              c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")] <- NA

ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$agebmi > ALL.LIFESTYLE$AGE.ANY_SN), c("Not_obese_yn")] <- NA


## Imputation of missing
# ALL.LIFESTYLE[c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "Not_Obese_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")]
# CASES.ALL.LIFESTYLE <- ALL.LIFESTYLE[ALL.LIFESTYLE$CACO == 1,c("smoker_former_or_never_yn", "PhysicalActivity_yn", "Not_obese_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")]
# CASES.ALL.LIFESTYLE$Missing.Values <- rowSums(!is.na(CASES.ALL.LIFESTYLE))
# ## If 50 percent or less missing items (i.e, <= 2 missing items)
# 
# 
# CASES.ALL.LIFESTYLE$goodSum[CASES.ALL.LIFESTYLE$Missing.Values <= 2 & CASES.ALL.LIFESTYLE$Missing.Values >= 1] <- rowSums(CASES.ALL.LIFESTYLE[1:5], na.rm = T)[CASES.ALL.LIFESTYLE$Missing.Values <= 2 & CASES.ALL.LIFESTYLE$Missing.Values >= 1]
# CASES.ALL.LIFESTYLE$goodSumMean <- CASES.ALL.LIFESTYLE$goodSum 
# CASES.ALL.LIFESTYLE$SCORE[CASES.ALL.LIFESTYLE$Missing.Values <= 2 & CASES.ALL.LIFESTYLE$Missing.Values >= 1] <- CASES.ALL.LIFESTYLE$goodSum[CASES.ALL.LIFESTYLE$Missing.Values <= 2 & CASES.ALL.LIFESTYLE$Missing.Values >= 1]

# # Number of non-missing items
# CASES.ALL.LIFESTYLE$Non.Missing.Values <- 5-CASES.ALL.LIFESTYLE$Missing.Values 



## Addd lifestyle to Pheno
# PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])

# Count missing
PHENO.ANY_SN$missing.lifestyles <- rowSums(is.na(PHENO.ANY_SN[c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")]))
table(PHENO.ANY_SN$missing.lifestyles)
# 0    1    2    3    4    5 
# 2670  567   46   54   78  986 

## Relevel 5 lifestyle variables
# PHENO.ANY_SN$smoker_never_yn[is.na(PHENO.ANY_SN$smoker_never_yn)] <- "Unknown"; PHENO.ANY_SN$smoker_never_yn <- factor(PHENO.ANY_SN$smoker_never_yn, level = c(0, "Unknown", 1)) 
PHENO.ANY_SN$smoker_former_or_never_yn[is.na(PHENO.ANY_SN$smoker_former_or_never_yn)] <- "Unknown"
PHENO.ANY_SN$smoker_former_or_never_yn <- factor(PHENO.ANY_SN$smoker_former_or_never_yn, level = c(1, 0, "Unknown")) 

PHENO.ANY_SN$PhysicalActivity_yn[is.na(PHENO.ANY_SN$PhysicalActivity_yn)] <- "Unknown"
PHENO.ANY_SN$PhysicalActivity_yn <- factor(PHENO.ANY_SN$PhysicalActivity_yn, level = c(1, 0, "Unknown")) 

PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn[is.na(PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn)] <- "Unknown"
PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn <- factor(PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn, level = c(1, 0, "Unknown")) 

PHENO.ANY_SN$HEALTHY_Diet_yn[is.na(PHENO.ANY_SN$HEALTHY_Diet_yn)] <- "Unknown"
PHENO.ANY_SN$HEALTHY_Diet_yn <- factor(PHENO.ANY_SN$HEALTHY_Diet_yn, level = c(1, 0, "Unknown")) 

PHENO.ANY_SN$Not_obese_yn[is.na(PHENO.ANY_SN$Not_obese_yn)] <- "Unknown";
PHENO.ANY_SN$Not_obese_yn <- factor(PHENO.ANY_SN$Not_obese_yn, level = c(1, 0, "Unknown")) 


#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/final.5.Q_header2_SJLIFE_only", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AMR", "SAS", "AFR")])


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

dat_all = PHENO.ANY_SN
fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
                Pleiotropy_PRSWEB_PRS.tertile.category + 
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category + 
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + HEALTHY_Diet_yn + Not_obese_yn +
                EAS + AMR + SAS + AFR,
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

dat_tx$maxsegrtdose.category = dat_tx$maxabdrtdose.category =  dat_tx$maxchestrtdose.category = dat_tx$epitxn_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.284
# 0.283 (Without diet)
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
# 0.081
# 0.081 (Without diet)
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all
dat_lifestyle$smoker_former_or_never_yn = dat_lifestyle$PhysicalActivity_yn = dat_lifestyle$NOT_RiskyHeavyDrink_yn = dat_lifestyle$HEALTHY_Diet_yn = dat_lifestyle$Not_obese_yn = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# 0.301
# 0.305 (Without diet)


