#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
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
#################### 
## Thyroid cancer ##
####################
table(subneo$diag)


THYROIDcancer <- subneo[grepl("thyroid", subneo$diag, ignore.case = T),]
THYROIDcancer <- setDT(THYROIDcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(THYROIDcancer)
# Removing samples with SNs within 5 years of childhood cancer
THYROIDcancer <- THYROIDcancer[!THYROIDcancer$sjlid %in% subneo.within5$sjlid,]
nrow(THYROIDcancer)
# 86
PHENO.ANY_SN$THYROIDcancer <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% THYROIDcancer$sjlid, 0, 1))
table(THYROIDcancer$diaggrp)

#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% THYROIDcancer$sjlid, 0, 1))

## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- THYROIDcancer$gradedt[match(ALL.LIFESTYLE$SJLIFEID, THYROIDcancer$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- THYROIDcancer$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, THYROIDcancer$sjlid)]

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("smoker_former_or_never_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("PhysicalActivity_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("NOT_RiskyHeavyDrink_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEALTHY_Diet_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEALTHY_Diet_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$Not_obese_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("Not_obese_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEI2015_TOTAL_SCORE_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEI2015_TOTAL_SCORE")] <- NA

#############################
## Addd lifestyle to Pheno ##
#############################
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2015_TOTAL_SCORE", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])
# PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])

# Count missing
PHENO.ANY_SN$missing.lifestyles <- rowSums(is.na(PHENO.ANY_SN[c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")]))
table(PHENO.ANY_SN$missing.lifestyles)
# 0    1    2    3    4    5 
# 4333    2   27    9    4   26

## Relevel 6 lifestyle variables
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

## Distribution of overlapping "Unknown" lifestyle factors
source("https://raw.githubusercontent.com/achalneupane/Achal_St_Jude/main/rcodes/attributable_fraction_R_codes/YS_V11_without_diet/get_missing_combination.R")
get_missing_combinations(PHENO.ANY_SN[c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])


#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


###########################################
## Check data in each category/cross tab ##
###########################################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- CROSS_CASES.df[c("THYROIDcancer", "smoker_former_or_never_yn", "PhysicalActivity_yn",
                                   "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, THYROIDcancer = "THYROIDcancer", 
                               smoker_former_or_never_yn = "smoker_former_or_never_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               NOT_RiskyHeavyDrink_yn = "NOT_RiskyHeavyDrink_yn", HEALTHY_Diet_yn = "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(THYROIDcancer, list(smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Not_obese_yn))))

cc.THYROID <- as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(THYROIDcancer, list(smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Not_obese_yn))))

rownames(cc.THYROID) <- NULL 
# View(cc)

##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = THYROIDcancer ~ Thyroid_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxneckrtdose.category + epitxn_dose_5.category + 
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + HEALTHY_Diet_yn + Not_obese_yn +
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


dat_tx$maxneckrtdose.category =
dat_tx$epitxn_dose_5.category = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.5552
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"
dat_plp.prs$Thyroid_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# 0.345
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$smoker_former_or_never_yn =
dat_lifestyle$PhysicalActivity_yn =
dat_lifestyle$NOT_RiskyHeavyDrink_yn =
dat_lifestyle$HEALTHY_Diet_yn =
dat_lifestyle$Not_obese_yn = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# 0.936
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxneckrtdose.category =
  dat_tx.plp.prs.lifestyle$epitxn_dose_5.category = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Thyroid_PRS.tertile.category = "1st"


## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn =
  dat_tx.plp.prs.lifestyle$PhysicalActivity_yn =
  dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn =
  dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn =
  dat_tx.plp.prs.lifestyle$Not_obese_yn = "1"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
# 0.993
THYROID.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
THYROID.res
# 0.555 0.341 0.703 0.971

# PAFs can sum to more than 1 because some individuals with more than one risk
# factor can have disease prevented in more than one way, and the prevented
# cases of these individuals could be counted more than once. (source:  DOI: 10.1016/j.amepre.2003.12.007)