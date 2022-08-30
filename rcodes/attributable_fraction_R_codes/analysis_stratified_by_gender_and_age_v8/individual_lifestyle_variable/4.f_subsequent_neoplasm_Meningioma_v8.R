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
###############
## Meningioma
###############
MENINGIOMA <- subneo[grepl("meningioma", subneo$diag, ignore.case = T),]
MENINGIOMA <- setDT(MENINGIOMA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(MENINGIOMA)
# 149
# Removing samples with SNs within 5 years of childhood cancer
MENINGIOMA <- MENINGIOMA[!MENINGIOMA$sjlid %in% subneo.within5$sjlid,]
nrow(MENINGIOMA)
# 149
PHENO.ANY_SN$MENINGIOMA <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% MENINGIOMA$sjlid, 0, 1))
table(MENINGIOMA$diaggrp)


#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% MENINGIOMA$sjlid, 0, 1))

## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- MENINGIOMA$gradedt[match(ALL.LIFESTYLE$SJLIFEID, MENINGIOMA$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- MENINGIOMA$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, MENINGIOMA$sjlid)]

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$agesurvey > ALL.LIFESTYLE$AGE.ANY_SN),
              c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")] <- NA


ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$agebmi > ALL.LIFESTYLE$AGE.ANY_SN), c("Not_obese_yn")] <- NA



##################################
## Imputation of missing values ##
##################################
# wanted items
wanted.cols <- c("smoker_former_or_never_yn", "PhysicalActivity_yn", "Not_obese_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")
ALL.LIFESTYLE.IMPUTE <- ALL.LIFESTYLE[c("SJLIFEID", wanted.cols)]
## Number of items available to get scores
ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score <- rowSums(!is.na(ALL.LIFESTYLE.IMPUTE[wanted.cols]))

# If all 5 items are available; they were simply summed up to get the score 
ALL.LIFESTYLE.IMPUTE$SCORE [ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 5] <-  rowSums(ALL.LIFESTYLE.IMPUTE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 5, wanted.cols])

## For missing items, per Qi: "For people with missing items, mean of non-missing items were taken and then multiplied by the # of items to get the score"
ALL.LIFESTYLE.IMPUTE$mean_of_non_missing_items <- rowMeans(ALL.LIFESTYLE[wanted.cols], na.rm = T)
ALL.LIFESTYLE.IMPUTE$imputedSCORE <- ALL.LIFESTYLE.IMPUTE$mean_of_non_missing_items * ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score

# replace the score values with imputed score for missing ones; this was applied to rows where items available (non-missing ones) are >= 3 & < 5 
ALL.LIFESTYLE.IMPUTE$SCORE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score >= 3 & ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score < 5] <- ALL.LIFESTYLE.IMPUTE$imputedSCORE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score >= 3 & ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score < 5]

# Define favorable, intermediate and unfavorable----
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS[ALL.LIFESTYLE.IMPUTE$SCORE>= 3] <- "favorable" # favorable [at least 3 of the five health lifestyle factors]
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS[ALL.LIFESTYLE.IMPUTE$SCORE == 2] <- "intermediate" # intermediate [two healthy lifestyle factors] 
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS[ALL.LIFESTYLE.IMPUTE$SCORE <= 1] <- "unfavorable" # unfavorable [no or only one healthy lifestyle factor] 

## Add to PHENO
PHENO.ANY_SN$LIFESTYLE_STATUS <- ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE.IMPUTE$SJLIFEID)]

# relevel factors
PHENO.ANY_SN$LIFESTYLE_STATUS[is.na(PHENO.ANY_SN$LIFESTYLE_STATUS)] <- "Unknown"
PHENO.ANY_SN$LIFESTYLE_STATUS <- factor(PHENO.ANY_SN$LIFESTYLE_STATUS, level = c("favorable", "intermediate", "unfavorable", "Unknown"))


## Repeat for composite without diet
# wanted items
wanted.cols <- c("smoker_former_or_never_yn", "PhysicalActivity_yn", "Not_obese_yn", "NOT_RiskyHeavyDrink_yn")
ALL.LIFESTYLE.IMPUTE <- ALL.LIFESTYLE[c("SJLIFEID", wanted.cols)]
## Number of items available to get scores
ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score <- rowSums(!is.na(ALL.LIFESTYLE.IMPUTE[wanted.cols]))

# If all 4 (after removing HEALTHY_Diet_yn) items are available; they were simply summed up to get the score 
ALL.LIFESTYLE.IMPUTE$SCORE [ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 4] <-  rowSums(ALL.LIFESTYLE.IMPUTE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 4, wanted.cols])

## For missing items, per Qi: "For people with missing items, mean of non-missing items were taken and then multiplied by the # of items to get the score"
ALL.LIFESTYLE.IMPUTE$mean_of_non_missing_items <- rowMeans(ALL.LIFESTYLE[wanted.cols], na.rm = T)
ALL.LIFESTYLE.IMPUTE$imputedSCORE <- ALL.LIFESTYLE.IMPUTE$mean_of_non_missing_items * ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score

# replace the score values with imputed score for missing ones; this was applied to rows where items available (non-missing ones) are >= 3 & < 5 
ALL.LIFESTYLE.IMPUTE$SCORE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score >= 2 & ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score < 4] <- ALL.LIFESTYLE.IMPUTE$imputedSCORE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score >= 2 & ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score < 4]

# Define favorable, intermediate and unfavorable----
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[ALL.LIFESTYLE.IMPUTE$SCORE>= 3] <- "favorable" # favorable [at least 3 of the five health lifestyle factors]
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[ALL.LIFESTYLE.IMPUTE$SCORE == 2] <- "intermediate" # intermediate [two healthy lifestyle factors] 
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[ALL.LIFESTYLE.IMPUTE$SCORE <= 1] <- "unfavorable" # unfavorable [no or only one healthy lifestyle factor] 

## Add to PHENO
PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET <- ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE.IMPUTE$SJLIFEID)]

# relevel factors
PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET[is.na(PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET)] <- "Unknown"
PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET <- factor(PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET, level = c("favorable", "intermediate", "unfavorable", "Unknown"))


############################
## Add lifestyle to Pheno ##
############################
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])
# PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])

# Count missing
PHENO.ANY_SN$missing.lifestyles <- rowSums(is.na(PHENO.ANY_SN[c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")]))
table(PHENO.ANY_SN$missing.lifestyles)
# 0    1    2    3    4    5 
# 2861  593   50   56   75  766 

## Relevel 6 lifestyle variables
PHENO.ANY_SN$smoker_never_yn[is.na(PHENO.ANY_SN$smoker_never_yn)] <- "Unknown"
PHENO.ANY_SN$smoker_never_yn <- factor(PHENO.ANY_SN$smoker_never_yn, level = c(1, 0, "Unknown")) 

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
## Create HEI tertiles ##
#########################
HEI.to.categorize <- c("HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")

## Tertile categories
for(i in 1:length(HEI.to.categorize)){
  TERT = unname(quantile(PHENO.ANY_SN[HEI.to.categorize[i]][PHENO.ANY_SN[HEI.to.categorize[i]] !=0], c(1/3, 2/3, 1), na.rm = T))
  if(sum(duplicated(TERT)) > 0) next
  print (HEI.to.categorize[i])
  print(TERT)
  
  PHENO.ANY_SN$HEI.tmp.tert.category <- as.character(cut(PHENO.ANY_SN[,HEI.to.categorize[i]], breaks = c(0, TERT),
                                                         labels = c("1st", "2nd", "3rd"),
                                                         include.lowest = TRUE))
  
  PHENO.ANY_SN$HEI.tmp.tert.category[is.na(PHENO.ANY_SN$HEI.tmp.tert.category)] <- "Unknown"
  PHENO.ANY_SN$HEI.tmp.tert.category <- factor(PHENO.ANY_SN$HEI.tmp.tert.category, levels = c("3rd", "2nd", "1st", "Unknown"))
  colnames(PHENO.ANY_SN)[colnames(PHENO.ANY_SN) == "HEI.tmp.tert.category"] <- paste0(HEI.to.categorize[i], ".tertile.category")
}

table(PHENO.ANY_SN$HEI2015_TOTAL_SCORE.tertile.category)
# 3rd     2nd     1st Unknown 
# 1160    1159    1161     921
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

CROSS_CASES.df <- CROSS_CASES.df[c("MENINGIOMA", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn",
                                   "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, MENINGIOMA = "MENINGIOMA", smoker_never_yn = "smoker_never_yn", 
                               smoker_former_or_never_yn = "smoker_former_or_never_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               NOT_RiskyHeavyDrink_yn = "NOT_RiskyHeavyDrink_yn", HEALTHY_Diet_yn = "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(MENINGIOMA, list(smoker_never_yn , smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Not_obese_yn))))




#################
## MODEL TESTS ##
#################

###########################
## 1. Qin baseline model ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(MENINGIOMA ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(MENINGIOMA ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)

##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = MENINGIOMA ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
                Meningioma_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + epitxn_dose_5.category + 
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + HEALTHY_Diet_yn + Not_obese_yn +
                EAS + AMR + SAS + AFR,
              family = binomial,
              data = dat_all)

# HEI2015_TOTAL_SCORE.tertile.category

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

dat_tx$maxsegrtdose.category [!grepl("Unknown", dat_tx$maxsegrtdose.category)] =
dat_tx$epitxn_dose_5.category [!grepl("Unknown", dat_tx$epitxn_dose_5.category)] = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all.male = sum(dat_all$pred_all[dat_all$gender == "Male"], na.rm = TRUE) # subset by gender
N_all.female = sum(dat_all$pred_all[dat_all$gender == "Female"], na.rm = TRUE) # subset by gender
## subset by age at diagnosis group
# median(dat_all$AGE_AT_LAST_CONTACT.cs1)
N_all.lt.35 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE) # subset by age 35
N_all.gteq.35 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE) # subset by age 35

N_all.lt.45 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 < 45], na.rm = TRUE) # subset by age 45
N_all.gteq.45 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 45], na.rm = TRUE) # subset by age 45


## Male
N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.male = (N_all.male - N_no_tx) / N_all.male
round(af_by_tx.male,3)

## Female
N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Female"], na.rm = TRUE)
af_by_tx.female = (N_all.female - N_no_tx) / N_all.female
round(af_by_tx.female,3)

## < 35
N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_tx.lt.35 = (N_all.lt.35 - N_no_tx) / N_all.lt.35
round(af_by_tx.lt.35,3)

## >= 35
N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_tx.gteq.35 = (N_all.gteq.35 - N_no_tx) / N_all.gteq.35
round(af_by_tx.gteq.35,3)

## < 45
N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 < 45], na.rm = TRUE)
af_by_tx.lt.45 = (N_all.lt.45 - N_no_tx) / N_all.lt.45
round(af_by_tx.lt.45,3)

## >= 45
N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 45], na.rm = TRUE)
af_by_tx.gteq.45 = (N_all.gteq.45 - N_no_tx) / N_all.gteq.45
round(af_by_tx.gteq.45,3)
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"
dat_plp.prs$dat_prs$Meningioma_PRS.tertile.category = "1st"

## Male
dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs[dat_all$gender == "Male"], na.rm = TRUE)
af_by_plp.prs.male = (N_all.male - N_no_plp.prs) / N_all.male
round(af_by_plp.prs.male,3)

## Female
dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs[dat_all$gender == "Female"], na.rm = TRUE)
af_by_plp.prs.female = (N_all.female - N_no_plp.prs) / N_all.female
round(af_by_plp.prs.female,3)

## < 35
dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_plp.prs.lt.35 = (N_all.lt.35 - N_no_plp.prs) / N_all.lt.35
round(af_by_plp.prs.lt.35,3)


## >= 35
dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_plp.prs.gteq.35 = (N_all.gteq.35 - N_no_plp.prs) / N_all.gteq.35
round(af_by_plp.prs.gteq.35,3)

## < 45
dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs[dat_all$AGE_AT_LAST_CONTACT.cs1 < 45], na.rm = TRUE)
af_by_plp.prs.lt.45 = (N_all.lt.45 - N_no_plp.prs) / N_all.lt.45
round(af_by_plp.prs.lt.45,3)


## >= 45
dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 45], na.rm = TRUE)
af_by_plp.prs.gteq.45 = (N_all.gteq.45 - N_no_plp.prs) / N_all.gteq.45
round(af_by_plp.prs.gteq.45,3)
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

## Male
dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$gender == "Male"], na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category.male = (N_all.male - N_no_favorable_lifestyle.category) / N_all.male
round(af_by_N_no_favorable_lifestyle.category.male,3)

## Female
dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$gender == "Female"], na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category.female = (N_all.female - N_no_favorable_lifestyle.category) / N_all.female
round(af_by_N_no_favorable_lifestyle.category.female,3)

## < 35
dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category.lt.35 = (N_all.lt.35 - N_no_favorable_lifestyle.category) / N_all.lt.35
round(af_by_N_no_favorable_lifestyle.category.lt.35,3)

## >= 35
dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category.gteq.35 = (N_all.gteq.35 - N_no_favorable_lifestyle.category) / N_all.gteq.35
round(af_by_N_no_favorable_lifestyle.category.gteq.35,3)

## < 45
dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 < 45], na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category.lt.45 = (N_all.lt.45 - N_no_favorable_lifestyle.category) / N_all.lt.45
round(af_by_N_no_favorable_lifestyle.category.lt.45,3)

## >= 45
dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 45], na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category.gteq.45 = (N_all.gteq.45 - N_no_favorable_lifestyle.category) / N_all.gteq.45
round(af_by_N_no_favorable_lifestyle.category.gteq.45,3)


setNames(cbind.data.frame(af_by_tx.female, af_by_plp.prs.female, af_by_N_no_favorable_lifestyle.category.female, # female
                          af_by_tx.male, af_by_plp.prs.male, af_by_N_no_favorable_lifestyle.category.male, # male
                          af_by_tx.lt.35, af_by_plp.prs.lt.35, af_by_N_no_favorable_lifestyle.category.lt.35, # less than 35
                          af_by_tx.gteq.35, af_by_plp.prs.gteq.35,af_by_N_no_favorable_lifestyle.category.gteq.35, # 35 or above
                          af_by_tx.lt.45, af_by_plp.prs.lt.45, af_by_N_no_favorable_lifestyle.category.lt.45, # less than 45
                          af_by_tx.gteq.45, af_by_plp.prs.gteq.45,af_by_N_no_favorable_lifestyle.category.gteq.45 # 45 or above
), c("Treatment_female", "Genetics_female", "Lifestyle_female",
     "Treatment_male", "Genetics_male", "Lifestyle_male",
     "Treatment_<35", "Genetics_<35", "Lifestyle_<35",
     "Treatment_>=35", "Genetics_>=35", "Lifestyle_>=35",
     "Treatment_<45", "Genetics_<45", "Lifestyle_<45",
     "Treatment_>=45", "Genetics_>=45", "Lifestyle_>=45"))
