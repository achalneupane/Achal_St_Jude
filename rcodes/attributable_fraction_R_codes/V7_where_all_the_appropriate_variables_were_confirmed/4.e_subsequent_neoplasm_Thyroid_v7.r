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
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$agesurvey > ALL.LIFESTYLE$AGE.ANY_SN),
              c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")] <- NA


ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$agebmi > ALL.LIFESTYLE$AGE.ANY_SN), c("Not_obese_yn")] <- NA



## Addd lifestyle to Pheno
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])
# PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")])

# Count missing
PHENO.ANY_SN$missing.lifestyles <- rowSums(is.na(PHENO.ANY_SN[c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Not_obese_yn")]))
table(PHENO.ANY_SN$missing.lifestyles)
# 0    1    2    3    4    5 
# 2880  590   50   56   76  749 

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
# 1165    1165    1167     904 

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

CROSS_CASES.df <- CROSS_CASES.df[c("THYROIDcancer", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn",
                                   "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, THYROIDcancer = "THYROIDcancer", smoker_never_yn = "smoker_never_yn", 
                               smoker_former_or_never_yn = "smoker_former_or_never_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               NOT_RiskyHeavyDrink_yn = "NOT_RiskyHeavyDrink_yn", HEALTHY_Diet_yn = "HEALTHY_Diet_yn", Not_obese_yn = "Not_obese_yn")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(THYROIDcancer, list(smoker_never_yn , smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Not_obese_yn))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(THYROIDcancer, list(smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Not_obese_yn))))

rownames(cc) <- NULL 
View(cc)



#################
## MODEL TESTS ##
#################

###########################
## 1. Qin baseline model ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(THYROIDcancer ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxneckrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(THYROIDcancer ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxneckrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)

##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = THYROIDcancer ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
                Thyroid_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxneckrtdose.category + epitxn_dose_5.category + 
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + HEALTHY_Diet_yn + Not_obese_yn +
                EAS + AMR + SAS + AFR,
              family = binomial,
              data = dat_all)



summary(fit_all)

# ## With HEI score
# dat_all = PHENO.ANY_SN
# fit_all = glm(formula = THYROIDcancer ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
#                 Thyroid_PRS.tertile.category +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 AGE_AT_DIAGNOSIS + gender + maxneckrtdose.category + epitxn_dose_5.category + 
#                 smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + HEI2015_TOTAL_SCORE.tertile.category + Not_obese_yn +
#                 EAS + AMR + SAS + AFR,
#               family = binomial,
#               data = dat_all)
# 
# 
# 
# summary(fit_all)

##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

###############
## Treatment ##
###############

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all


dat_tx$maxneckrtdose.category [!grepl("Unknown", dat_tx$maxneckrtdose.category)] =
dat_tx$epitxn_dose_5.category [!grepl("Unknown", dat_tx$epitxn_dose_5.category)] = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.549
# 0.55 (Without diet)
# 0.555 [With HEI 2015 (instead of Diet) and other 4 lifestyle factors]
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"
dat_plp.prs$Thyroid_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# 0.416
# 0.417 (Without diet)
# 0.416 [With HEI 2015 (instead of Diet) and other 4 lifestyle factors]
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

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# 0.094
# -0.031 (Without diet)
# -0.138 [With HEI 2015 (instead of Diet) and other 4 lifestyle factors]

#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxneckrtdose.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$maxneckrtdose.category)] =
  dat_tx.plp.prs.lifestyle$epitxn_dose_5.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$epitxn_dose_5.category)] = "None"


## Nullify Genetics
dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Thyroid_PRS.tertile.category = "1st"


## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn)] =
  dat_tx.plp.prs.lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$PhysicalActivity_yn)] =
  dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn)] =
  dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn)] =
  dat_tx.plp.prs.lifestyle$Not_obese_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Not_obese_yn)] = "1"


# ## HEI
# dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn)] =
# dat_tx.plp.prs.lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$PhysicalActivity_yn)] =
# dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn)] =
# dat_tx.plp.prs.lifestyle$Not_obese_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Not_obese_yn)] = "1"
# dat_tx.plp.prs.lifestyle$HEI2015_TOTAL_SCORE.tertile.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$HEI2015_TOTAL_SCORE.tertile.category)] = "3rd"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
# 0.777 With diet
# 0.764 Without diet

