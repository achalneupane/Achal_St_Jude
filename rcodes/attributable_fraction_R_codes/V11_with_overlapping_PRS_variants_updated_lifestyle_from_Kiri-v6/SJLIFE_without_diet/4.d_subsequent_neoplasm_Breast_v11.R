#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11-6.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)
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
## GET SN 18 or older
subneo <- subneo[subneo$AGE.ANY_SN >= 18,]

ALL.LIFESTYLE[c("PhysicalActivity_yn_agesurvey", "smoker_former_or_never_yn_agesurvey", "NOT_RiskyHeavyDrink_yn_agesurvey",
                "Not_obese_yn_agesurvey", "HEALTHY_Diet_yn_agesurvey")] <- sapply(ALL.LIFESTYLE[c("PhysicalActivity_yn_agesurvey", "smoker_former_or_never_yn_agesurvey", "NOT_RiskyHeavyDrink_yn_agesurvey",
                                                                                                  "Not_obese_yn_agesurvey", "HEALTHY_Diet_yn_agesurvey")], floor)

ALL.LIFESTYLE$SURVEY_MIN <- apply(ALL.LIFESTYLE[c("PhysicalActivity_yn_agesurvey", "smoker_former_or_never_yn_agesurvey", "NOT_RiskyHeavyDrink_yn_agesurvey",
                                                  "Not_obese_yn_agesurvey")], 1, min)

# Anyone before the first survey, remove them
subneo$survey_First <- ALL.LIFESTYLE$SURVEY_MIN[match(subneo$sjlid,ALL.LIFESTYLE$SJLIFEID)]
table(subneo$survey_First > subneo$AGE.ANY_SN)
# FALSE  TRUE 
# 981   613 

ALL.LIFESTYLE$PhysicalActivity_yn[which(ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey != ALL.LIFESTYLE$SURVEY_MIN)] <- NA
ALL.LIFESTYLE$Current_smoker_yn[which(ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey != ALL.LIFESTYLE$SURVEY_MIN)] <- NA
ALL.LIFESTYLE$RiskyHeavyDrink_yn[which(ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey != ALL.LIFESTYLE$SURVEY_MIN)] <- NA
ALL.LIFESTYLE$Obese_yn[which(ALL.LIFESTYLE$Not_obese_yn_agesurvey != ALL.LIFESTYLE$SURVEY_MIN)] <- NA
# ALL.LIFESTYLE$HEALTHY_Diet_yn[which(ALL.LIFESTYLE$HEALTHY_Diet_yn_agesurvey != ALL.LIFESTYLE$SURVEY_MIN)] <- NA

subneo <- subneo[!subneo$survey_First > subneo$AGE.ANY_SN,]
##############

BREASTcancer <- subneo[grepl("breast", subneo$diag, ignore.case = T),]
BREASTcancer <- setDT(BREASTcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(BREASTcancer)
# 58
table(BREASTcancer$diaggrp)
# Removing samples with SNs within 5 years of childhood cancer
BREASTcancer <- BREASTcancer[!BREASTcancer$sjlid %in% subneo.within5$sjlid,]
nrow(BREASTcancer)
# 57
PHENO.ANY_SN$BREASTcancer <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% BREASTcancer$sjlid, 0, 1))

#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% BREASTcancer$sjlid, 0, 1))

## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- BREASTcancer$gradedt[match(ALL.LIFESTYLE$SJLIFEID, BREASTcancer$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- BREASTcancer$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, BREASTcancer$sjlid)]

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("Current_smoker_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("PhysicalActivity_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("RiskyHeavyDrink_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEALTHY_Diet_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEALTHY_Diet_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$Not_obese_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("Obese_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEI2015_TOTAL_SCORE_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEI2015_TOTAL_SCORE")] <- NA

############################
## Add lifestyle to Pheno ##
############################
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2015_TOTAL_SCORE", "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Obese_yn")])
# PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("smoker_never_yn", "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Obese_yn")])

# Count missing
PHENO.ANY_SN$missing.lifestyles <- rowSums(is.na(PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Obese_yn")]))
table(PHENO.ANY_SN$missing.lifestyles)

## Relevel 6 lifestyle variables
PHENO.ANY_SN$Current_smoker_yn[is.na(PHENO.ANY_SN$Current_smoker_yn)] <- "Unknown"
PHENO.ANY_SN$Current_smoker_yn <- factor(PHENO.ANY_SN$Current_smoker_yn, level = c("No", "Yes", "Unknown")) 

PHENO.ANY_SN$PhysicalActivity_yn[is.na(PHENO.ANY_SN$PhysicalActivity_yn)] <- "Unknown"
PHENO.ANY_SN$PhysicalActivity_yn <- factor(PHENO.ANY_SN$PhysicalActivity_yn, level = c("Yes", "No", "Unknown")) 

PHENO.ANY_SN$RiskyHeavyDrink_yn[is.na(PHENO.ANY_SN$RiskyHeavyDrink_yn)] <- "Unknown"
PHENO.ANY_SN$RiskyHeavyDrink_yn <- factor(PHENO.ANY_SN$RiskyHeavyDrink_yn, level = c("No", "Yes", "Unknown")) 

PHENO.ANY_SN$HEALTHY_Diet_yn[is.na(PHENO.ANY_SN$HEALTHY_Diet_yn)] <- "Unknown"
PHENO.ANY_SN$HEALTHY_Diet_yn <- factor(PHENO.ANY_SN$HEALTHY_Diet_yn, level = c("Yes", "No", "Unknown")) 

PHENO.ANY_SN$Obese_yn[is.na(PHENO.ANY_SN$Obese_yn)] <- "Unknown";
PHENO.ANY_SN$Obese_yn <- factor(PHENO.ANY_SN$Obese_yn, level = c("No", "Yes", "Unknown")) 

#########################
## Create HEI tertiles ##
#########################
HEI.to.categorize <- c("HEI2015_TOTAL_SCORE")

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

CROSS_CASES.df <- CROSS_CASES.df[c("BREASTcancer", "Current_smoker_yn", "PhysicalActivity_yn",
                                   "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", Obese_yn = "Obese_yn")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, BREASTcancer = "BREASTcancer", 
                               Current_smoker_yn = "Current_smoker_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               RiskyHeavyDrink_yn = "RiskyHeavyDrink_yn", HEALTHY_Diet_yn = "HEALTHY_Diet_yn", Obese_yn = "Obese_yn")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(BREASTcancer, list(Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Obese_yn))))

cc <- (as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(BREASTcancer, list(Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Obese_yn)))))

rownames(cc) <- NULL 
# View(cc)

############################
## Attributable Fractions ##
############################

## -------------------------------------- PRS 2019
dat_all = PHENO.ANY_SN
fit_all = glm(formula = BREASTcancer ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
                Mavaddat_2019_ER_POS_Breast_PRS.tertile.category +
                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                Mavaddat_2019_ER_NEG_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2+ 
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category +
                Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
                EAS + AFR,
              family = binomial,
              data = dat_all)

summary(fit_all)

# ## WIth Diet
# fit_all = glm(formula = BREASTcancer ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
#                 Mavaddat_2019_ER_POS_Breast_PRS.tertile.category +
#                 Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
#                 Mavaddat_2019_ER_NEG_Breast_PRS.tertile.category +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2+ 
#                 AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
#                 maxchestrtdose.category + anthra_jco_dose_5.category +
#                 Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + HEALTHY_Diet_yn + Obese_yn +
#                 EAS + AFR,
#               family = binomial,
#               data = dat_all)
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

dat_tx$maxchestrtdose.category =
dat_tx$anthra_jco_dose_5.category = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx, 3)
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"
dat_plp.prs$Mavaddat_2019_ER_POS_Breast_PRS.tertile.category = dat_plp.prs$Mavaddat_2019_ER_NEG_Breast_PRS.tertile.category = dat_plp.prs$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$Current_smoker_yn = "No"
dat_lifestyle$PhysicalActivity_yn = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn = "No"
# dat_lifestyle$HEALTHY_Diet_yn = "Yes"
dat_lifestyle$Obese_yn = "No"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)

#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxchestrtdose.category =
  dat_tx.plp.prs.lifestyle$anthra_jco_dose_5.category = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Mavaddat_2019_ER_POS_Breast_PRS.tertile.category = dat_tx.plp.prs.lifestyle$Mavaddat_2019_ER_NEG_Breast_PRS.tertile.category = dat_tx.plp.prs.lifestyle$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_lifestyle$Current_smoker_yn = "No"
dat_lifestyle$PhysicalActivity_yn = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn = "No"
# dat_lifestyle$HEALTHY_Diet_yn = "Yes"
dat_lifestyle$Obese_yn = "No"



dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)

BREAST.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
BREAST.res
# 0.487  0.318 -0.304  0.560
