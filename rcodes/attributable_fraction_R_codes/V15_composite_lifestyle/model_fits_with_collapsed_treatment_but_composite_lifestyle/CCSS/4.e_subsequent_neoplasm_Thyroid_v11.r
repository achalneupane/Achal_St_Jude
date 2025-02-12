rm()
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")

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

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$ccssid))
# 26


#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
THYROIDcancer <- subneo[grepl("thyroid", subneo$groupdx3, ignore.case = T),]
THYROIDcancer <- setDT(THYROIDcancer)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(THYROIDcancer)
# 167

## Remove SNs if younger than 18 **
dim(PHENO.ANY_SN)
# 7943   50

PHENO.ANY_SN$AGE.ANY_SN <- THYROIDcancer$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, THYROIDcancer$ccssid)]
if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
  PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
}

dim(PHENO.ANY_SN)
## 7934 51 ** END

# Removing samples with SN within the 5 years of childhood cancer **
sum(PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid,]
dim(PHENO.ANY_SN)
# 7912 ** END

## CA CO status
PHENO.ANY_SN$THYROIDcancer <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% THYROIDcancer$ccssid, 0, 1))
table(PHENO.ANY_SN$THYROIDcancer)
# 0    1 
# 7755  157


######################### **
############################################################################################
############################################################################################
#################################### Work for V11-4-v2 #####################################
############################################################################################
############################################################################################
# Yutaka's email 02/23/2023: The number of SN should decrease
# in the new analysis because those who developed SN before the first adult
# survey should be out of the analysis.  Also, we should use the "any missing in
# the lifestyle variables" rather than using the individual missing (with
# missing combined to the reference in each lifestyle variable).

## remove all lifestyle missing
PHENO.ANY_SN <- PHENO.ANY_SN[!(PHENO.ANY_SN$Current_smoker_yn == "Unknown" &
                                 PHENO.ANY_SN$PhysicalActivity_yn == "Unknown" &
                                 PHENO.ANY_SN$RiskyHeavyDrink_yn == "Unknown" &
                                 PHENO.ANY_SN$Obese_yn == "Unknown" ),]

dim(PHENO.ANY_SN)
# [1] 7828   52

sum((PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey >= 18|
       PHENO.ANY_SN$PhysicalActivity_yn_agesurvey >= 18|
       PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey >= 18|
       PHENO.ANY_SN$Not_obese_yn_agesurvey >= 18), na.rm = T)
# 7828


PHENO.ANY_SN <- PHENO.ANY_SN[which(PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$PhysicalActivity_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$Not_obese_yn_agesurvey >= 18),]

cols <- c(
  "smoker_former_or_never_yn_agesurvey",
  "PhysicalActivity_yn_agesurvey",
  "NOT_RiskyHeavyDrink_yn_agesurvey",
  "Not_obese_yn_agesurvey"
)



## round to nearest integer
# saved.cc <- ALL.LIFESTYLE[, c("SJLIFEID", cols)]
# ALL.LIFESTYLE[, cols] <- apply(PHENO.ANY_SN[, cols], 2, round)
library(matrixStats)
PHENO.ANY_SN$survey_min <- rowMins(as.matrix(PHENO.ANY_SN[, cols]), na.rm = TRUE)
# cc <- PHENO.ANY_SN[, c("SJLIFEID", cols, "survey_min")]
PHENO.ANY_SN$Current_smoker_yn [which(PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"
PHENO.ANY_SN$PhysicalActivity_yn [which(PHENO.ANY_SN$PhysicalActivity_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"
PHENO.ANY_SN$RiskyHeavyDrink_yn [which(PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"
PHENO.ANY_SN$Obese_yn [which(PHENO.ANY_SN$Not_obese_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"

## Remove SN cases if the diagnosis date is prior to the youngest adult survey date
PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$survey_min > PHENO.ANY_SN$AGE.ANY_SN),]
dim(PHENO.ANY_SN)
# 7806  53
######################### ** END


##################################
## Imputation of missing values ##
##################################
# wanted items
wanted.cols <- wanted.cols <- c("Current_smoker_yn", "PhysicalActivity_yn", "Obese_yn", "RiskyHeavyDrink_yn")

ALL.LIFESTYLE.IMPUTE <- PHENO.ANY_SN[c("ccssid", wanted.cols)]

## create composite without diet

ALL.LIFESTYLE.IMPUTE[c("Current_smoker_yn", "PhysicalActivity_yn", "Obese_yn", "RiskyHeavyDrink_yn")] <- lapply(ALL.LIFESTYLE.IMPUTE[c("Current_smoker_yn", "PhysicalActivity_yn", "Obese_yn", "RiskyHeavyDrink_yn")], function(x) ifelse(x == "Unknown", NA, as.character(x)))

## Making healthy outcomes as 1 else 0
ALL.LIFESTYLE.IMPUTE$Current_smoker_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$Current_smoker_yn == "No", 1, 0)
ALL.LIFESTYLE.IMPUTE$RiskyHeavyDrink_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$RiskyHeavyDrink_yn == "No", 1, 0)
ALL.LIFESTYLE.IMPUTE$Obese_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$Obese_yn == "No", 1, 0)
ALL.LIFESTYLE.IMPUTE$PhysicalActivity_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$PhysicalActivity_yn == "Yes", 1, 0)

## Number of items available to get scores
ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score <- rowSums(!is.na(ALL.LIFESTYLE.IMPUTE[wanted.cols]))

# If all 4 (after removing HEALTHY_Diet_yn) items are available; they were simply summed up to get the score 
ALL.LIFESTYLE.IMPUTE$SCORE [ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 4] <-  rowSums(ALL.LIFESTYLE.IMPUTE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 4, wanted.cols])

## For missing items, per Qi: "For people with missing items, mean of non-missing items were taken and then multiplied by the # of items to get the score"
ALL.LIFESTYLE.IMPUTE$mean_of_non_missing_items <- rowMeans(ALL.LIFESTYLE.IMPUTE[wanted.cols], na.rm = T)
ALL.LIFESTYLE.IMPUTE$imputedSCORE <- ALL.LIFESTYLE.IMPUTE$mean_of_non_missing_items * ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score

# replace the score values with imputed score for missing ones; this was applied to rows where items available (non-missing ones) are >= 2 & < 4 
ALL.LIFESTYLE.IMPUTE$SCORE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score >= 2 & ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score < 4] <- ALL.LIFESTYLE.IMPUTE$imputedSCORE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score >= 2 & ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score < 4]

# Define favorable, intermediate and unfavorable----
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[ALL.LIFESTYLE.IMPUTE$SCORE>= 3] <- "favorable" # favorable [at least 3 of the five health lifestyle factors]
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[ALL.LIFESTYLE.IMPUTE$SCORE == 2] <- "intermediate" # intermediate [two healthy lifestyle factors] 
ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[ALL.LIFESTYLE.IMPUTE$SCORE <= 1] <- "unfavorable" # unfavorable [no or only one healthy lifestyle factor] 

## Add to PHENO
PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET <- ALL.LIFESTYLE.IMPUTE$LIFESTYLE_STATUS_WO_DIET[match(PHENO.ANY_SN$ccssid, ALL.LIFESTYLE.IMPUTE$ccssid)]

# relevel factors
PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET[is.na(PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET)] <- "Unknown"
PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET <- factor(PHENO.ANY_SN$LIFESTYLE_STATUS_WO_DIET, level = c("favorable", "intermediate", "unfavorable", "Unknown"))

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxneckrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_tx_missing)

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
ethnicity.admixture$INDIVIDUAL <- sapply(strsplit(ethnicity.admixture$INDIVIDUAL,"_"), `[`, 1)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


#####################################################
## Drop Unknown level from the tx factor variables ##
#####################################################
## Recode tx variables to fit the model for missingness
## Missing tx
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxsegrtdose.category <- droplevels(PHENO.ANY_SN$maxsegrtdose.category)

PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxneckrtdose.category <- droplevels(PHENO.ANY_SN$maxneckrtdose.category)

PHENO.ANY_SN$maxabdrtdose.category[PHENO.ANY_SN$maxabdrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxabdrtdose.category <- droplevels(PHENO.ANY_SN$maxabdrtdose.category)

PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxchestrtdose.category <- droplevels(PHENO.ANY_SN$maxchestrtdose.category)

PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxpelvisrtdose.category <- droplevels(PHENO.ANY_SN$maxpelvisrtdose.category)

PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$epitxn_dose_5.category <- droplevels(PHENO.ANY_SN$epitxn_dose_5.category)

PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$anthra_jco_dose_5.category <- droplevels(PHENO.ANY_SN$anthra_jco_dose_5.category)

PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$aa_class_dose_5.category <- droplevels(PHENO.ANY_SN$aa_class_dose_5.category)


################
## Cross tabs ##
################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- PHENO.ANY_SN

CROSS_CASES.df <- CROSS_CASES.df[,c("THYROIDcancer", "maxneckrtdose.category", "epitxn_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, THYROIDcancer = "THYROIDcancer", 
                               maxneckrtdose.category = "maxneckrtdose.category", epitxn_dose_5.category = "epitxn_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(THYROIDcancer, list(maxneckrtdose.category, epitxn_dose_5.category))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(THYROIDcancer, list(maxneckrtdose.category, epitxn_dose_5.category))))

rownames(cc) <- NULL 
cc


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.THYROIDcancer.V15_with_composite_lifestyle.Rdata")

