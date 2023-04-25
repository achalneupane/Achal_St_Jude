# Yutakas email on 03/08/2023: I think this result looks good to me.  The only
# thing is that large (positive or negative) coefficients are clustered in
# "Unknown" treatment categories.  This is presumably due to the same reason as
# the lifestyle unknowns where "unknowns" of different items are overlapping.
# In the case of treatment, if the survivor's medical record abstraction is
# missing, all treatment variables will be unknown.  Thus, I suggest you do the
# same for the treatment "Unknowns" as you did lifestyle "Unknowns", i.e.,
# create an indicator variable "Any tx Unknown" and put everyone with "Unknown"
# at any of the treatment variables there.  Then, you will remove "Unknown" from
# each of the treatment variables and we will not have the large coefficients (I
# think).
#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
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
#############
## Any SNs ##
#############
# Get THYROIDcancer for the first time and Age at First THYROIDcancer.
# For this, I will first sort the table by date
library(data.table)
THYROIDcancer <- subneo[grepl("thyroid", subneo$diag, ignore.case = T),]
THYROIDcancer <- setDT(THYROIDcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
THYROIDcancer <- THYROIDcancer[!THYROIDcancer$sjlid %in% subneo.within5$sjlid,]
nrow(THYROIDcancer)
# 86

## Remove THYROIDcancer if younger than 18
PHENO.ANY_SN$AGE.ANY_SN <- THYROIDcancer$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, THYROIDcancer$sjlid)]
if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
}

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]

PHENO.ANY_SN$THYROIDcancer <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% THYROIDcancer$sjlid, 0, 1))

table(PHENO.ANY_SN$THYROIDcancer)
# 0    1 
# 4293 81

#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% THYROIDcancer$sjlid, 0, 1))



## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- THYROIDcancer$gradedt[match(ALL.LIFESTYLE$SJLIFEID, THYROIDcancer$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- THYROIDcancer$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, THYROIDcancer$sjlid)]


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

ALL.LIFESTYLE <- ALL.LIFESTYLE[!(ALL.LIFESTYLE$Current_smoker_yn == "Unknown" &
                                   ALL.LIFESTYLE$PhysicalActivity_yn == "Unknown" &
                                   ALL.LIFESTYLE$RiskyHeavyDrink_yn == "Unknown" &
                                   ALL.LIFESTYLE$Obese_yn == "Unknown" ),]

dim(ALL.LIFESTYLE)
# [1] 3692   25

sum((ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey >= 18|
       ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey >= 18|
       ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey >= 18|
       ALL.LIFESTYLE$Not_obese_yn_agesurvey >= 18), na.rm = T)
# 3692


# ALL.LIFESTYLE.test <- ALL.LIFESTYLE


ALL.LIFESTYLE <- ALL.LIFESTYLE[which(ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey >= 18 |
                                       ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey >= 18 |
                                       ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey >= 18 |
                                       ALL.LIFESTYLE$Not_obese_yn_agesurvey >= 18),]

cols <- c(
  "smoker_former_or_never_yn_agesurvey",
  "PhysicalActivity_yn_agesurvey",
  "NOT_RiskyHeavyDrink_yn_agesurvey",
  "Not_obese_yn_agesurvey"
)

## round to nearest integer
# saved.cc <- ALL.LIFESTYLE[, c("SJLIFEID", cols)]
ALL.LIFESTYLE[, cols] <- apply(ALL.LIFESTYLE[, cols], 2, round)
library(matrixStats)
ALL.LIFESTYLE$survey_min <- rowMins(as.matrix(ALL.LIFESTYLE[, cols]), na.rm = TRUE)
# cc <- ALL.LIFESTYLE[, c("SJLIFEID", cols, "survey_min")]
ALL.LIFESTYLE$Current_smoker_yn [which(ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey != ALL.LIFESTYLE$survey_min)] <- "Unknown"
ALL.LIFESTYLE$PhysicalActivity_yn [which(ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey != ALL.LIFESTYLE$survey_min)] <- "Unknown"
ALL.LIFESTYLE$RiskyHeavyDrink_yn [which(ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey != ALL.LIFESTYLE$survey_min)] <- "Unknown"
ALL.LIFESTYLE$Obese_yn [which(ALL.LIFESTYLE$Not_obese_yn_agesurvey != ALL.LIFESTYLE$survey_min)] <- "Unknown"


to.remove <- ALL.LIFESTYLE$SJLIFEID[which(ALL.LIFESTYLE$survey_min > ALL.LIFESTYLE$AGE.ANY_SN)]
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% to.remove,]

sum(PHENO.ANY_SN$sjlid %in% ALL.LIFESTYLE$SJLIFEID)
# 3649

## Remove any samples that do not have lifestyle
PHENO.ANY_SN  <- PHENO.ANY_SN[PHENO.ANY_SN$sjlid %in% ALL.LIFESTYLE$SJLIFEID,]

#############################
## Addd lifestyle to Pheno ##
#############################
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2015_TOTAL_SCORE", "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Obese_yn")])


##################################
## Imputation of missing values ##
##################################
# wanted items
wanted.cols <- c("Current_smoker_yn", "PhysicalActivity_yn", "Obese_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")

ALL.LIFESTYLE.IMPUTE <- ALL.LIFESTYLE[c("SJLIFEID", wanted.cols)]

ALL.LIFESTYLE.IMPUTE[c("Current_smoker_yn", "PhysicalActivity_yn", "Obese_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")] <- lapply(ALL.LIFESTYLE.IMPUTE[c("Current_smoker_yn", "PhysicalActivity_yn", "Obese_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")], function(x) ifelse(x == "Unknown", NA, as.character(x)))

## Making healthy outcomes as 1 else 0
ALL.LIFESTYLE.IMPUTE$Current_smoker_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$Current_smoker_yn == "No", 1, 0)
ALL.LIFESTYLE.IMPUTE$RiskyHeavyDrink_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$RiskyHeavyDrink_yn == "No", 1, 0)
ALL.LIFESTYLE.IMPUTE$Obese_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$Obese_yn == "No", 1, 0)
ALL.LIFESTYLE.IMPUTE$PhysicalActivity_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$PhysicalActivity_yn == "Yes", 1, 0)
ALL.LIFESTYLE.IMPUTE$HEALTHY_Diet_yn <- ifelse(ALL.LIFESTYLE.IMPUTE$HEALTHY_Diet_yn == "Yes", 1, 0)

## Number of items available to get scores
ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score <- rowSums(!is.na(ALL.LIFESTYLE.IMPUTE[wanted.cols]))

# If all 5 items are available; they were simply summed up to get the score 
ALL.LIFESTYLE.IMPUTE$SCORE [ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 5] <-  rowSums(ALL.LIFESTYLE.IMPUTE[ALL.LIFESTYLE.IMPUTE$number_of_items_to_get_score == 5, wanted.cols])

## For missing items, per Qi: "For people with missing items, mean of non-missing items were taken and then multiplied by the # of items to get the score"
ALL.LIFESTYLE.IMPUTE$mean_of_non_missing_items <- rowMeans(ALL.LIFESTYLE.IMPUTE[wanted.cols], na.rm = T)
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
wanted.cols <- wanted.cols <- c("Current_smoker_yn", "PhysicalActivity_yn", "Obese_yn", "RiskyHeavyDrink_yn")
ALL.LIFESTYLE.IMPUTE <- ALL.LIFESTYLE[c("SJLIFEID", wanted.cols)]

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
## Add missing treatments ##
############################
## Add any missing to each tx variable
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxneckrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_tx_missing)

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


############################################################
## Drop Unknown level from the lifestyle factor variables ##
############################################################
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


# Create a cross-tabulation table between maxsegrtdose and maxchedtrtdose for cases
cases_table <- table(Max_SegmentedRT_Dose = PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$THYROIDcancer == 1],
                     Max_ChestRT_Dose = PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$THYROIDcancer == 1])

# Create a cross-tabulation table between maxsegrtdose and maxchedtrtdose for controls
control_table <- table(Max_SegmentedRT_Dose = PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$THYROIDcancer == 0],
                       Max_ChestRT_Dose = PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$THYROIDcancer == 0])


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_diet.THYROIDcancer.V15_with_composite_lifestyle.Rdata")
