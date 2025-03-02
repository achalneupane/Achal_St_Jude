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
################
## Meningioma ##
################
# Get MENINGIOMA for the first time and Age at First MENINGIOMA.
# For this, I will first sort the table by date
library(data.table)
MENINGIOMA <- subneo[grepl("meningioma", subneo$diag, ignore.case = T),]
MENINGIOMA <- setDT(MENINGIOMA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
MENINGIOMA <- MENINGIOMA[!MENINGIOMA$sjlid %in% subneo.within5$sjlid,]
nrow(MENINGIOMA)
nrow(MENINGIOMA)
# 149

## Remove MENINGIOMA if younger than 18
PHENO.ANY_SN$AGE.ANY_SN <- MENINGIOMA$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, MENINGIOMA$sjlid)]
if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
}

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]

PHENO.ANY_SN$MENINGIOMA <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% MENINGIOMA$sjlid, 0, 1))

table(PHENO.ANY_SN$MENINGIOMA)
# 0    1 
# 4230 139

#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% MENINGIOMA$sjlid, 0, 1))



## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- MENINGIOMA$gradedt[match(ALL.LIFESTYLE$SJLIFEID, MENINGIOMA$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- MENINGIOMA$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, MENINGIOMA$sjlid)]


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
# 3627

## Remove any samples that do not have lifestyle
PHENO.ANY_SN  <- PHENO.ANY_SN[PHENO.ANY_SN$sjlid %in% ALL.LIFESTYLE$SJLIFEID,]

#############################
## Addd lifestyle to Pheno ##
#############################
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2015_TOTAL_SCORE", "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Obese_yn")])


## Add any missing to each lifestyle variable
# PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn")]
PHENO.ANY_SN$any_lifestyle_missing <- apply(PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_lifestyle_missing  <- factor(ifelse(PHENO.ANY_SN$any_lifestyle_missing == FALSE, "No", "Yes"))

########################################
## Do the same for missing treatments ##
########################################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


############################################################
## Drop Unknown level from the lifestyle factor variables ##
############################################################
## Recode lifestyle variables to fit the model for missingness
## Missing lifestyle
PHENO.ANY_SN$Current_smoker_yn[PHENO.ANY_SN$Current_smoker_yn == "Unknown"] <- "No"
PHENO.ANY_SN$Current_smoker_yn <- droplevels(PHENO.ANY_SN$Current_smoker_yn)

PHENO.ANY_SN$PhysicalActivity_yn[PHENO.ANY_SN$PhysicalActivity_yn == "Unknown"] <- "Yes"
PHENO.ANY_SN$PhysicalActivity_yn <- droplevels(PHENO.ANY_SN$PhysicalActivity_yn)

PHENO.ANY_SN$RiskyHeavyDrink_yn[PHENO.ANY_SN$RiskyHeavyDrink_yn == "Unknown"] <- "No"
PHENO.ANY_SN$RiskyHeavyDrink_yn <- droplevels(PHENO.ANY_SN$RiskyHeavyDrink_yn)

PHENO.ANY_SN$Obese_yn[PHENO.ANY_SN$Obese_yn == "Unknown"] <- "No"
PHENO.ANY_SN$Obese_yn <- droplevels(PHENO.ANY_SN$Obese_yn)


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

CROSS_CASES.df <- CROSS_CASES.df[,c("MENINGIOMA", "maxsegrtdose.category", "epitxn_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, MENINGIOMA = "MENINGIOMA", 
                               maxsegrtdose.category = "maxsegrtdose.category", epitxn_dose_5.category = "epitxn_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(MENINGIOMA, list(maxsegrtdose.category, epitxn_dose_5.category))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(MENINGIOMA, list(maxsegrtdose.category, epitxn_dose_5.category))))

rownames(cc) <- NULL 
cc


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_diet.MENINGIOMA.V14-4-3.Rdata")
