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
MENINGIOMA <- subneo[grepl("meningioma", subneo$groupdx3, ignore.case = T),]
MENINGIOMA <- setDT(MENINGIOMA)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(MENINGIOMA)
# 256

## Remove SNs if younger than 18 **
dim(PHENO.ANY_SN)
# 7943   50

PHENO.ANY_SN$AGE.ANY_SN <- MENINGIOMA$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, MENINGIOMA$ccssid)]
if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
  PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
}

dim(PHENO.ANY_SN)
## 7932 51 ** END

# Removing samples with SN within the 5 years of childhood cancer **
sum(PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid)
# 25
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid,]
dim(PHENO.ANY_SN)
# 7907 ** END

## CA CO status
PHENO.ANY_SN$MENINGIOMA <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% MENINGIOMA$ccssid, 0, 1))
table(PHENO.ANY_SN$MENINGIOMA)
# 0    1 
# 7663  244


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
# [1] 7823   52

sum((PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey >= 18|
       PHENO.ANY_SN$PhysicalActivity_yn_agesurvey >= 18|
       PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey >= 18|
       PHENO.ANY_SN$Not_obese_yn_agesurvey >= 18), na.rm = T)
# 7823


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
# 7796  53
######################### ** END


## Add any missing to each lifestyle variable
# PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn")]
PHENO.ANY_SN$any_lifestyle_missing <- apply(PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_lifestyle_missing  <- factor(ifelse(PHENO.ANY_SN$any_lifestyle_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_lifestyle_missing)
# No  Yes 
# 69 7727
########################################
## Do the same for missing treatments ##
########################################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_tx_missing)
# No  Yes 
# 7155  641
#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
ethnicity.admixture$INDIVIDUAL <- sapply(strsplit(ethnicity.admixture$INDIVIDUAL,"_"), `[`, 1)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


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
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.MENINGIOMA.V14.Rdata")

