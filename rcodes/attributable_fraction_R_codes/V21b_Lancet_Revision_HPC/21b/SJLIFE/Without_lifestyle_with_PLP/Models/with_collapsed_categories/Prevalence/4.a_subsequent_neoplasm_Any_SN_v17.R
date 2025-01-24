## **
# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
rm(list=ls())
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
library(survival)

# subneo <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/subneo.sas7bdat")
# subneo <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_07_22_2023_subneo.txt", header = T, stringsAsFactors = F)

subneo <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/subneoplasms.sas7bdat')


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

subneo$gradedt <- as.Date(subneo$gradedt, "%m/%d/%Y") ## **
subneo$AGE.ANY_SN <- time_length(interval(as.Date(subneo$DOB), as.Date(subneo$gradedt)), "years")

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
# Get SNs for the first time and Age at First SN. Note: I am only keeping first event SN 5 years post diagnosis
# For this, I will first sort the table by date
library(data.table)

## Remove SNs as cases that are within 5 years of primary diagnosis
subneo <- subneo[!subneo$sjlid %in% subneo.within5$sjlid,]
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


PHENO.ANY_SN$gradedt <- ANY_SNs$gradedt[match(PHENO.ANY_SN$sjlid, ANY_SNs$sjlid)]
PHENO.ANY_SN$AGE.ANY_SN <- ANY_SNs$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, ANY_SNs$sjlid)]

# ## Remove SNs if younger than 18 **
# if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
# PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
# }

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]


PHENO.ANY_SN$ANY_SNs <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))

table(PHENO.ANY_SN$ANY_SNs)
# 0    1 
# 3774  605 


#################
## missingness ##
#################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

PHENO.ANY_SN$any_chemo_missing <- apply(PHENO.ANY_SN[c("epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_chemo_missing  <- factor(ifelse(PHENO.ANY_SN$any_chemo_missing == FALSE, "No", "Yes"))

PHENO.ANY_SN$any_rt_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_rt_missing  <- factor(ifelse(PHENO.ANY_SN$any_rt_missing == FALSE, "No", "Yes"))

PHENO.ANY_SN$any_tx_missing <- factor(PHENO.ANY_SN$any_tx_missing, levels = c("No", "Yes"))
PHENO.ANY_SN$any_chemo_missing <- factor(PHENO.ANY_SN$any_chemo_missing, levels = c("No", "Yes"))
PHENO.ANY_SN$any_rt_missing <- factor(PHENO.ANY_SN$any_rt_missing, levels = c("No", "Yes"))

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


############################################################
## Drop Unknown level from the lifestyle factor variables ##
############################################################
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
## Prevalence ##
################

## Qin without Zhaoming

# Output of the table function
table_data <- table(PHENO.ANY_SN$Qin_without_Zhaoming_vars_carriers, PHENO.ANY_SN$ANY_SNs)

# Calculate the total number of individuals with and without Any_SNs
total_no_any_sns <- sum(table_data[, "0"])
total_with_any_sns <- sum(table_data[, "1"])

# Calculate the number of carriers with and without Any_SNs
carriers_no_any_sns <- table_data["Y", "0"]
carriers_with_any_sns <- table_data["Y", "1"]

# Calculate prevalence of carriers among those with and without Any_SNs
prevalence_carriers_no_any_sns <- (carriers_no_any_sns / total_no_any_sns) * 100
prevalence_carriers_with_any_sns <- (carriers_with_any_sns / total_with_any_sns) * 100

# Print results
round(prevalence_carriers_with_any_sns,2)
round(prevalence_carriers_no_any_sns,2)

# prevalence all
prevalence.counts <- sum(PHENO.ANY_SN$Qin_without_Zhaoming_vars.Non.Ref.Counts > 0)
table(ifelse(PHENO.ANY_SN$Qin_without_Zhaoming_vars.Non.Ref.Counts > 0, "Y", "N"))
round((prop.test(prevalence.counts, nrow(PHENO.ANY_SN))$estimate*100), 2)

## Zhaoming ##

# Output of the table function
table_data <- table(PHENO.ANY_SN$Zhaoming_carriers, PHENO.ANY_SN$ANY_SNs)

# Calculate the total number of individuals with and without Any_SNs
total_no_any_sns <- sum(table_data[, "0"])
total_with_any_sns <- sum(table_data[, "1"])

# Calculate the number of carriers with and without Any_SNs
carriers_no_any_sns <- table_data["Y", "0"]
carriers_with_any_sns <- table_data["Y", "1"]

# Calculate prevalence of carriers among those with and without Any_SNs
prevalence_carriers_no_any_sns <- (carriers_no_any_sns / total_no_any_sns) * 100
prevalence_carriers_with_any_sns <- (carriers_with_any_sns / total_with_any_sns) * 100

# Print results
round(prevalence_carriers_with_any_sns,2)
round(prevalence_carriers_no_any_sns,2)

# prevalence all
prevalence.counts <- sum(PHENO.ANY_SN$Zhaoming_Non.Ref.Counts > 0)
table(ifelse(PHENO.ANY_SN$Zhaoming_Non.Ref.Counts > 0, "Y", "N"))
round((prop.test(prevalence.counts, nrow(PHENO.ANY_SN))$estimate*100), 2)

