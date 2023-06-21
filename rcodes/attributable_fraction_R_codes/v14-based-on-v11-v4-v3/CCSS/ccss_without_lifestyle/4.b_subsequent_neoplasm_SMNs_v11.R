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
SMNs <- subneo[!grepl("skin", subneo$groupdx3, ignore.case = T),]
SMNs <- setDT(SMNs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]

## Remove SNs if younger than 18 **
dim(PHENO.ANY_SN)
# 7943   50

PHENO.ANY_SN$AGE.ANY_SN <- SMNs$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, SMNs$ccssid)]
if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
  PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
}

dim(PHENO.ANY_SN)
## 7870 51 ** END

# Removing samples with SN within the 5 years of childhood cancer **
sum(PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid)
# 7
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid,]
dim(PHENO.ANY_SN)
# 7863 ** END

## CA CO status
PHENO.ANY_SN$SMNs <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% SMNs$ccssid, 0, 1))
table(PHENO.ANY_SN$SMNs)
# 0    1 
# 6307 1556 



######################### **

########################################
## Do the same for missing treatments ##
########################################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_tx_missing)


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


##########################################
## Find out benign SMNs and remove them ##
##########################################
# # based on Yadav's email on 03/09/2023, I am removing all benign diagnoses
## This file is from Kyla
KIRI.ccss <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Kyla/combinedsn_final_02_17_2023.csv", header = T, sep = ",", stringsAsFactors = F)
dim(KIRI.ccss)
## Keep non-missing candxo3
KIRI.ccss <- KIRI.ccss[!is.na(KIRI.ccss$candxo3),]
# KIRI.ccss <- KIRI.ccss[KIRI.ccss$candxo3 !="",]
KIRI.ccss <- KIRI.ccss[KIRI.ccss$d_candx !="",]
dim(KIRI.ccss)
KIRI.ccss$SN_diagnosis_date <- as.Date(KIRI.ccss$d_candx, format = "%d%b%Y")
KIRI.ccss$SN_diagnosis_date <- format(KIRI.ccss$SN_diagnosis_date, "%m-%d-%Y") # 06-30-2008
KIRI.ccss$KEY <- paste(KIRI.ccss$ccssid, KIRI.ccss$SN_diagnosis_date, sep = ":")

PHENO.ANY_SN$SN_diagnosis_date <- as.Date(PHENO.ANY_SN$d_candx, format = "%d%b%Y")
PHENO.ANY_SN$SN_diagnosis_date  <- format(PHENO.ANY_SN$SN_diagnosis_date, "%m-%d-%Y") # "08-23-2016"
PHENO.ANY_SN$ccssid <- gsub("_.*","",PHENO.ANY_SN$ccssid)
PHENO.ANY_SN$KEY <- paste(PHENO.ANY_SN$ccssid, PHENO.ANY_SN$SN_diagnosis_date, sep = ":")


table(PHENO.ANY_SN$KEY %in% KIRI.ccss$KEY)
# FALSE  TRUE 
# 6326   1537



SMNs <- PHENO.ANY_SN[PHENO.ANY_SN$SMNs == 1,]
KIRI.ccss <- KIRI.ccss[KIRI.ccss$KEY %in% SMNs$KEY,]

KIRI.ccss <- KIRI.ccss[grepl("\\/0|\\/1", KIRI.ccss$candxo3),]

PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% KIRI.ccss$ccssid,]
table(PHENO.ANY_SN$SMNs)
# 6307 1276
##########################################

################
## Cross tabs ##
################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- PHENO.ANY_SN

CROSS_CASES.df <- CROSS_CASES.df[,c("SMNs", "maxsegrtdose.category", "maxchestrtdose.category", "maxabdrtdose.category", "epitxn_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, SMNs = "SMNs", 
                               maxsegrtdose.category = "maxsegrtdose.category", maxchestrtdose.category = "maxchestrtdose.category", 
                               maxabdrtdose.category = "maxabdrtdose.category", epitxn_dose_5.category = "epitxn_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(SMNs, list(maxsegrtdose.category, maxchestrtdose.category, maxabdrtdose.category, epitxn_dose_5.category))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(SMNs, list(maxsegrtdose.category, maxchestrtdose.category, maxabdrtdose.category, epitxn_dose_5.category))))

rownames(cc) <- NULL 
cc


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.SMNs_without_lifestyle.V14.Rdata")

