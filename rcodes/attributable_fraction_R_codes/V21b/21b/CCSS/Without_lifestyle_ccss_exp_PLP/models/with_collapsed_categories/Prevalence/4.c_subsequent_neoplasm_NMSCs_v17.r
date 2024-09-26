# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
rm(list=ls())
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

## Add PLP and subset ccss exp
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss_exp_P_LP_zhaoming_qin_v17.Rdata")
PHENO.ANY_SN <- PHENO.ANY_SN[PHENO.ANY_SN$ccssid %in% ccss_exp.samples$V2,]
PHENO.ANY_SN$Zhaoming_carriers <- ccss_exp.samples$Zhaoming_carriers[match(PHENO.ANY_SN$ccssid, ccss_exp.samples$V2)]
PHENO.ANY_SN$Qin_without_Zhaoming_vars_carriers <- ccss_exp.samples$Qin_without_Zhaoming_vars_carriers[match(PHENO.ANY_SN$ccssid, ccss_exp.samples$V2)]


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
## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)


## Read NMSC data from Qi
data1 = read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/sns2022.sas7bdat")
data1=as.data.frame(data1)
# data1$ccssid <- paste0(data1$ccssid, "_", data1$ccssid)
data1$KEY <- paste0(data1$ccssid,":",data1$d_candx)


#########################
## Subsequent neoplasm ##
#########################
## ADD NMSC from Qi
subneo$d_candx <- as.Date(subneo$d_candx, format = "%d%b%Y")
subneo$KEY <- paste0(subneo$ccssid,":",subneo$d_candx)
table(subneo$KEY %in% data1$KEY)
# FALSE  TRUE 
# 6307  3440 
table(data1$KEY %in% subneo$KEY)
# FALSE  TRUE 
# 4629  4434 
subneo$nmsc <- data1$nmsc[match(subneo$KEY, data1$KEY)]

# cc <- cbind.data.frame(subneo$KEY, subneo$nmsc, subneo$AGE.ANY_SN, subneo$groupdx3)

# Now get age of SN after first cancer
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

#########################
## Keep malignant only ##
#########################
subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)

# This will include basal cell, squamous cell and melanoma
NMSCs <- subneo[which(subneo$nmsc ==1),]
# cc <- cbind.data.frame(NMSC$KEY, NMSC$nmsc, NMSC$AGE.ANY_SN, NMSC$groupdx3)

subneo.within5 <- NMSCs[NMSCs$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

# NMSCs <- subneo[grepl("skin", subneo$groupdx3, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(NMSCs)
# 775
## Remove SNs if younger than 18 **
dim(PHENO.ANY_SN)
# 7943   50

PHENO.ANY_SN$AGE.ANY_SN <- NMSCs$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)]
# if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
#   PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
# }


# "a_dx"  : Primary cancer diagnosis age
# "a_end" : age at last contact
# "d_candx" : Date when second cancer was diagnosed                                    
# "a_candx": Age at second cancer diagnosis


# dat[,c("ccssid","strokedt","event","dob","agelstcontact","agedx")]
NMSCs$gradeage <- NMSCs$gradedt
NMSCs$gradedt <- as.Date(NMSCs$d_candx, format = "%d%b%Y")
## Calculate DOB
NMSCs$dob <- NMSCs$gradedt - as.numeric(NMSCs$gradeage) * 365.2422
PHENO.ANY_SN$dob <- NMSCs$dob[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)] ## 2009-02-12
PHENO.ANY_SN$gradedt <- NMSCs$gradedt[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)] ## 2009-02-12



dim(PHENO.ANY_SN)
## 7943 ** END

# Removing samples with SN within the 5 years of childhood cancer **
sum(PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid)
# 1
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid,]
dim(PHENO.ANY_SN)
# 7942   51 ** END

## CA CO status
PHENO.ANY_SN$NMSCs <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% NMSCs$ccssid, 0, 1))
table(PHENO.ANY_SN$NMSCs)
# 0    1 
# 7168  774 


#################
## Missingness ##
#################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxpelvisrtdose.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_tx_missing)
# No  Yes 
# 7421  521 
PHENO.ANY_SN$any_rt_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxpelvisrtdose.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_rt_missing  <- factor(ifelse(PHENO.ANY_SN$any_rt_missing == FALSE, "No", "Yes"))

PHENO.ANY_SN$any_tx_missing <- factor(PHENO.ANY_SN$any_tx_missing, levels = c("No", "Yes"))
PHENO.ANY_SN$any_rt_missing <- factor(PHENO.ANY_SN$any_rt_missing, levels = c("No", "Yes"))
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


################
## Prevalence ##
################

## Qin without Zhaoming

# Output of the table function
table_data <- table(PHENO.ANY_SN$Qin_without_Zhaoming_vars_carriers, PHENO.ANY_SN$NMSCs)

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
prevalence.counts <- sum(PHENO.ANY_SN$Qin_without_Zhaoming_vars_carriers == "Y", na.rm = T)
round((prop.test(prevalence.counts, nrow(PHENO.ANY_SN))$estimate*100), 2)

## Zhaoming ##

# Output of the table function
table_data <- table(PHENO.ANY_SN$Zhaoming_carriers, PHENO.ANY_SN$NMSCs)

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
prevalence.counts <- sum(PHENO.ANY_SN$Zhaoming_carriers  == "Y", na.rm = T)
round((prop.test(prevalence.counts, nrow(PHENO.ANY_SN))$estimate*100), 2)
