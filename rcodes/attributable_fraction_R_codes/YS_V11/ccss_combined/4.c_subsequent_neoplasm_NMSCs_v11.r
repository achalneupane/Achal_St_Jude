## Load CCSS original
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_org_Genetic_data_P_LP_V11.Rdata")

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

## Read NMSC data from Qi
data1 = read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/sns2022.sas7bdat")
data1=as.data.frame(data1)
data1$ccssid <- paste0(data1$ccssid, "_", data1$ccssid)
data1$KEY <- paste0(data1$ccssid,":",data1$d_candx)


#########################
## Subsequent neoplasm ##
#########################
## ADD NMSC from Qi
subneo$d_candx <- as.Date(subneo$d_candx, format = "%d%b%Y")
subneo$KEY <- paste0(subneo$ccssid,":",subneo$d_candx)
table(subneo$KEY %in% data1$KEY)
# FALSE  TRUE 
# 3646  3023
table(data1$KEY %in% subneo$KEY)
# FALSE  TRUE 
# 5163  3900 
subneo$nmsc <- data1$nmsc[match(subneo$KEY, data1$KEY)]

# cc <- cbind.data.frame(subneo$KEY, subneo$nmsc, subneo$AGE.ANY_SN, subneo$groupdx3)

# Now get age of SN after first cancer
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1351

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$ccssid))
# 19

###########
## NMSCs ##
###########
# This will include basal cell, squamous cell and melanoma

NMSCs <- subneo[which((subneo$nmsc ==1| (subneo$nmsc == 2 & subneo$groupdx3 == "Skin"))),]
# cc <- cbind.data.frame(NMSC$KEY, NMSC$nmsc, NMSC$AGE.ANY_SN, NMSC$groupdx3)

# NMSCs <- subneo[grepl("skin", subneo$groupdx3, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(NMSCs)
# 678

# Removing samples with SNs within 5 years of childhood cancer
NMSCs <- NMSCs[!NMSCs$ccssid %in% subneo.within5$ccssid,]
nrow(NMSCs)
# 672
PHENO.ANY_SN$NMSC <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% NMSCs$ccssid, 0, 1))
table(PHENO.ANY_SN$NMSC)
# 0    1 
# 4335  672 

PHENO.ANY_SN$AGE.ANY_SN <- NMSCs$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)]
#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
PHENO.ANY_SN$CACO <- PHENO.ANY_SN$NMSC

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "smoker_former_or_never_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "PhysicalActivity_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "NOT_RiskyHeavyDrink_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "Not_obese_yn"] <- "Unknown"

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])

CCSS.org.ANY_SN <- PHENO.ANY_SN # CCSS org

# Load CCSS expansion
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_exp_Genetic_data_P_LP_V11.Rdata")

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
# 2661   417
table(data1$KEY %in% subneo$KEY)
# FALSE  TRUE 
# 8529   534
subneo$nmsc <- data1$nmsc[match(subneo$KEY, data1$KEY)]

# cc <- cbind.data.frame(subneo$KEY, subneo$nmsc, subneo$AGE.ANY_SN, subneo$groupdx3)

# Now get age of SN after first cancer
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 269

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$ccssid))
# 8


###########
## NMSCs ##
###########
# This will include basal cell, squamous cell and melanoma

NMSCs <- subneo[which((subneo$nmsc ==1| (subneo$nmsc == 2 & subneo$groupdx3 == "Skin"))),]
# cc <- cbind.data.frame(NMSC$KEY, NMSC$nmsc, NMSC$AGE.ANY_SN, NMSC$groupdx3)

# NMSCs <- subneo[grepl("skin", subneo$groupdx3, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(NMSCs)
# 97

# Removing samples with SNs within 5 years of childhood cancer
NMSCs <- NMSCs[!NMSCs$ccssid %in% subneo.within5$ccssid,]
nrow(NMSCs)
# 97
PHENO.ANY_SN$NMSC <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% NMSCs$ccssid, 0, 1))
table(PHENO.ANY_SN$NMSC)
# 0    1 
# 2839   97 

PHENO.ANY_SN$AGE.ANY_SN <- NMSCs$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)]
#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
PHENO.ANY_SN$CACO <- PHENO.ANY_SN$NMSC


## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "smoker_former_or_never_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "PhysicalActivity_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "NOT_RiskyHeavyDrink_yn"] <- "Unknown"
PHENO.ANY_SN[which(PHENO.ANY_SN$CACO == 1 & (PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), "Not_obese_yn"] <- "Unknown"


#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])

CCSS.exp.ANY_SN <- PHENO.ANY_SN ## CCSS exp
sum(colnames(CCSS.org.ANY_SN) == colnames(CCSS.exp.ANY_SN))
# 53
# Since columns are same, we can simply rbind the dataframes
PHENO.ANY_SN <- rbind.data.frame(CCSS.org.ANY_SN, CCSS.exp.ANY_SN)
rm(list=setdiff(ls(), c("PHENO.ANY_SN")))

table(PHENO.ANY_SN$CACO)
# 0    1 
# 7174  769 

## Distribution of overlapping "Unknown" lifestyle factors
source("https://raw.githubusercontent.com/achalneupane/Achal_St_Jude/main/rcodes/attributable_fraction_R_codes/YS_V11_without_diet/get_missing_combination.R")
get_missing_combinations(PHENO.ANY_SN[c("smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "Not_obese_yn")])

###########################################
## Check data in each category/cross tab ##
###########################################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- CROSS_CASES.df[c("NMSC", "smoker_former_or_never_yn", "PhysicalActivity_yn",
                                   "NOT_RiskyHeavyDrink_yn", Not_obese_yn = "Not_obese_yn")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, NMSC = "NMSC", 
                               smoker_former_or_never_yn = "smoker_former_or_never_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               NOT_RiskyHeavyDrink_yn = "NOT_RiskyHeavyDrink_yn", Not_obese_yn = "Not_obese_yn")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(NMSC, list(smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, Not_obese_yn))))

cc.NMSC <- as.data.frame(t(CROSS_CASES.df %>%
                             cross_cases(NMSC, list(smoker_former_or_never_yn, PhysicalActivity_yn, NOT_RiskyHeavyDrink_yn, Not_obese_yn))))


rownames(cc.NMSC) <- NULL 
##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = NMSC ~ BASALcell_PRS.tertile.category + SQUAMOUScell_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + Not_obese_yn +
                EAS + AFR,
              family = binomial,
              data = dat_all)

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

dat_tx$maxsegrtdose.category =
dat_tx$maxabdrtdose.category =
dat_tx$maxpelvisrtdose.category = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.297
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"
dat_plp.prs$SQUAMOUScell_PRS.tertile.category = dat_plp.prs$BASALcell_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# 0.189
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$smoker_former_or_never_yn =
dat_lifestyle$PhysicalActivity_yn =
dat_lifestyle$NOT_RiskyHeavyDrink_yn =
dat_lifestyle$Not_obese_yn = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# 0.311
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxsegrtdose.category =
  dat_tx.plp.prs.lifestyle$maxabdrtdose.category =
  dat_tx.plp.prs.lifestyle$maxpelvisrtdose.category = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$SQUAMOUScell_PRS.tertile.category = dat_tx.plp.prs.lifestyle$BASALcell_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn =
dat_tx.plp.prs.lifestyle$PhysicalActivity_yn =
dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn =
dat_tx.plp.prs.lifestyle$Not_obese_yn = "1"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
# 0.702

NMSC.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
NMSC.res
