setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/")
# CCSS.org.ANY_SN <- PHENO.ANY_SN # CCSS org
# CCSS.exp.ANY_SN <- PHENO.ANY_SN ## CCSS exp
# sum(colnames(CCSS.org.ANY_SN) == colnames(CCSS.exp.ANY_SN))
# # 54
# # Since columns are same, we can simply rbind the dataframes
# PHENO.ANY_SN <- rbind.data.frame(CCSS.org.ANY_SN, CCSS.exp.ANY_SN)
# rm(list=setdiff(ls(), c("PHENO.ANY_SN")))
# save.image("00.PHENO.ANY_SARCOMA_CCSS_combined_v11.Rdata")

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.PHENO.ANY_SARCOMA_CCSS_combined_v11.Rdata")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

table(PHENO.ANY_SN$CACO)
# 0    1 
# 7839  104

table(SARCOMA$ANY_SN_TYPE) # groupdx3


library(haven)
library(sas7bdat)
library(benchmarkme)
library(dplyr)
library(plyr)
library(data.table)
library (birk)
library(gtools)
library(stringr)
# library(tidyverse)
library(lubridate)
# install.packages("pdftools")
library(pdftools)

# RE: Kiri's email on 02-09-2023
# KIRI.ccss <- read_sas("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS 20200205/combined/others/combinedsn_final.sas7bdat")
# KIRI.ccss.format <- read_sas("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS 20200205/formats.sas7bcat")
# KIRI.ccss <- read_sas("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS 20200205/combined/others/combinedsn_final.sas7bdat", "Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS 20200205/formats.sas7bcat")
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
# 30  7913 



SARCOMA <- PHENO.ANY_SN[PHENO.ANY_SN$SARCOMA == 1,]
KIRI.ccss <- KIRI.ccss[KIRI.ccss$ccssid %in% SARCOMA$ccssid,]

SARCOMA$ANY_SN_TYPE_NEW <- KIRI.ccss$candxo3[match(SARCOMA$KEY, KIRI.ccss$KEY)]
# Not identified: 5146972 and 12127547

write.table(SARCOMA[c("ccssid", "SN_diagnosis_date", "ANY_SN_TYPE", "ANY_SN_TYPE_NEW")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/CCSS_combined_Sarcoma_cases_with_Unspecified_diagnosis.txt", col.names = T, sep = "\t", quote = F, row.names = F)



# SARCOMA$candxo3 <- KIRI.ccss$candxo3[match(SARCOMA$ccssid, KIRI.ccss$ccssid)]
# unique(SARCOMA$candxo3)
# # REF: https://seer.cancer.gov/icd-o-3/sitetype.icdo3.20220429.pdf
# 
# sarcoma.icdo3 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/V11_with_overlapping_PRS_variants_updated_lifestyle_from_Kiri-v2/CCSS_combined/site_types.txt", sep = "\t")
# sarcoma.icdo3$V1 <- gsub("\\/", ".", sarcoma.icdo3$V1)
# 
# SARCOMA$new_SN_diagnosis_labels <- sarcoma.icdo3$V2[match(SARCOMA$candxo3, sarcoma.icdo3$V1)]
# table(SARCOMA$ANY_SN_TYPE)
# 
# SN.question <- SARCOMA[grepl("Other soft tissue sarcomas|Unspecified", SARCOMA$ANY_SN_TYPE),]
# table(is.na(SN.question$new_SN_diagnosis_labels))
# # FALSE  TRUE 
# # 53    38 
# table(is.na(SARCOMA$candxo3))
# as.data.frame(table(SN.question$new_SN_diagnosis_labels))

###########################################
## Check data in each category/cross tab ##
###########################################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- CROSS_CASES.df[c("SARCOMA", "Current_smoker_yn", "PhysicalActivity_yn",
                                   "RiskyHeavyDrink_yn", "Obese_yn",  "aa_class_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, SARCOMA = "SARCOMA",  
                               Current_smoker_yn = "Current_smoker_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               RiskyHeavyDrink_yn = "RiskyHeavyDrink_yn", Obese_yn = "Obese_yn", aa_class_dose_5.category= "aa_class_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(SARCOMA, list( Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, Obese_yn, aa_class_dose_5.category))))

cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(SARCOMA, list( Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, Obese_yn, aa_class_dose_5.category))))

rownames(cc) <- NULL 
# View(cc)


##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SARCOMA ~ Sarcoma_Machiela_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + aa_class_dose_5.category +
                Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
                EAS + AFR,
              family = binomial,
              data = dat_all)


summary(fit_all)

(output <- summary(fit_all)$coefficients)
as.data.frame(apply(output, 2, formatC, format="f", digits=4))
# options(scipen=999)
estimate <- format(round(output[,1],3), nsmall = 3)
std.error <- format(round(output[,2],3), nsmall = 3)
# P.val <- formatC(output[,4], format="G", digits=3)
P.val <- output[,4]
P.val[P.val < 0.001] <- "<0.001"
P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
sarcoma.model <- (setNames(cbind.data.frame(estimate, std.error, P.val
), c("Estimate", "Std.error", "P")))

##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

###############
## Treatment ##
###############

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all

dat_tx$aa_class_dose_5.category = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.026
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_plp.prs$Sarcoma_Machiela_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
# -0.045
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
# 0.494
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$aa_class_dose_5.category = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Sarcoma_Machiela_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$Current_smoker_yn = "No"
dat_tx.plp.prs.lifestyle$PhysicalActivity_yn = "Yes"
dat_tx.plp.prs.lifestyle$RiskyHeavyDrink_yn = "No"
# dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn = "Yes"
dat_tx.plp.prs.lifestyle$Obese_yn = "No"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
# 0.491
SARCOMA.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
SARCOMA.res

all.res <- cbind.data.frame(SN=SN.res, SMN=SMN.res, NMSC=NMSC.res, BREAST=BREAST.res, THYROID=THYROID.res, MENINGIOMA=MENINGIOMA.res, SARCOMA=SARCOMA.res)
# View(all.res)
View(t(all.res))


ChangeNames <- function(x) {
  row.names(x)[grepl("PRS.tertile.category2nd", row.names(x))] <- "PRS.tertile.category2nd"
  row.names(x)[grepl("PRS.tertile.category3rd", row.names(x))] <- "PRS.tertile.category3rd"
  return(x)
}

df_list <- list(sn.model=sn.model, smn.model=smn.model, nmsc.model=nmsc.model, breast.model=breast.model
                ,thyroid.model=thyroid.model, meningioma.model = meningioma.model, sarcoma.model = sarcoma.model)
df_list <- lapply(df_list, ChangeNames)
list2env(df_list ,.GlobalEnv)
merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by, all = TRUE)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}

df <- merge.all(sn.model, smn.model,nmsc.model, breast.model, thyroid.model, meningioma.model, sarcoma.model)
df <- df[!grepl("AGE_AT_LAST_CONTACT", row.names(df)),]
dim(df)

target <- c("(Intercept)", "genderFemale", "AGE_AT_DIAGNOSIS5-9", "AGE_AT_DIAGNOSIS10-14", "AGE_AT_DIAGNOSIS>=15", "AFR", "EAS",
            "PRS.tertile.category2nd","PRS.tertile.category3rd",
            "epitxn_dose_5.category1st","epitxn_dose_5.category2nd","epitxn_dose_5.category3rd","epitxn_dose_5.categoryUnknown",
            "aa_class_dose_5.category1st","aa_class_dose_5.category2nd","aa_class_dose_5.category3rd","aa_class_dose_5.categoryUnknown",
            "anthra_jco_dose_5.category1st","anthra_jco_dose_5.category2nd","anthra_jco_dose_5.category3rd","anthra_jco_dose_5.categoryUnknown",
            "maxabdrtdose.category0-30", "maxabdrtdose.category>=30","maxabdrtdose.categoryUnknown",
            "maxchestrtdose.category0-20","maxchestrtdose.category>=20", "maxchestrtdose.categoryUnknown",
            "maxneckrtdose.category0-11", "maxneckrtdose.category11-20","maxneckrtdose.category20-30","maxneckrtdose.category>=30", "maxneckrtdose.categoryUnknown",
            "maxpelvisrtdose.category0-20", "maxpelvisrtdose.category>=20","maxpelvisrtdose.categoryUnknown",
            "maxsegrtdose.category0-18", "maxsegrtdose.category18-30", "maxsegrtdose.category>=30","maxsegrtdose.categoryUnknown",
            "Current_smoker_ynYes", "Current_smoker_ynUnknown", 
            "RiskyHeavyDrink_ynYes", "RiskyHeavyDrink_ynUnknown",
            "PhysicalActivity_ynNo", "PhysicalActivity_ynUnknown",
            "Obese_ynYes", "Obese_ynUnknown")

df <- df[match(target, rownames(df)),]
