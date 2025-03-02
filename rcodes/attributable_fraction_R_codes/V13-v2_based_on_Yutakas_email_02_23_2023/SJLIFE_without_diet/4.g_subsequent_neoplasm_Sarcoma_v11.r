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
## Sarcoma ##
#############
SARCOMA <- subneo[grepl("sarcoma", subneo$diag, ignore.case = T),]
SARCOMA <- setDT(SARCOMA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(SARCOMA)
# 35
# Removing samples with SNs within 5 years of childhood cancer
SARCOMA <- SARCOMA[!SARCOMA$sjlid %in% subneo.within5$sjlid,]
nrow(SARCOMA)
# 32

## Samples with either retinoblastoma or germ cell tumor
# retino.germ <- c("SJL5002218", "SJL5004118", "SJL5013318", "SJL5013918", "SJL5015118", "SJL5016518", "SJL5025005", "SJL5029618", "SJL5029918", "SJL5033218")
# SARCOMA <- SARCOMA[!SARCOMA$sjlid %in% retino.germ,]

PHENO.ANY_SN$SARCOMA <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% SARCOMA$sjlid, 0, 1))

#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% SARCOMA$sjlid, 0, 1))

## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- SARCOMA$gradedt[match(ALL.LIFESTYLE$SJLIFEID, SARCOMA$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- SARCOMA$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, SARCOMA$sjlid)]

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("Current_smoker_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("PhysicalActivity_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("RiskyHeavyDrink_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEALTHY_Diet_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEALTHY_Diet_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$Not_obese_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("Obese_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEI2015_TOTAL_SCORE_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEI2015_TOTAL_SCORE")] <- NA

#############################
## Addd lifestyle to Pheno ##
#############################
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
# 3rd     2nd     1st Unknown 
# 1175    1175    1176     875 

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

CROSS_CASES.df <- CROSS_CASES.df[c("SARCOMA", "Current_smoker_yn", "PhysicalActivity_yn",
                                   "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", Obese_yn = "Obese_yn", aa_class_dose_5.category= "aa_class_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, SARCOMA = "SARCOMA",  
                               Current_smoker_yn = "Current_smoker_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
                               RiskyHeavyDrink_yn = "RiskyHeavyDrink_yn", HEALTHY_Diet_yn = "HEALTHY_Diet_yn", Obese_yn = "Obese_yn", aa_class_dose_5.category= "aa_class_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(SARCOMA, list( Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Obese_yn, aa_class_dose_5.category))))

cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(SARCOMA, list( Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Obese_yn, aa_class_dose_5.category))))

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
# 0.369
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
# 0.067
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
# 0.774
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
# 0.877
SARCOMA.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
SARCOMA.res
# 0.373 0.023 0.541 0.744

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
            "maxabdrtdose.category>0-<30", "maxabdrtdose.category>=30","maxabdrtdose.categoryUnknown",
            "maxchestrtdose.category>0-<20","maxchestrtdose.category>=20", "maxchestrtdose.categoryUnknown",
            "maxneckrtdose.category>0-<11", "maxneckrtdose.category>=11-<20","maxneckrtdose.category>=20-<30","maxneckrtdose.category>=30", "maxneckrtdose.categoryUnknown",
            "maxpelvisrtdose.category>0-<20", "maxpelvisrtdose.category>=20","maxpelvisrtdose.categoryUnknown",
            "maxsegrtdose.category>0-<18", "maxsegrtdose.category>=18-<30", "maxsegrtdose.category>=30","maxsegrtdose.categoryUnknown",
            "Current_smoker_ynYes", "Current_smoker_ynUnknown", 
            "RiskyHeavyDrink_ynYes", "RiskyHeavyDrink_ynUnknown",
            "PhysicalActivity_ynNo", "PhysicalActivity_ynUnknown",
            "Obese_ynYes", "Obese_ynUnknown")

df <- df[match(target, rownames(df)),]
