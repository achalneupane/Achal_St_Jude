#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_categories_v2.RDATA")
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
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


# Removing samples with with SN within the 5 years of childhood cancer
ANY_SNs <- ANY_SNs[!ANY_SNs$sjlid %in% subneo.within5$sjlid,]
dim(ANY_SNs)
# 605

PHENO.ANY_SN$ANY_SN <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))

#############################
## Add Lifestyle variables ##
#############################
## For each samples, get habits immediately after 18 years of age in agesurvey

# adlthabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.sas7bdat")
adlthabits <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.txt", header = T, sep ="\t")
head(adlthabits)
# remove duplicated rows
adlthabits <- distinct(adlthabits)
## Get DOB
adlthabits$DOB <- PHENO.ANY_SN$dob [match(adlthabits$SJLIFEID, PHENO.ANY_SN$sjlid)]
adlthabits <- adlthabits[!is.na(adlthabits$DOB),]
# change the format of dates YYYY-MM-DD
adlthabits$datecomp <- paste(sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 3), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 1), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 2), sep = "-")
adlthabits$agesurvey <- time_length(interval(as.Date(adlthabits$DOB), as.Date(adlthabits$datecomp)), "years")
# adlthabits$agesurvey <- floor(adlthabits$agesurvey_diff_of_dob_datecomp)

samples.sjlife <- unique(adlthabits$SJLIFEID)
length(samples.sjlife)
# 3571


lifestyle <- {}
for (i in 1:length(samples.sjlife)){
  print(paste0("Doing ", i))
  dat <- adlthabits[adlthabits$SJLIFEID == samples.sjlife[i],]
  if (max(dat$agesurvey) >= 18){
    print("YES")
    dat2 <- dat[dat$agesurvey >= 18,]
    lifestyle.tmp <- dat2[which(dat2$agesurvey == min(dat2$agesurvey)),] # Keep the earliest age after 18 years
  }
  lifestyle <- rbind.data.frame(lifestyle, lifestyle.tmp)
}





# lifestyle$STATUS <- ifelse(lifestyle$SJLIFEID %in% ANY_SNs$sjlid, "case", "control")
# table(lifestyle$STATUS)
# # case control 
# # 592    2978 
# 
# ## Lifestyle.control (Keep the latest variables)
# lifestyle.control <- lifestyle[lifestyle$STATUS == "control",]
# 
# 
# ## Lifestyle.case (non-missing variables before SN diag date)
# lifestyle.case <- lifestyle[lifestyle$STATUS == "case",]


sum(duplicated(lifestyle$SJLIFEID))
# 2
lifestyle$SJLIFEID[duplicated(lifestyle$SJLIFEID)]
# "SJL1080201" "SJL5359215"
## Remove duplicate row
lifestyle <- lifestyle[!duplicated(lifestyle$SJLIFEID),]

## Add all samples
# lifestyle <- cbind.data.frame(wgspop[,1:2], lifestyle[match(wgspop$MRN, lifestyle$mrn), ])
# lifestyle <- lifestyle[-c(3,4)]
# tt <- lifestyle
## Recode categorical variables
lifestyle$relation[lifestyle$relation == 1] <- "Self"
lifestyle$relation[lifestyle$relation == 2] <- "Parent"
lifestyle$relation[lifestyle$relation == 3] <- "Other"

## Recode smoker
lifestyle$smoker[lifestyle$smoker == 1] <- "Past"
lifestyle$smoker[lifestyle$smoker == 2] <- "Current"
lifestyle$smoker[lifestyle$smoker == 3] <- "Never"
lifestyle$smoker_current_yn <- factor(ifelse(lifestyle$smoker != "Current", 0, 1))
lifestyle$smoker_ever_yn <- factor(ifelse(grepl("Never", lifestyle$smoker), 0, 1))

## Recode 1/2 or 0/1 to 0(N) and 1 (Y)
lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink|ltpaw|wtlt|vpa10|yoga",
                colnames(lifestyle))][lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink|ltpaw|wtlt|vpa10|yoga",
                                                      colnames(lifestyle))] == 1 ] <- 1

lifestyle[grepl("nopa|ltpa|ltpaw|wtlt|vpa10|yoga", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa|ltpaw|wtlt|vpa10|yoga", colnames(lifestyle))] == 2 ] <- 0

lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 0 ] <- 0



#######################
## Adolescent habits ##
#######################

adolhabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adolhabits.sas7bdat")
head(adolhabits)

# ###############
# ## Adult BMI ##
# ###############
# Keep the earliest age after 18 years, same as in lifestyle
# adultbmi <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.sas7bdat")
adultbmi <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.txt", sep = "\t", header = T)
head(adultbmi)
adultbmi$AGE <- as.numeric(adultbmi$AGE) ## This age seems to be wrong; for example, SJL1287901 age

# Remove those that is missing DateVisitStart
sum(adultbmi$DateVisitStart == "")
# 149
adultbmi <- adultbmi[adultbmi$DateVisitStart != "",]

## Keep only those BMI data only for the samples in phenotype 
adultbmi <- adultbmi[adultbmi$sjlid %in% PHENO.ANY_SN$sjlid,]; dim(adultbmi)

## Add DOB
adultbmi$DOB <- PHENO.ANY_SN$dob[match(adultbmi$sjlid, PHENO.ANY_SN$sjlid)]

adultbmi$DateVisitStart_edited <-  paste(sapply(strsplit(adultbmi$DateVisitStart, "\\/"), `[`, 3), sapply(strsplit(adultbmi$DateVisitStart, "\\/"), `[`, 1), sapply(strsplit(adultbmi$DateVisitStart, "\\/"), `[`, 2), sep ="-")

# adultbmi$AGE_at_Visit <- floor(time_length(interval(as.Date(adultbmi$DOB), as.Date(adultbmi$DateVisitStart_edited)), "years"))
# table(adultbmi$AGE_at_Visit == adultbmi$AGE) # 125 $AGE seem to be wrong, so calculating the Age as AGE_at_Visit below
# WRONG.AGE <- adultbmi[which(adultbmi$AGE_at_Visit != adultbmi$AGE),]

adultbmi$AGE_at_Visit <- time_length(interval(as.Date(adultbmi$DOB), as.Date(adultbmi$DateVisitStart_edited)), "years")

## Keep the earliest age after 18
samples.sjlife <- unique(adultbmi$sjlid)
length(samples.sjlife)
# 3640


BMI <- {}
for (i in 1:length(samples.sjlife)){
  print(paste0("Doing ", i))
  dat <- adultbmi[adultbmi$sjlid == samples.sjlife[i],]
  if (max(dat$AGE_at_Visit, na.rm = T) >= 18){
    print("YES")
    dat2 <- dat[dat$AGE_at_Visit >= 18,]
    BMI.tmp <- dat2[which(dat2$AGE_at_Visit == min(dat2$AGE_at_Visit, na.rm = T)),] # Keep the earliest age after 18 years
  }
  BMI <- rbind.data.frame(BMI, BMI.tmp)
}

save.BIM <- BMI

BMI <- distinct(BMI)
sum(duplicated(BMI$sjlid))
# 0
BMI$sjlid[duplicated(BMI$sjlid)]


#######################
## Physical Activity ##
#######################

# 1. Create Physical activity based variables every week
# Criteria: 
#   - Physical activity at least 20 mins, three times a week (lifestyle$pa20 > 3). All seem to be doing at least once.
#   - Light physical activity every week (ltpaw == "Y") 
# as.data.frame(ftable(ltpa=lifestyle$ltpa, wtlt=lifestyle$wtlt, vpa10=lifestyle$vpa10, pa20=lifestyle$pa20, yoga=lifestyle$yoga))  
physical.activity <- as.data.frame(ftable(ltpa=lifestyle$ltpa, wtlt=lifestyle$wtlt, vpa10=lifestyle$vpa10, pa20=lifestyle$pa20, yoga=lifestyle$yoga))

# If ltpa, wtlt, vpa10 and yoga, ny three or more have 1, then Yes for PhysicalActivity_YN
lifestyle$PhysicalActivity_YN <- factor(ifelse (rowSums(lifestyle[c("ltpa", "wtlt", "vpa10", "yoga")])>=3, "Y", "N"))

#############
## Obesity ##
#############

BMI$Obesity_YN <- factor(ifelse(BMI$BMI < 30, "N", "Y"))

#############
## Alcohol ##
#############
# If binge, heavy or risky drink any is 1, then Yes for RiskyHeavyDrink_YN
lifestyle$RiskyHeavyDrink_YN <- factor(ifelse (rowSums(lifestyle[c("bingedrink", "heavydrink", "riskydrink")])>=1, "Y", "N"))


##########
## Diet ##
##########
colnames(adultbmi)
# VEGSRV"               
# [19] "FRUITSRV"              "GRAINSRV"              "MEATSRV"               "WGRAINS"               "DAIRYSRV"              "FATSRV"               
# [25] "DT_SODI"               "DT_TFAT"               "DT_CARB"               "DT_PROT"               "M_EGG"                 "AV_TOT_S"             
# [31] "AF_TOT_S"              "R_MEAT_S"              "A_NUT_S"               "A_BEAN_S"
adultbmi$VEGSRV_YN <- factor(ifelse(adultbmi$VEGSRV >= 3, 1, 0)) # Veggie 
adultbmi$FRUITSRV_YN <- factor(ifelse(adultbmi$FRUITSRV >= 3, 1, 0)) # fruits
adultbmi$WGRAINS_YN <- factor(ifelse(adultbmi$WGRAINS >= 3, 1, 0)) # whole grains
adultbmi$DAIRYSRV_YN <- factor(ifelse(adultbmi$DAIRYSRV >= 2.5, 1, 0)) # Dairy
adultbmi$GRAINSRV_YN <- factor(ifelse(adultbmi$GRAINSRV <= 1.5, 1, 0)) # Dairy 


########################################
## Merge BMI, Lifestyle and Phenotype ##
########################################


## Add BMI and Nutrition to Lifestyle. Here, I am extracting the the BMI and Nutrition for the samples
lifestyle$BMI_KEY <- paste(lifestyle$SJLIFEID, lifestyle$agesurvey_floor, sep = ":")

length(unique(adultbmi$sjlid))
# 3640
adultbmi$BMI_KEY <- paste(adultbmi$sjlid, adultbmi$AGE, sep = ":")
sum(lifestyle$BMI_KEY %in% adultbmi$BMI_KEY)
# 2964
## samples that did not match by corresponding age
cc <- lifestyle[!lifestyle$BMI_KEY %in% adultbmi$BMI_KEY,]
lifestyle <- cbind.data.frame(lifestyle, adultbmi[match(lifestyle$BMI_KEY, adultbmi$BMI_KEY), c("BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")])


#########################
## Extract Ethnicities ##
#########################
PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]
PHENO.ANY_SN.AFR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'AFR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]


#################
## MODEL TESTS ##
#################

###########################
## 1. Qin baseline model ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(ANY_SN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

# ## SJLIFE (EUR)
# mod1.EUR <- glm(ANY_SN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
# summary(mod1.EUR)


######################################
## Attributable fraction of Any SNs ##
######################################

dat_all = PHENO.ANY_SN
fit_all = glm(formula = ANY_SN ~ Zhaoming_carriers + Qin_without_Zhaoming_vars_carriers + 
                Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category +
                maxchestrtdose.category + epitxn_dose_5.category, family = binomial,
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

dat_tx$maxsegrtdose.category = dat_tx$maxabdrtdose.category =  dat_tx$maxchestrtdose.category = dat_tx$epitxn_dose_5.category = "None"
dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.455

##########
## P/LP ##
##########
## P/LP Zhaoming and Qin without Zhaoming
dat_plp = dat_all
dat_plp$Zhaoming_carriers = dat_plp$Qin_without_Zhaoming_vars_carriers = "N"

dat_all$pred_no_plp = predict(fit_all, newdata = dat_plp, type = "response")
N_no_plp = sum(dat_all$pred_no_plp, na.rm = TRUE)
af_by_plp_Zhaoming = (N_all - N_no_plp) / N_all
round(af_by_plp_Zhaoming,3)
# 0.022

#########
## PRS ##
#########
dat_prs = dat_all
dat_prs$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

dat_all$pred_no_Pleiotropy_PRSWEB_PRS.tertile.category = predict(fit_all, newdata = dat_prs, type = "response")
N_no_pred_no_Pleiotropy_PRSWEB_PRS.tertile.category = sum(dat_all$pred_no_Pleiotropy_PRSWEB_PRS.tertile.category, na.rm = TRUE)
af_by_N_no_pred_no_Pleiotropy_PRSWEB_PRS.tertile.category = (N_all - N_no_pred_no_Pleiotropy_PRSWEB_PRS.tertile.category) / N_all
round(af_by_N_no_pred_no_Pleiotropy_PRSWEB_PRS.tertile.category, 3)
# 0.085
