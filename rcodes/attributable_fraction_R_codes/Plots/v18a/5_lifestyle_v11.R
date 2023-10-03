
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

#########################
## Load Phenotype data ##
#########################
rm(list=ls())
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_categories_v11.RDATA")

##########################################################################################################
## Recode Physical activity and Smoker based on what Kiri Ness' suggested (Using SAS code leanmass.sas) ##
##########################################################################################################
# get additional adult health habits
adult_habbits <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/adult_healthhabits.sas7bdat')

adult_habbits$sjlid <- adult_habbits$SJLIFEID

# remove duplicated rows
adult_habbits <- distinct(adult_habbits)
## Get DOB
adult_habbits$DOB <- PHENO.ANY_SN$dob [match(adult_habbits$SJLIFEID, PHENO.ANY_SN$sjlid)]
# adult_habbits <- adult_habbits[!is.na(adult_habbits$DOB),]
# change the format of dates YYYY-MM-DD
adult_habbits$agesurvey2 <- time_length(interval(as.Date(adult_habbits$DOB), as.Date(adult_habbits$datecomp)), "years")
adult_habbits$agesurvey2[is.na(adult_habbits$agesurvey2)] <- adult_habbits$agesurvey[is.na(adult_habbits$agesurvey2)]
adult_habbits$agesurvey <- adult_habbits$agesurvey2

lifestyle <- adult_habbits

lifestyle$sjlid <- lifestyle$SJLIFEID 
lifestyle <- lifestyle[lifestyle$sjlid %in% PHENO.ANY_SN$sjlid,]


## 1.-------------- Physical activity

## First, work on vpa10
# if vpa10=1 and vpadays=. then do vpadays=1; end;
lifestyle$vpadays[which(lifestyle$vpa10 == 1 & is.na(lifestyle$vpadays))] <- 1

# if vpa10=1 and vpamin=. then do vpamin=10; end;
lifestyle$vpamin[which(lifestyle$vpa10 == 1 & is.na(lifestyle$vpamin))] <- 10
# cc <- lifestyle[which(lifestyle$vpa10 == 1 & is.na(lifestyle$vpamin)),]
max(lifestyle$vpamin, na.rm = T)
# 720
# if vpamin>360 then do vpamin=360; end; /*Cap six hours per day*/; So, cap this to six hours only
lifestyle$vpamin[which(lifestyle$vpamin > 360)] <- 360

lifestyle$wvpa <- NA
lifestyle$wvpa[which(lifestyle$vpa10 == 1)] <- lifestyle$vpadays[which(lifestyle$vpa10 == 1)] * lifestyle$vpamin[which(lifestyle$vpa10 == 1)]
lifestyle$wvpa[which(lifestyle$vpa10 == 2)] <- 0

# if wvpa=. then do;---
# if vpa10 in (.,2) and nopa=1 and pa20 in (.,0) then do wvpa=0; end;
# if vpa10 in (.,2) and pa20 not in (.,0) then do wvpa=pa20*20; end;
# if vpa10 in (.,2) and nopa=2 and pa20 in (.,0) then do wvpa=0; end;
lifestyle$wvpa[which(is.na(lifestyle$wvpa) & (is.na(lifestyle$vpa10) | lifestyle$vpa10 ==2) & lifestyle$nopa ==1 & (is.na(lifestyle$pa20) | lifestyle$pa20 ==0))] <- 0

index <- which(is.na(lifestyle$wvpa) & (is.na(lifestyle$vpa10) | lifestyle$vpa10 ==2) & (!is.na(lifestyle$pa20) | lifestyle$pa20 !=0))
lifestyle$wvpa[index] <- lifestyle$pa20[index] *20

lifestyle$wvpa[which(is.na(lifestyle$wvpa) & (is.na(lifestyle$vpa10) | lifestyle$vpa10 ==2) & lifestyle$nopa ==2 & (is.na(lifestyle$pa20) | lifestyle$pa20 ==0))] <- 0

## Now work on mpa
# if mpa10=. and (mpadays ne . or mpamin ne .) then do mpa10=1; end;
lifestyle$mpa10[which(is.na(lifestyle$mpa10) & (!is.na(lifestyle$mpadays) | !is.na(lifestyle$mpamin)))] <- 1
# if mpa10=1 and mpadays=. then do mpadays=1; end;
lifestyle$mpadays[which((lifestyle$mpa10==1) & (is.na(lifestyle$mpadays)))] <- 1
# if mpa10=1 and mpamin=. then do mpamin=10; end;
lifestyle$mpamin[which((lifestyle$mpa10==1) & (is.na(lifestyle$mpamin)))] <- 1
# if mpamin>360 then do mpamin=360; end; /*Cap six hours per day*/
lifestyle$mpamin[which(lifestyle$mpamin > 360)] <- 360

# if mpa10=1 then do wmpa=mpadays*mpamin; end;
lifestyle$wmpa <- NA
lifestyle$wmpa[which(lifestyle$mpa10 == 1)] <- lifestyle$mpadays[which(lifestyle$mpa10 == 1)] * lifestyle$mpamin[which(lifestyle$mpa10 == 1)] 
# if mpa10=2 then do wmpa=0; end;
lifestyle$wmpa[which(lifestyle$mpa10 == 2)] <- 0
# if wmpa=. then do;---
# if mpa10 in (.,2) and nopa=1 then do wmpa=0; end;
# if wvpa ne . then do wmpa=0; end;
lifestyle$wmpa[which(is.na(lifestyle$wmpa) & ((is.na(lifestyle$mpa10)|lifestyle$mpa10 == 2) & lifestyle$nopa == 1))] <- 0
lifestyle$wmpa[which(is.na(lifestyle$wmpa) & !is.na(lifestyle$wvpa)) ] <- 0

# mvpawk=(wmpa)+(wvpa*2); if mvpawk>2520 then mvpawk=2520; /*cap six hours per day*/
lifestyle$mvpawk <- lifestyle$wmpa + (lifestyle$wvpa *2)
lifestyle$mvpawk[which(lifestyle$mvpawk > 2520)] <- 2520
lifestyle$CDC_PA <- ifelse(lifestyle$mvpawk < 150, 0,1) # Check this line


## Physical activity
MET_iid_dob_18 = subset(lifestyle, agesurvey >= 18)
MET_iid_dob_18 <- MET_iid_dob_18[!is.na(MET_iid_dob_18$CDC_PA), ]
MET_iid_dob_18_sorted = MET_iid_dob_18[order(MET_iid_dob_18$sjlid, MET_iid_dob_18$agesurvey, decreasing = FALSE),]
MET_iid_dob_18_uniq = MET_iid_dob_18_sorted[!duplicated(MET_iid_dob_18_sorted$sjlid),]

MET_iid_dob_18_uniq$PhysicalActivity_yn <- as.numeric(ifelse(MET_iid_dob_18_uniq$CDC_PA == 1, 1, 0))



## 2.-------------- Smoker status
# Get evsm, smnow, smnvr, cigmo, cigd, smyr
lifestyle$evsm <- adult_habbits$evsm[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$smnow <- adult_habbits$smnow[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$cigmo <- adult_habbits$cigmo[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$cigd <- adult_habbits$cigd[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$smyr <- adult_habbits$smyr[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]



# if evsm=2 then do ; smnow=2 ; cigmo=2 ; cigd=0 ; smyr=0 ; end ;
lifestyle$smnow[which(lifestyle$evsm == 2)] <- 2
lifestyle$cigmo[which(lifestyle$evsm == 2)] <- 2
lifestyle$cigd[which(lifestyle$evsm == 2)] <- 0
lifestyle$smyr[which(lifestyle$evsm == 2)] <- 0

# lifestyle$smnvr <- NA
# # if smnvr=1 then do ; evsm=2 ; smnow=2 ; cigmo=2 ; cigd=0 ; smyr=0 ; end ;
# lifestyle$evsm[which(lifestyle$smnvr == 1)] <- 2
# lifestyle$smnow[which(lifestyle$smnvr == 1)] <- 2
# lifestyle$cigmo[which(lifestyle$smnvr == 1)] <- 2
# lifestyle$cigd[which(lifestyle$smnvr == 1)] <- 0
# lifestyle$smyr[which(lifestyle$smnvr == 1)] <- 0

# if evsm=. and smnow=1 then do ; evsm=1 ; end ;
lifestyle$evsm[which(is.na(lifestyle$evsm) & lifestyle$smnow == 1)] <- 1

# if evsm=1 and smnow=. and cigmo=1 then smnow=1
lifestyle$smnow[which(lifestyle$evsm == 1 & is.na(lifestyle$smnow))] <- 1

# if evsm=1 and smnow=. and cigmo=2 then smnow=2
lifestyle$smnow[which(lifestyle$evsm == 1 & is.na(lifestyle$smnow) & lifestyle$cigmo == 2)] <- 2


lifestyle$smkStat <- NA
lifestyle$smkStat[which(lifestyle$evsm == 1 & lifestyle$smnow == 2)] <- 1 # Former
lifestyle$smkStat[which(lifestyle$evsm == 1 & lifestyle$smnow == 1)] <- 2 # Current
lifestyle$smkStat[which(lifestyle$evsm == 2)] <- 3 # Never

table(lifestyle$smkStat)

## smoker status
smk_iid_dob_18 = subset(lifestyle, agesurvey >= 18)
smk_iid_dob_18 <- smk_iid_dob_18[!is.na(smk_iid_dob_18$smkStat), ]
smk_iid_dob_18_sorted = smk_iid_dob_18[order(smk_iid_dob_18$sjlid, smk_iid_dob_18$agesurvey, decreasing = FALSE),]
smk_iid_dob_18_uniq = smk_iid_dob_18_sorted[!duplicated(smk_iid_dob_18_sorted$sjlid),]

smk_iid_dob_18_uniq$smoker_former_or_never_yn <- as.numeric(ifelse(smk_iid_dob_18_uniq$smkStat != 2, 1, 0))
smk_iid_dob_18_uniq$smoker_never_yn <- as.numeric(ifelse(smk_iid_dob_18_uniq$smkStat == 3, 1, 0))



## 3.-----------Drinking
drk_iid_dob_18 = subset(lifestyle, agesurvey >= 18)
# If heavy or risky drink any is 1, then Yes for RiskyHeavyDrink_yn
drk_iid_dob_18$NOT_RiskyHeavyDrink_yn <- as.numeric(ifelse (rowSums(drk_iid_dob_18[c("heavydrink", "riskydrink")])==0, 1, 0))
drk_iid_dob_18 <- drk_iid_dob_18[!is.na(drk_iid_dob_18$NOT_RiskyHeavyDrink_yn), ]
drk_iid_dob_18_sorted = drk_iid_dob_18[order(drk_iid_dob_18$sjlid, drk_iid_dob_18$agesurvey, decreasing = FALSE),]
drk_iid_dob_18_uniq = drk_iid_dob_18_sorted[!duplicated(drk_iid_dob_18_sorted$sjlid),]


## 4.-----------Adult BMI
# This BMI from Yadav has more samples, so using this for BMI only
library(sas7bdat)
bmi = read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/function_combo_basic.sas7bdat')
iid = read.table('Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/4401_attributable_fraction_ids.txt', header = FALSE)
bmi_iid = subset(bmi, sjlid %in% iid$V1, select = c('sjlid', 'assmntdate', 'BMIadj'))
demo = read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat')
dob = demo[c('sjlid', 'dob')]
bmi_iid_dob = merge(bmi_iid, dob, by='sjlid')
bmi_iid_dob$assmntdate = bmi_iid_dob$assmntdate
bmi_iid_dob$agebmi = as.numeric(difftime(as.Date(bmi_iid_dob$assmntdate), as.Date(bmi_iid_dob$dob))/365.25)
bmi_iid_dob_18 = subset(bmi_iid_dob, agebmi>=18)
bmi_iid_dob_18_sorted = bmi_iid_dob_18[order(bmi_iid_dob_18$sjlid, bmi_iid_dob_18$agebmi, decreasing = FALSE),]
bmi_iid_dob_18_uniq = bmi_iid_dob_18_sorted[!duplicated(bmi_iid_dob_18_sorted$sjlid),]

bmi_iid_dob_18_uniq$Not_obese_yn <- as.numeric(ifelse(bmi_iid_dob_18_uniq$BMIadj < 30, 1, 0))

## Part 2.------------Adult habits from Siddhant
adlthabits <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.txt", header = T, sep ="\t")
head(adlthabits)
adlthabits$sjlid <- adlthabits$SJLIFEID
## Fix data format
adlthabits$datecomp <- paste(sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 3), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 1), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 2), sep = "-")

# remove duplicated rows
adlthabits <- distinct(adlthabits)
## Get DOB
adlthabits$DOB <- PHENO.ANY_SN$dob [match(adlthabits$SJLIFEID, PHENO.ANY_SN$sjlid)]
# adlthabits <- adlthabits[!is.na(adlthabits$DOB),]
# change the format of dates YYYY-MM-DD
adlthabits$agesurvey2 <- time_length(interval(as.Date(adlthabits$DOB), as.Date(adlthabits$datecomp)), "years")
adlthabits$agesurvey2[is.na(adlthabits$agesurvey2)] <- adlthabits$agesurvey[is.na(adlthabits$agesurvey2)]
adlthabits$agesurvey <- adlthabits$agesurvey2

adlthabits <- adlthabits[adlthabits$sjlid %in% PHENO.ANY_SN$sjlid,]


## 2.-------------- Smoker status ## from Siddhant
smk_iid_dob_18.2 = subset(adlthabits, agesurvey >= 18)
smk_iid_dob_18.2 <- smk_iid_dob_18.2[!is.na(smk_iid_dob_18.2$smoker), ]
smk_iid_dob_18_sorted.2 = smk_iid_dob_18.2[order(smk_iid_dob_18.2$sjlid, smk_iid_dob_18.2$agesurvey, decreasing = FALSE),]
smk_iid_dob_18_uniq.2 = smk_iid_dob_18_sorted.2[!duplicated(smk_iid_dob_18_sorted.2$sjlid),]

## Recode smoker
# adlthabits$smoker[adlthabits$smoker == 1] <- "Past"
# adlthabits$smoker[adlthabits$smoker == 2] <- "Current"
# adlthabits$smoker[adlthabits$smoker == 3] <- "Never"
smk_iid_dob_18_uniq.2$smoker_former_or_never_yn <- as.numeric(ifelse(smk_iid_dob_18_uniq.2$smoker != 2, 1, 0))
smk_iid_dob_18_uniq.2$smoker_never_yn <- as.numeric(ifelse(smk_iid_dob_18_uniq.2$smoker == 3, 1, 0))

table(smk_iid_dob_18_uniq.2$smoker)
table(smk_iid_dob_18_uniq$smkStat)

# PHENO.ANY_SN$smk_iid_dob_18_uniq.2 <- smk_iid_dob_18_uniq.2$smoker[match(PHENO.ANY_SN$sjlid, smk_iid_dob_18_uniq.2$sjlid)]
# PHENO.ANY_SN$smk_iid_dob_18_uniq <- smk_iid_dob_18_uniq$smkStat[match(PHENO.ANY_SN$sjlid, smk_iid_dob_18_uniq$sjlid)]

## smk_iid_dob_18_uniq is what I created and smk_iid_dob_18_uniq.2 is received from Siddhant; Use smk_iid_dob_18_uniq.2
# table(PHENO.ANY_SN$smk_iid_dob_18_uniq.2 == PHENO.ANY_SN$smk_iid_dob_18_uniq) 


## 3.-----------Drinking ## from Siddhant
adlthabits[grepl("bingedrink|heavydrink|riskydrink", colnames(adlthabits))][adlthabits[grepl("bingedrink|heavydrink|riskydrink", colnames(adlthabits))] == 0 ] <- 0

drk_iid_dob_18.2 = subset(adlthabits, agesurvey >= 18)
# If heavy or risky drink any is 1, then Yes for RiskyHeavyDrink_yn
drk_iid_dob_18.2$NOT_RiskyHeavyDrink_yn <- as.numeric(ifelse (rowSums(drk_iid_dob_18.2[c("heavydrink", "riskydrink")])==0, 1, 0))
drk_iid_dob_18.2 <- drk_iid_dob_18.2[!is.na(drk_iid_dob_18.2$NOT_RiskyHeavyDrink_yn), ]
drk_iid_dob_18_sorted.2 = drk_iid_dob_18[order(drk_iid_dob_18.2$sjlid, drk_iid_dob_18.2$agesurvey, decreasing = FALSE),]
drk_iid_dob_18_uniq.2 = drk_iid_dob_18_sorted.2[!duplicated(drk_iid_dob_18_sorted.2$sjlid),]

# PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn.2 <- drk_iid_dob_18_uniq.2$NOT_RiskyHeavyDrink_yn[match(PHENO.ANY_SN$sjlid, drk_iid_dob_18_uniq.2$sjlid)]
# PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn <- drk_iid_dob_18_uniq$NOT_RiskyHeavyDrink_yn[match(PHENO.ANY_SN$sjlid, drk_iid_dob_18_uniq$sjlid)]

## drk_iid_dob_18_uniq is what I created and drk_iid_dob_18_uniq.2 is received from Siddhant; Use smk_iid_dob_18_uniq.2
# table(PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn.2 == PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn) 

##############################
## More lifestyle variables ##
##############################

## 5.------------Diet
# Keep the earliest age after 18 years, same as in lifestyle; Extracting diet here
# adultdiet <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.sas7bdat")
adultdiet <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.txt", sep = "\t", header = T)
head(adultdiet)
# Remove those that is missing DateVisitStart
sum(adultdiet$DateVisitStart == "")
# 149
adultdiet <- adultdiet[adultdiet$DateVisitStart != "",]

adultdiet$DateVisitStart <-  paste(sapply(strsplit(adultdiet$DateVisitStart, "\\/"), `[`, 3), sapply(strsplit(adultdiet$DateVisitStart, "\\/"), `[`, 1), sapply(strsplit(adultdiet$DateVisitStart, "\\/"), `[`, 2), sep ="-")

## There are some missing variables in the data that Siddhant provided, so extracting those from the original SAS data
original.adultdiet <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/ffq_grams.sas7bdat")
colnames(original.adultdiet)

dim(original.adultdiet)
# 10807   681

# original.adultdiet seems to have leading zeros in date, so removing those zeros
original.adultdiet$DateVisitStart <- gsub("-0", "-", original.adultdiet$DateVisitStart)
original.adultdiet$KEY <- paste(original.adultdiet$sjlid, original.adultdiet$DateVisitStart, sep = ":")
adultdiet$KEY <- paste(adultdiet$sjlid, adultdiet$DateVisitStart, sep = ":") # Key to search in original.adultdiet

sum(adultdiet$KEY %in% original.adultdiet$KEY) # Does match all samples Keys by date
# 6037
dim(adultdiet)
# 6037  35

original.adultdiet <- original.adultdiet[original.adultdiet$KEY %in% adultdiet$KEY,]

# keep this in the same order as adultdiet
original.adultdiet <- original.adultdiet[match(adultdiet$KEY, original.adultdiet$KEY),]


# table(original.adultdiet$KEY == adultdiet$KEY)
##########################################################
## Now adding a few more variables that Siddhant missed ##
##########################################################
adultdiet$NUTSFREQ <- as.numeric(original.adultdiet$NUTSFREQ) # 5 indicates 1.0 serving/frequency per Week
## Fish is missing; see below for details
# adultdiet$ANYFISHFREQ <- original.adultdiet$ANYFISHFREQ # 3 indicates 2 serving days last week
adultdiet$SOFTDRINKSFREQ <-  as.numeric(original.adultdiet$SOFTDRINKSFREQ) # Sugary beverage; 5 indicates 1.0 per Week; 1 Never
# adultdiet$BOLOGNAFREQ <- as.numeric(original.adultdiet$BOLOGNAFREQ) # Type LunchMeat: low-fat/turkey, reg.
adultdiet$MIXEDBEEFPORKFREQ <- as.numeric(original.adultdiet$MIXEDBEEFPORKFREQ)

adultdiet$NOTFRIEDFISHFREQ <- as.numeric(original.adultdiet$NOTFRIEDFISHFREQ)

# Processed meat: hot dogs, ham, bacon, sausage
adultdiet$HOTDOGFREQ <- as.numeric(original.adultdiet$HOTDOGFREQ)
adultdiet$BACONFREQ <- as.numeric(original.adultdiet$BACONFREQ)
adultdiet$SAUSAGEFREQ <- as.numeric(original.adultdiet$SAUSAGEFREQ)
adultdiet$G_NWHL <- as.numeric(original.adultdiet$G_NWHL)


adultdiet$AGE <- as.numeric(adultdiet$AGE) ## This age seems to be wrong; for example, SJL1287901 age


## Add DOB
adultdiet$DOB <- PHENO.ANY_SN$dob[match(adultdiet$sjlid, PHENO.ANY_SN$sjlid)]

adultdiet$AGE_at_Visit <- time_length(interval(as.Date(adultdiet$DOB), as.Date(adultdiet$DateVisitStart)), "years")

## Keep the earliest age after 18
samples.sjlife <- unique(adultdiet$sjlid)

# Only interested in those present in the phenotype data
sum(samples.sjlife%in% PHENO.ANY_SN$sjlid)
# 3574
samples.sjlife <- samples.sjlife[(samples.sjlife%in% PHENO.ANY_SN$sjlid)]

length(samples.sjlife)
# 3574


##########
## Diet ## 
##########
colnames(adultdiet)
# > colnames(adultdiet)
# [1] "mrn"                 "sjlid"               "DateVisitStart"      "visittype"           "AGE"                 "SEX"                 "DIET"                
# [8] "HEI2005_TOTAL_SCORE" "HEI2010_TOTAL_SCORE" "HEI2015_TOTAL_SCORE" "AHEI_VEGS"           "AHEI_FRUITS"         "AHEI_SUGBEVS"        "AHEI_NUTLEGS"       
# [15] "AHEI_RMEATS"         "AHEI_TRFATPCT"       "AHEI2010"            "VEGSRV"              "FRUITSRV"            "GRAINSRV"            "MEATSRV"            
# [22] "WGRAINS"             "DAIRYSRV"            "FATSRV"              "DT_SODI"             "DT_TFAT"             "DT_CARB"             "DT_PROT"            
# [29] "M_EGG"               "AV_TOT_S"            "AF_TOT_S"            "R_MEAT_S"            "A_NUT_S"             "A_BEAN_S"            "KEY"                
# [36] "NUTSFREQ"            "SOFTDRINKSFREQ"      "DOB"                 "AGE_at_Visit"   


DIET <- adultdiet

#----------------------------1. Fruits
DIET$FRUITSRV_yn <- as.numeric(ifelse(DIET$FRUITSRV >= 3, 1, 0)) # fruits; daily servings

#----------------------------2. Nuts
DIET$NUTSFREQ_yn <- as.numeric(ifelse(DIET$NUTSFREQ >= 5, 1, 0)) # NUTSFREQ 1 or more serving per week; indicated by 5 or more

#----------------------------3. Veggies
DIET$VEGSRV_yn <- as.numeric(ifelse(DIET$VEGSRV >= 3, 1, 0)) # Veggie; daily servings

#----------------------------4. Whole grains
DIET$WGRAINS_yn <- as.numeric(ifelse(DIET$WGRAINS >= 3, 1, 0)) # whole grains; daily servings

#----------------------------5. Fish 
DIET$NOTFRIEDFISHFREQ_yn <- as.numeric(ifelse(DIET$NOTFRIEDFISHFREQ >= 6, 1, 0)) # Other fish; 6 indicates 2 times per week; 1 Never

#----------------------------6. Dairy
DIET$DAIRYSRV_yn <- as.numeric(ifelse(DIET$DAIRYSRV >= 2.5, 1, 0)) # Dairy; daily servings

#----------------------------9. Unprocessed red meats
DIET$MIXEDBEEFPORKFREQ_yn <- as.numeric(ifelse(DIET$MIXEDBEEFPORKFREQ <= 5, 1, 0)) # 5 indicates 1 serving per week; 1 Never


#----------------------------10. Trans fat
DIET$DT_TFAT_cohort_median_yn <- as.numeric(ifelse(DIET$DT_TFAT <= median(DIET$DT_TFAT, na.rm = T), 1, 0)) # Fat serving less than or equal to cohort median
#----------------------------11. Sugary beverage
DIET$SOFTDRINKSFREQ_yn <- as.numeric(ifelse(DIET$SOFTDRINKSFREQ <= 5, 1, 0)) # Sugary beverage; 5 indicates 1.0 per Week; 1 Never
#----------------------------12. Sodium
DIET$DT_SODI_yn <- as.numeric(ifelse(DIET$DT_SODI <= 2000, 1, 0))

## Define Healthy diet
DIET$HEALTHY_Diet_yn <- ifelse(rowSums(DIET[,c("FRUITSRV_yn", "NUTSFREQ_yn", "VEGSRV_yn", "WGRAINS_yn", "NOTFRIEDFISHFREQ_yn", "DAIRYSRV_yn", "MIXEDBEEFPORKFREQ_yn", "DT_TFAT_cohort_median_yn", "SOFTDRINKSFREQ_yn", "DT_SODI_yn")]) >= 5, 1,0)

DIET$AGE_at_Visit[is.na(DIET$AGE_at_Visit)] <- DIET$AGE[is.na(DIET$AGE_at_Visit)]

DIET$agesurvey <- as.numeric(DIET$AGE_at_Visit)

## Extract the unique ones by keeping the earliest age after 18 years

diet_iid_dob_18 = subset(DIET, agesurvey >= 18)
# If heavy or risky drink any is 1, then Yes for RiskyHeavyDrink_yn
diet_iid_dob_18 <- diet_iid_dob_18[!is.na(diet_iid_dob_18$HEALTHY_Diet_yn), ]
diet_iid_dob_18_sorted = diet_iid_dob_18[order(diet_iid_dob_18$sjlid, diet_iid_dob_18$agesurvey, decreasing = FALSE),]
diet_iid_dob_18_uniq = diet_iid_dob_18_sorted[!duplicated(diet_iid_dob_18_sorted$sjlid),]



HEI2015_iid_dob_18 = subset(DIET, agesurvey >= 18)
# If heavy or risky drink any is 1, then Yes for RiskyHeavyDrink_yn
HEI2015_iid_dob_18 <- HEI2015_iid_dob_18[!is.na(HEI2015_iid_dob_18$HEI2015_TOTAL_SCORE), ]
HEI2015_iid_dob_18_sorted = HEI2015_iid_dob_18[order(HEI2015_iid_dob_18$sjlid, HEI2015_iid_dob_18$agesurvey, decreasing = FALSE),]
HEI2015_iid_dob_18_uniq = HEI2015_iid_dob_18_sorted[!duplicated(HEI2015_iid_dob_18_sorted$sjlid),]




sjlid <- PHENO.ANY_SN$sjlid


# Creating a dataframe with only required lifestyle variables
ALL.LIFESTYLE <- cbind.data.frame(SJLIFEID = sjlid, PhysicalActivity_yn = MET_iid_dob_18_uniq$PhysicalActivity_yn[match(sjlid, MET_iid_dob_18_uniq$SJLIFEID)])
ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey <- MET_iid_dob_18_uniq$agesurvey[match(ALL.LIFESTYLE$SJLIFEID, MET_iid_dob_18_uniq$SJLIFEID)]
ALL.LIFESTYLE$PhysicalActivity_yn[is.na(ALL.LIFESTYLE$PhysicalActivity_yn)] <- "Unknown"
ALL.LIFESTYLE$PhysicalActivity_yn <- factor(ALL.LIFESTYLE$PhysicalActivity_yn, level = c(1, 0, "Unknown")) 
table(ALL.LIFESTYLE$PhysicalActivity_yn)
# 1       0 Unknown 
# 1868    1696     837

ALL.LIFESTYLE$smoker_former_or_never_yn <- smk_iid_dob_18_uniq.2$smoker_former_or_never_yn[match(ALL.LIFESTYLE$SJLIFEID, smk_iid_dob_18_uniq.2$SJLIFEID)]
ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey <- smk_iid_dob_18_uniq.2$agesurvey[match(ALL.LIFESTYLE$SJLIFEID, smk_iid_dob_18_uniq.2$SJLIFEID)]
ALL.LIFESTYLE$smoker_former_or_never_yn[is.na(ALL.LIFESTYLE$smoker_former_or_never_yn)] <- "Unknown"
ALL.LIFESTYLE$smoker_former_or_never_yn <- factor(ALL.LIFESTYLE$smoker_former_or_never_yn, level = c(1, 0, "Unknown")) 
table(ALL.LIFESTYLE$smoker_former_or_never_yn)
# 1       0 Unknown 
# 2784     770     847 

ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn <- drk_iid_dob_18_uniq.2$NOT_RiskyHeavyDrink_yn[match(ALL.LIFESTYLE$SJLIFEID, drk_iid_dob_18_uniq.2$SJLIFEID)]
ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey <- drk_iid_dob_18_uniq.2$agesurvey[match(ALL.LIFESTYLE$SJLIFEID, drk_iid_dob_18_uniq.2$SJLIFEID)]
ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn[is.na(ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn)] <- "Unknown"
ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn <- factor(ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn, level = c(1, 0, "Unknown")) 
table(ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn)
# 1       0 Unknown 
# 2210    1306     885

ALL.LIFESTYLE$Not_obese_yn <- bmi_iid_dob_18_uniq$Not_obese_yn[match(ALL.LIFESTYLE$SJLIFEID, bmi_iid_dob_18_uniq$sjlid)]
ALL.LIFESTYLE$Not_obese_yn_agesurvey <- bmi_iid_dob_18_uniq$agebmi[match(ALL.LIFESTYLE$SJLIFEID, bmi_iid_dob_18_uniq$sjlid)]
ALL.LIFESTYLE$Not_obese_yn[is.na(ALL.LIFESTYLE$Not_obese_yn)] <- "Unknown"
ALL.LIFESTYLE$Not_obese_yn <- factor(ALL.LIFESTYLE$Not_obese_yn, level = c(1, 0, "Unknown")) 
table(ALL.LIFESTYLE$Not_obese_yn)
# 1       0 Unknown 
# 2398    1274     729

ALL.LIFESTYLE$HEALTHY_Diet_yn <- diet_iid_dob_18_uniq$HEALTHY_Diet_yn[match(ALL.LIFESTYLE$SJLIFEID, diet_iid_dob_18_uniq$sjlid)]
ALL.LIFESTYLE$HEALTHY_Diet_yn_agesurvey <- diet_iid_dob_18_uniq$agesurvey[match(ALL.LIFESTYLE$SJLIFEID, diet_iid_dob_18_uniq$sjlid)]
ALL.LIFESTYLE$HEALTHY_Diet_yn[is.na(ALL.LIFESTYLE$HEALTHY_Diet_yn)] <- "Unknown"
ALL.LIFESTYLE$HEALTHY_Diet_yn <- factor(ALL.LIFESTYLE$HEALTHY_Diet_yn, level = c(1, 0, "Unknown")) 
table(ALL.LIFESTYLE$HEALTHY_Diet_yn)
# 1       0 Unknown 
# 306    2850    1245

ALL.LIFESTYLE$HEI2015_TOTAL_SCORE <- HEI2015_iid_dob_18_uniq$HEI2015_TOTAL_SCORE[match(ALL.LIFESTYLE$SJLIFEID, HEI2015_iid_dob_18_uniq$sjlid)]
ALL.LIFESTYLE$HEI2015_TOTAL_SCORE_agesurvey <- HEI2015_iid_dob_18_uniq$agesurvey[match(ALL.LIFESTYLE$SJLIFEID, HEI2015_iid_dob_18_uniq$sjlid)]


# # Yadav's emial on 7/19/2022: In RT dose variables, you may see values of 20 cGy
# # (centigray) and 200 cGy which are scatter doses than actual doses to that body
# # region. For example, if someone was treated with 2000 cGy (or 20 Gray) of
# # radiation to pelvis, it is possible some radiation will get scattered to
# # nearby body regions such as leg or abdomen or even chest. These doses are
# # indicated as scatter doses and a value of 20 cGy indicates scatter low (or SL
# # in some cases) and 200 cGy indicates scatter high (or SH in some cases). For
# # the above example, doses to abdomen would probably be scatter high because it
# # is close to pelvis (actual region of irradiation) but doses to chest would
# # probably be scatter low. Both scatter low and high doses should not be
# # considered as RT-exposed. When you categorize RT dose variables, values <=200
# # cGy (that is, scatter doses and 0 dose) should be categorized as no exposure
# # (or 0 dose).
# 
# # ## Do not: Re-classify categorical variables based on Yadav's email above!!
# # ## maxsegrtdose
# # epsilon <- 1e-10
# # PHENO.ANY_SN$maxsegrtdose.gray <- PHENO.ANY_SN$maxsegrtdose/100
# # # Define the breaks for the categories
# # breaks <- c(-epsilon, 0, 18-epsilon, 30-epsilon, Inf)
# # # Define the labels for the categories
# # labels <- c("None", ">0-<18", ">=18-<30", ">=30")
# # # Cut the variable into categories
# # PHENO.ANY_SN$maxsegrtdose.category <- cut(PHENO.ANY_SN$maxsegrtdose.gray, breaks = breaks, labels = labels, include.lowest = TRUE)
# # levels(PHENO.ANY_SN$maxsegrtdose.category) <- c(levels(PHENO.ANY_SN$maxsegrtdose.category), "Unknown")
# # PHENO.ANY_SN$maxsegrtdose.category [is.na(PHENO.ANY_SN$maxsegrtdose)] <- "Unknown"
# # # gg <- cbind.data.frame(PHENO.ANY_SN$maxsegrtdose, PHENO.ANY_SN$maxsegrtdose.category, PHENO.ANY_SN$maxsegrtdose.gray.category)
# # 
# # 
# # ## maxabdrtdose
# # epsilon <- 1e-10
# # PHENO.ANY_SN$maxabdrtdose.gray <- PHENO.ANY_SN$maxabdrtdose/100
# # # Define the breaks for the categories
# # breaks <- c(-epsilon, 0, 30-epsilon, Inf)
# # # Define the labels for the categories
# # labels <- c("None", ">0-<30", ">=30")
# # # Cut the variable into categories
# # PHENO.ANY_SN$maxabdrtdose.category <- cut(PHENO.ANY_SN$maxabdrtdose.gray, breaks = breaks, labels = labels, include.lowest = TRUE)
# # levels(PHENO.ANY_SN$maxabdrtdose.category) <- c(levels(PHENO.ANY_SN$maxabdrtdose.category), "Unknown")
# # PHENO.ANY_SN$maxabdrtdose.category [is.na(PHENO.ANY_SN$maxabdrtdose)] <- "Unknown"
# # gg <- cbind.data.frame(PHENO.ANY_SN$maxabdrtdose, PHENO.ANY_SN$maxabdrtdose.category)
# # 
# # ## maxpelvisrtdose
# # epsilon <- 1e-10
# # PHENO.ANY_SN$maxpelvisrtdose.gray <- PHENO.ANY_SN$maxpelvisrtdose/100
# # # Define the breaks for the categories
# # breaks <- c(-epsilon, 0, 20-epsilon, Inf)
# # # Define the labels for the categories
# # labels <- c("None", ">0-<20", ">=20")
# # # Cut the variable into categories
# # PHENO.ANY_SN$maxpelvisrtdose.category <- cut(PHENO.ANY_SN$maxpelvisrtdose.gray, breaks = breaks, labels = labels, include.lowest = TRUE)
# # levels(PHENO.ANY_SN$maxpelvisrtdose.category) <- c(levels(PHENO.ANY_SN$maxpelvisrtdose.category), "Unknown")
# # PHENO.ANY_SN$maxpelvisrtdose.category [is.na(PHENO.ANY_SN$maxpelvisrtdose)] <- "Unknown"
# # 
# # ## maxchestrtdose
# # epsilon <- 1e-10
# # PHENO.ANY_SN$maxchestrtdose.gray <- PHENO.ANY_SN$maxchestrtdose/100
# # # Define the breaks for the categories
# # breaks <- c(-epsilon, 0, 20-epsilon, Inf)
# # # Define the labels for the categories
# # labels <- c("None", ">0-<20", ">=20")
# # # Cut the variable into categories
# # PHENO.ANY_SN$maxchestrtdose.category <- cut(PHENO.ANY_SN$maxchestrtdose.gray, breaks = breaks, labels = labels, include.lowest = TRUE)
# # levels(PHENO.ANY_SN$maxchestrtdose.category) <- c(levels(PHENO.ANY_SN$maxchestrtdose.category), "Unknown")
# # PHENO.ANY_SN$maxchestrtdose.category [is.na(PHENO.ANY_SN$maxchestrtdose)] <- "Unknown"
# # 
# # 
# # ## maxneckrtdose
# # epsilon <- 1e-10
# # PHENO.ANY_SN$maxneckrtdose.gray <- PHENO.ANY_SN$maxneckrtdose/100
# # # Define the breaks for the categories
# # breaks <- c(-epsilon, 0, 11-epsilon, 20-epsilon, 30-epsilon, Inf)
# # # Define the labels for the categories
# # labels <- c("None", ">0-<11", ">=11-<20", ">=20-<30", ">=30")
# # # Cut the variable into categories
# # PHENO.ANY_SN$maxneckrtdose.category <- cut(PHENO.ANY_SN$maxneckrtdose.gray, breaks = breaks, labels = labels, include.lowest = TRUE)
# # levels(PHENO.ANY_SN$maxneckrtdose.category) <- c(levels(PHENO.ANY_SN$maxneckrtdose.category), "Unknown")
# # PHENO.ANY_SN$maxneckrtdose.category [is.na(PHENO.ANY_SN$maxneckrtdose)] <- "Unknown"


## Recheck chemo and radiation variables
## read from datafreeze
chemo.2020 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat')
# chemo.2020$match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)
# chemo.2018 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20181231/Clinical Data/chemosum_dose.sas7bdat')

radiation.2020 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat')
# radiation.2018 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20181231/Clinical Data/radiation_dosimetry.sas7bdat')

# ## Also check sex and dob 2020 data
demographics <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat')
PHENO.ANY_SN$gender <- demographics$gender[match(PHENO.ANY_SN$sjlid, demographics$sjlid)]
PHENO.ANY_SN$gender <- factor(PHENO.ANY_SN$gender, levels = c("Male", "Female"))

# table(PHENO.ANY_SN$gender == PHENO.ANY_SN$gender.2) ## this looks correct
PHENO.ANY_SN$dob <- demographics$dob[match(PHENO.ANY_SN$sjlid, demographics$sjlid)]
# table(PHENO.ANY_SN$dob == PHENO.ANY_SN$dob.2) ## this looks correct

# # ## Also check sex and dob 2018 data
# demographics <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20181231/Clinical Data/demographics.sas7bdat')
# PHENO.ANY_SN$gender.2 <- demographics$gender[match(PHENO.ANY_SN$sjlid, demographics$sjlid)]
# table(PHENO.ANY_SN$gender == PHENO.ANY_SN$gender.2) ## this looks correct
# # ## Also check dob
# PHENO.ANY_SN$dob.2 <- demographics$dob[match(PHENO.ANY_SN$sjlid, demographics$sjlid)]
# table(PHENO.ANY_SN$dob == PHENO.ANY_SN$dob.2) ## this looks correct


# ## Also check sex and dob 2020 data
diagnosis <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat')
PHENO.ANY_SN$agedx <- diagnosis$agedx[match(PHENO.ANY_SN$sjlid, diagnosis$sjlid)]
# table(PHENO.ANY_SN$agedx == PHENO.ANY_SN$agedx.2) ## six were incorrect
# ## Also check diagdt
PHENO.ANY_SN$diagdt <- diagnosis$diagdt[match(PHENO.ANY_SN$sjlid, diagnosis$sjlid)]
# table(PHENO.ANY_SN$diagdt == PHENO.ANY_SN$diagdt.2) ## six were incorrect

# # ## Also check sex and dob 2018 data
# PHENO.ANY_SN$diagdt.3 <- diagnosis$diagdt[match(PHENO.ANY_SN$sjlid, diagnosis$sjlid)]
# table(PHENO.ANY_SN$diagdt.2 == PHENO.ANY_SN$diagdt.3) ## this looks correct


PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 0 & PHENO.ANY_SN$agedx < 5 ] <- "0-4"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 5 & PHENO.ANY_SN$agedx < 10 ] <- "5-9"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 10 & PHENO.ANY_SN$agedx < 15 ] <- "10-14"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 15 ] <- ">=15"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- factor(PHENO.ANY_SN$AGE_AT_DIAGNOSIS, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference


## Age at last contact
lstcondt <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data/lstcondt.sas7bdat')
PHENO.ANY_SN$agelstcontact <- lstcondt$agelstcontact[match(PHENO.ANY_SN$sjlid, lstcondt$sjlid)]
# table(PHENO.ANY_SN$agelstcontact.2 == PHENO.ANY_SN$agelstcontact) # This looks correct

## Remove unnecessary columns
PHENO.ANY_SN <- PHENO.ANY_SN[!colnames(PHENO.ANY_SN) %in% c("aa_class_dose_any_yn", "aa_class_dose_any.category", "aa_class_dose_5_yn", "cisplat_dose_any_yn", "cisplat_dose_any.category", "cisplateq_dose_5_yn", "cisplateq_dose_5.category",
  "aa_hvymtl_dose_any_yn", "aa_hvymtl_dose_any.category", "aa_hvymtl_dose_5_yn", "aa_hvymtl_dose_5.category", "anthra_jco_dose_any_yn", "anthra_jco_dose_any.category",
  "anthra_jco_dose_5_yn", "epitxn_dose_any_yn", "epitxn_dose_any.category", "epitxn_dose_5_yn",
  "aa_class_dose_any", "aa_hvymtl_dose_5", "aa_hvymtl_dose_any", "cisplat_dose_5", "cisplat_dose_any", "cisplateq_dose_5", "cisplateq_dose_any", 
  "epitxn_dose_any", "anthra_cog_dose_5", "anthra_cog_dose_any", "anthra_jco_dose_any", "carbo_dose_any", "carbo_dose_5" )]


## MaxsegRT
PHENO.ANY_SN$maxsegrtdose <- radiation.2020$maxsegrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxsegrtdose.category <- cut(PHENO.ANY_SN$maxsegrtdose, breaks = c(0, 200, 1799, 2999, max(PHENO.ANY_SN$maxsegrtdose, na.rm = T)),
                                          labels = c("None", ">0-<18", ">=18-<30", ">=30"),
                                          include.lowest = TRUE)
levels(PHENO.ANY_SN$maxsegrtdose.category) <- c(levels(PHENO.ANY_SN$maxsegrtdose.category), "Unknown")
PHENO.ANY_SN$maxsegrtdose.category [is.na(PHENO.ANY_SN$maxsegrtdose)] <- "Unknown"



## Maxchest
PHENO.ANY_SN$maxchestrtdose <- radiation.2020$maxchestrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxchestrtdose.category <- cut(PHENO.ANY_SN$maxchestrtdose, breaks = c(0, 200, 1999, max(PHENO.ANY_SN$maxchestrtdose, na.rm = T)),
                                            labels = c("None", ">0-<20", ">=20"),
                                            include.lowest = TRUE)
levels(PHENO.ANY_SN$maxchestrtdose.category) <- c(levels(PHENO.ANY_SN$maxchestrtdose.category), "Unknown")
PHENO.ANY_SN$maxchestrtdose.category [is.na(PHENO.ANY_SN$maxchestrtdose)] <- "Unknown"


## Maxneck
PHENO.ANY_SN$maxneckrtdose <- radiation.2020$maxneckrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxneckrtdose.category <- cut(PHENO.ANY_SN$maxneckrtdose, breaks = c(0, 200, 1099, 1999, 2999, max(PHENO.ANY_SN$maxneckrtdose, na.rm = T)),
                                           labels = c("None", ">0-<11", ">=11-<20", ">=20-<30", ">=30"),
                                           include.lowest = TRUE)
levels(PHENO.ANY_SN$maxneckrtdose.category) <- c(levels(PHENO.ANY_SN$maxneckrtdose.category), "Unknown")
PHENO.ANY_SN$maxneckrtdose.category [is.na(PHENO.ANY_SN$maxneckrtdose)] <- "Unknown"


## Maxabdomen
PHENO.ANY_SN$maxabdrtdose <- radiation.2020$maxabdrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxabdrtdose.category <- cut(PHENO.ANY_SN$maxabdrtdose, breaks = c(0, 200, 2999, max(PHENO.ANY_SN$maxabdrtdose, na.rm = T)),
                                          labels = c("None", ">0-<30", ">=30"),
                                          include.lowest = TRUE)
levels(PHENO.ANY_SN$maxabdrtdose.category) <- c(levels(PHENO.ANY_SN$maxabdrtdose.category), "Unknown")
PHENO.ANY_SN$maxabdrtdose.category [is.na(PHENO.ANY_SN$maxabdrtdose)] <- "Unknown"


## Maxpelvis
PHENO.ANY_SN$maxpelvisrtdose <- radiation.2020$maxpelvisrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxpelvisrtdose.category <- cut(PHENO.ANY_SN$maxpelvisrtdose, breaks = c(0, 200, 1999, max(PHENO.ANY_SN$maxpelvisrtdose, na.rm = T)),
                                             labels = c("None", ">0-<20", ">=20"),
                                             include.lowest = TRUE)
levels(PHENO.ANY_SN$maxpelvisrtdose.category) <- c(levels(PHENO.ANY_SN$maxpelvisrtdose.category), "Unknown")
PHENO.ANY_SN$maxpelvisrtdose.category [is.na(PHENO.ANY_SN$maxpelvisrtdose)] <- "Unknown"


## aa_class_dose_5 **
PHENO.ANY_SN$aa_class_dose_5 <- chemo.2020$alkylating_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)]
TERT = unname(quantile(PHENO.ANY_SN$aa_class_dose_5[PHENO.ANY_SN$aa_class_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$aa_class_dose_5.category <- cut(PHENO.ANY_SN$aa_class_dose_5, breaks = c(-Inf, 0, TERT),
                                             labels = c("None", "1st", "2nd", "3rd"),
                                             include.lowest = TRUE)
# levels(PHENO.ANY_SN$aa_class_dose_5.category) <- c(levels(PHENO.ANY_SN$aa_class_dose_5.category), "Unknown")
# PHENO.ANY_SN$aa_class_dose_5.category [is.na(PHENO.ANY_SN$aa_class_dose_5.category)] <- "Unknown"


## anthra_jco_dose_5 **
PHENO.ANY_SN$anthra_jco_dose_5 <- chemo.2020$anthracyclines_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)]
TERT = unname(quantile(PHENO.ANY_SN$anthra_jco_dose_5[PHENO.ANY_SN$anthra_jco_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$anthra_jco_dose_5.category <- cut(PHENO.ANY_SN$anthra_jco_dose_5, breaks = c(-Inf, 0, TERT),
                                               labels = c("None", "1st", "2nd", "3rd"),
                                               include.lowest = TRUE)
# levels(PHENO.ANY_SN$anthra_jco_dose_5.category) <- c(levels(PHENO.ANY_SN$anthra_jco_dose_5.category), "Unknown")
# PHENO.ANY_SN$anthra_jco_dose_5.category [is.na(PHENO.ANY_SN$anthra_jco_dose_5.category)] <- "Unknown"


## epitxn_dose_5 **
PHENO.ANY_SN$epitxn_dose_5 <- chemo.2020$epipodophyllotoxins_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)]
TERT = unname(quantile(PHENO.ANY_SN$epitxn_dose_5[PHENO.ANY_SN$epitxn_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$epitxn_dose_5.category <- cut(PHENO.ANY_SN$epitxn_dose_5, breaks = c(-Inf, 0, TERT),
                                           labels = c("None", "1st", "2nd", "3rd"),
                                           include.lowest = TRUE)
# levels(PHENO.ANY_SN$epitxn_dose_5.category) <- c(levels(PHENO.ANY_SN$epitxn_dose_5.category), "Unknown")
# PHENO.ANY_SN$epitxn_dose_5.category [is.na(PHENO.ANY_SN$epitxn_dose_5.category)] <- "Unknown"



save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11_exactly_as_Qin.RDATA")

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")