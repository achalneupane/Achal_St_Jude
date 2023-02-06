
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

MET_iid_dob_18_uniq$PhysicalActivity_yn <- ifelse(MET_iid_dob_18_uniq$CDC_PA == 1, "Yes", "No")



## 2.-------------- Smoker status
# Get evsm, smnow, smnvr, cigmo, cigd, smyr
lifestyle$evsm <- adult_habbits$evsm[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$smnow <- adult_habbits$smnow[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$cigmo <- adult_habbits$cigmo[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$cigd <- adult_habbits$cigd[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$smyr <- adult_habbits$smyr[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$smnvr <- NA


# if evsm=2 then do ; smnow=2 ; cigmo=2 ; cigd=0 ; smyr=0 ; end ;
lifestyle$smnow[which(lifestyle$evsm == 2)] <- 2
lifestyle$cigmo[which(lifestyle$evsm == 2)] <- 2
lifestyle$cigd[which(lifestyle$evsm == 2)] <- 0
lifestyle$smyr[which(lifestyle$evsm == 2)] <- 0

# if smnvr=1 then do ; evsm=2 ; smnow=2 ; cigmo=2 ; cigd=0 ; smyr=0 ; end ;
lifestyle$evsm[which(lifestyle$smnvr == 1)] <- 2
lifestyle$smnow[which(lifestyle$smnvr == 1)] <- 2
lifestyle$cigmo[which(lifestyle$smnvr == 1)] <- 2
lifestyle$cigd[which(lifestyle$smnvr == 1)] <- 0
lifestyle$smyr[which(lifestyle$smnvr == 1)] <- 0

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

smk_iid_dob_18_uniq$current_yn <- ifelse(smk_iid_dob_18_uniq$smkStat != 2, "No", "Yes")
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


# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")