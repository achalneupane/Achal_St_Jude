
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
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/3_PRS_scores_categories_v2.RDATA")
#############################
## Add Lifestyle variables ##
#############################
## For each samples, get habits immediately after 18 years of age in agesurvey
# adlthabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.sas7bdat")
adlthabits <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.txt", header = T, sep ="\t")
head(adlthabits)
## Fix data format
adlthabits$datecomp <- paste(sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 3), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 1), sapply(strsplit(adlthabits$datecomp, "\\/"), `[`, 2), sep = "-")

# remove duplicated rows
adlthabits <- distinct(adlthabits)
## Get DOB
adlthabits$DOB <- PHENO.ANY_SN$dob [match(adlthabits$SJLIFEID, PHENO.ANY_SN$sjlid)]
adlthabits <- adlthabits[!is.na(adlthabits$DOB),]
# change the format of dates YYYY-MM-DD
adlthabits$agesurvey <- time_length(interval(as.Date(adlthabits$DOB), as.Date(adlthabits$datecomp)), "years")
# adlthabits$agesurvey <- floor(adlthabits$agesurvey_diff_of_dob_datecomp)

samples.sjlife <- unique(adlthabits$SJLIFEID)
# Only interested in those present in the phenotype data
sum(samples.sjlife%in% PHENO.ANY_SN$sjlid)
# 3571
samples.sjlife <- samples.sjlife[(samples.sjlife%in% PHENO.ANY_SN$sjlid)]

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


sum(duplicated(lifestyle$SJLIFEID))
# 2
lifestyle$SJLIFEID[duplicated(lifestyle$SJLIFEID)]
# "SJL1080201" "SJL5359215"
## Remove duplicate row
lifestyle <- lifestyle[!duplicated(lifestyle$SJLIFEID),]

##################################
## Recode categorical variables ##
##################################
lifestyle$relation[lifestyle$relation == 1] <- "Self"
lifestyle$relation[lifestyle$relation == 2] <- "Parent"
lifestyle$relation[lifestyle$relation == 3] <- "Other"

## Recode smoker
lifestyle$smoker[lifestyle$smoker == 1] <- "Past"
lifestyle$smoker[lifestyle$smoker == 2] <- "Current"
lifestyle$smoker[lifestyle$smoker == 3] <- "Never"
lifestyle$smoker_former_or_never_yn <- as.numeric(ifelse(lifestyle$smoker != "Current", 1, 0))
lifestyle$smoker_never_yn <- as.numeric(ifelse(lifestyle$smoker == "Never", 1, 0))


## Recode 1/2 or 0/1 to 0 (N) and 1 (Y)
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

###############
## Adult BMI ##
###############
# This BMI has more samples, so using this for BMI only
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

##############################
## More lifestyle variables ##
##############################

# Keep the earliest age after 18 years, same as in lifestyle
# adultbmi <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.sas7bdat")
adultbmi <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.txt", sep = "\t", header = T)
head(adultbmi)
# Remove those that is missing DateVisitStart
sum(adultbmi$DateVisitStart == "")
# 149
adultbmi <- adultbmi[adultbmi$DateVisitStart != "",]

adultbmi$DateVisitStart <-  paste(sapply(strsplit(adultbmi$DateVisitStart, "\\/"), `[`, 3), sapply(strsplit(adultbmi$DateVisitStart, "\\/"), `[`, 1), sapply(strsplit(adultbmi$DateVisitStart, "\\/"), `[`, 2), sep ="-")



## There are some missing variables in the data that Siddhant provided, so extracting those from the original SAS data
original.adultbmi <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/ffq_grams.sas7bdat")
colnames(original.adultbmi)

dim(original.adultbmi)
# utils::View(original.adultbmi)
# original.adultbmi seems to have leading zeros in date, so removing those zeors
original.adultbmi$DateVisitStart <- gsub("-0", "-", original.adultbmi$DateVisitStart)
original.adultbmi$KEY <- paste(original.adultbmi$sjlid, original.adultbmi$DateVisitStart, sep = ":")
adultbmi$KEY <- paste(adultbmi$sjlid, adultbmi$DateVisitStart, sep = ":") # Key to search in original.adultbmi

sum(adultbmi$KEY %in% original.adultbmi$KEY) # Does match all samples Keys by date
# 6037
dim(adultbmi)
# 6037  35

original.adultbmi <- original.adultbmi[original.adultbmi$KEY %in% adultbmi$KEY,]

# keep this in the same order as adultbmi
original.adultbmi <- original.adultbmi[match(adultbmi$KEY, original.adultbmi$KEY),]

# ## check few variables for sanity check between the two datasets
# table(original.adultbmi$VEGSRV == adultbmi$VEGSRV)
# table(original.adultbmi$WGRAINS == adultbmi$WGRAINS) 
# table(original.adultbmi$M_EGG == adultbmi$M_EGG)
# 
# original.adultbmi$sjlid[which(original.adultbmi$M_EGG != adultbmi$M_EGG)]
# 
# ## Variables for SJL5248101 seem to be inconsistent between the Siddhant's adultbmi data and SJLIFE data
# adultbmi$M_EGG[grepl("SJL5248101", adultbmi$sjlid)]
# original.adultbmi$M_EGG[grepl("SJL5248101", original.adultbmi$sjlid)]

table(original.adultbmi$KEY == adultbmi$KEY)
##########################################################
## Now adding a few more variables that Siddhant missed ##
##########################################################
adultbmi$NUTSFREQ <- as.numeric(original.adultbmi$NUTSFREQ) # 5 indicates 1.0 serving/frequency per Week
## Fish is missing; see below for details
# adultbmi$ANYFISHFREQ <- original.adultbmi$ANYFISHFREQ # 3 indicates 2 serving days last week
adultbmi$SOFTDRINKSFREQ <-  as.numeric(original.adultbmi$SOFTDRINKSFREQ) # Sugary beverage; 5 indicates 1.0 per Week; 1 Never
# adultbmi$BOLOGNAFREQ <- as.numeric(original.adultbmi$BOLOGNAFREQ) # Type LunchMeat: low-fat/turkey, reg.
adultbmi$MIXEDBEEFPORKFREQ <- as.numeric(original.adultbmi$MIXEDBEEFPORKFREQ)

adultbmi$NOTFRIEDFISHFREQ <- as.numeric(original.adultbmi$NOTFRIEDFISHFREQ)

# Processed meat: hot dogs, ham, bacon, sausage
adultbmi$HOTDOGFREQ <- as.numeric(original.adultbmi$HOTDOGFREQ)
adultbmi$BACONFREQ <- as.numeric(original.adultbmi$BACONFREQ)
adultbmi$SAUSAGEFREQ <- as.numeric(original.adultbmi$SAUSAGEFREQ)
adultbmi$G_NWHL <- as.numeric(original.adultbmi$G_NWHL)


adultbmi$AGE <- as.numeric(adultbmi$AGE) ## This age seems to be wrong; for example, SJL1287901 age


## Add DOB
adultbmi$DOB <- PHENO.ANY_SN$dob[match(adultbmi$sjlid, PHENO.ANY_SN$sjlid)]



# adultbmi$AGE_at_Visit <- floor(time_length(interval(as.Date(adultbmi$DOB), as.Date(adultbmi$DateVisitStart)), "years"))
# table(adultbmi$AGE_at_Visit == adultbmi$AGE) # 125 $AGE seem to be wrong, so calculating the Age as AGE_at_Visit below
# WRONG.AGE <- adultbmi[which(adultbmi$AGE_at_Visit != adultbmi$AGE),]

adultbmi$AGE_at_Visit <- time_length(interval(as.Date(adultbmi$DOB), as.Date(adultbmi$DateVisitStart)), "years")

## Keep the earliest age after 18
samples.sjlife <- unique(adultbmi$sjlid)

# Only interested in those present in the phenotype data
sum(samples.sjlife%in% PHENO.ANY_SN$sjlid)
# 3574
samples.sjlife <- samples.sjlife[(samples.sjlife%in% PHENO.ANY_SN$sjlid)]

length(samples.sjlife)
# 3574


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
nrow(BMI)
# 3543
#######################
## Physical Activity ##
#######################

# Use variable 'wtlt': I do activities to increase muscle strength, such as lifting weights or aerobics, once a week or more 
lifestyle$PhysicalActivity_yn <- as.numeric(ifelse(lifestyle$wtlt == 1, 1, 0))

#############
## Obesity ##
#############

# BMI$Not_obese_yn <- as.numeric(ifelse(BMI$BMI < 30, 1, 0))

# BMI file Yadav used has more samples
bmi_iid_dob_18_uniq$Not_obese_yn <- as.numeric(ifelse(bmi_iid_dob_18_uniq$BMIadj < 30, 1, 0))
#############
## Alcohol ##
#############
# If heavy or risky drink any is 1, then Yes for RiskyHeavyDrink_yn
lifestyle$NOT_RiskyHeavyDrink_yn <- as.numeric(ifelse (rowSums(lifestyle[c("heavydrink", "riskydrink")])==0, 1, 0))


##########
## Diet ## 
##########

## Fish is missing
colnames(adultbmi)
# > colnames(adultbmi)
# [1] "mrn"                 "sjlid"               "DateVisitStart"      "visittype"           "AGE"                 "SEX"                 "BMI"                
# [8] "HEI2005_TOTAL_SCORE" "HEI2010_TOTAL_SCORE" "HEI2015_TOTAL_SCORE" "AHEI_VEGS"           "AHEI_FRUITS"         "AHEI_SUGBEVS"        "AHEI_NUTLEGS"       
# [15] "AHEI_RMEATS"         "AHEI_TRFATPCT"       "AHEI2010"            "VEGSRV"              "FRUITSRV"            "GRAINSRV"            "MEATSRV"            
# [22] "WGRAINS"             "DAIRYSRV"            "FATSRV"              "DT_SODI"             "DT_TFAT"             "DT_CARB"             "DT_PROT"            
# [29] "M_EGG"               "AV_TOT_S"            "AF_TOT_S"            "R_MEAT_S"            "A_NUT_S"             "A_BEAN_S"            "KEY"                
# [36] "NUTSFREQ"            "SOFTDRINKSFREQ"      "DOB"                 "AGE_at_Visit"   

#----------------------------1. Fruits
BMI$FRUITSRV_yn <- as.numeric(ifelse(BMI$FRUITSRV >= 3, 1, 0)) # fruits; daily servings

#----------------------------2. Nuts
BMI$NUTSFREQ_yn <- as.numeric(ifelse(BMI$NUTSFREQ >= 5, 1, 0)) # NUTSFREQ 1 or more serving per week; indicated by 5 or more

#----------------------------3. Veggies
BMI$VEGSRV_yn <- as.numeric(ifelse(BMI$VEGSRV >= 3, 1, 0)) # Veggie; daily servings

#----------------------------4. Whole grains
BMI$WGRAINS_yn <- as.numeric(ifelse(BMI$WGRAINS >= 3, 1, 0)) # whole grains; daily servings

#----------------------------5. Fish 
## There is this variable for Fish in dictionary file for child only, but
#ffq_child.sas7bdat data doesn't have it (Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE
#Data Freeze\2 Final Data SJLIFE\20200430\Clinical Data) BMI$ANYFISHFREQ <-
#as.numeric(ifelse(BMI$ANYFISHFREQ >= 3, 1, 0)) # Any fish; 3 indicates 2 days last week; 1 Not eaten in the last week

BMI$NOTFRIEDFISHFREQ_yn <- as.numeric(ifelse(BMI$NOTFRIEDFISHFREQ >= 6, 1, 0)) # Other fish; 6 indicates 2 times per week; 1 Never

#----------------------------6. Dairy
BMI$DAIRYSRV_yn <- as.numeric(ifelse(BMI$DAIRYSRV >= 2.5, 1, 0)) # Dairy; daily servings

# #----------------------------7. Grains (using this as refined grains)
# BMI$G_NWHL_yn <- as.numeric(ifelse(BMI$G_NWHL <= 1.5, 1, 0)) # Not sure about this; frequency is not clear
# BMI$GRAINSRV <- as.numeric(ifelse(BMI$GRAINSRV <= 1.5, 1, 0)) # GRAINSRV; daily servings

# #----------------------------8. Processed meats
# doi: 10.1080/01635580701684872
# Processed meat includes bacon, ham (raw, smoked or cooked), heated sausages
# like hot-dogs (frankfurters), raw sausages (like salami), bologna, blood
# sausage (UK: black pudding), liver pâté (or liverwurst) and other pâtés and
# spread meat, luncheon meat and other cold cuts, canned meat, and corned beef

# No processed meats; also frequency is not available for week
# BMI$R_MEAT_S_yn <- as.numeric(ifelse(BMI$R_MEAT_S <= 1.5, 1, 0)) # Red meat servings daily
# BMI$BOLOGNAFREQ_yn <- as.numeric(ifelse(BMI$BOLOGNAFREQ <= 5)) # Lunch Meats, 5 indicates 1 serving per week; 1 Never

# Processed meat" hotdogs, bacon, sausage

# # Processed meat: hot dogs, ham, bacon, sausage;  5 indicates 1 serving per week; 1 Never
# BMI$HOTDOGFREQ_yn <- as.numeric(ifelse(BMI$HOTDOGFREQ <= 5, 1, 0))
# BMI$BACONFREQ_yn <- as.numeric(ifelse(BMI$BACONFREQ <= 5, 1, 0)) # 5 indicates 1 serving per week; 1 Never
# BMI$SAUSAGEFREQ_yn <- as.numeric(ifelse(BMI$SAUSAGEFREQ <= 5, 1, 0)) # 5 indicates 1 serving per week; 1 Never
# BMI$processed_meats_yn <- as.numeric(ifelse(rowSums(cbind(BMI$HOTDOGFREQ, BMI$BACONFREQ, BMI$SAUSAGEFREQ)) <= 15, 1, 0)) # 5 indicates 1 serving per week; 1 Never

#----------------------------9. Unprocessed red meats
BMI$MIXEDBEEFPORKFREQ_yn <- as.numeric(ifelse(BMI$MIXEDBEEFPORKFREQ <= 5, 1, 0)) # 5 indicates 1 serving per week; 1 Never


#----------------------------10. Trans fat
BMI$DT_TFAT_cohort_median_yn <- as.numeric(ifelse(BMI$DT_TFAT <= median(BMI$DT_TFAT, na.rm = T), 1, 0)) # Fat serving less than or equal to cohort median
#----------------------------11. Sugary beverage
BMI$SOFTDRINKSFREQ_yn <- as.numeric(ifelse(BMI$SOFTDRINKSFREQ <= 5, 1, 0)) # Sugary beverage; 5 indicates 1.0 per Week; 1 Never
#----------------------------12. Sodium
BMI$DT_SODI_yn <- as.numeric(ifelse(BMI$DT_SODI <= 2000, 1, 0))

BMI$HEALTHY_Diet_yn <- ifelse(rowSums(BMI[,c("FRUITSRV_yn", "NUTSFREQ_yn", "VEGSRV_yn", "WGRAINS_yn", "NOTFRIEDFISHFREQ_yn", "DAIRYSRV_yn", "MIXEDBEEFPORKFREQ_yn", "DT_TFAT_cohort_median_yn", "SOFTDRINKSFREQ_yn", "DT_SODI_yn")]) >= 5, 1,0)

 
######################################
## Merge BMI and Lifestyle datasets ##
######################################

ALL.LIFESTYLE <- merge(lifestyle, BMI, by.x = "SJLIFEID", by.y = "sjlid", all = T)

## Only Keep the ones that are needed
ALL.LIFESTYLE <- ALL.LIFESTYLE[c("agesurvey", "SJLIFEID", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")]

## Merge BMI from Yadav
ALL.LIFESTYLE <- merge(ALL.LIFESTYLE, bmi_iid_dob_18_uniq, by.x = "SJLIFEID", by.y = "sjlid", all = T)


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA")
