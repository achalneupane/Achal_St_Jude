
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

##########################################################################################################
## Recode Physical activity and Smoker based on what Kiri Ness' suggested (Using SAS code leanmass.sas) ##
##########################################################################################################

## 1.-------------- Physical activity
sum(is.na(lifestyle$vpa10))
# get additional adult health habits
adult_habbits <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/adult_healthhabits.sas7bdat')

# Get vpadays, vpamin, mpadays, mpamin
lifestyle$vpadays <- adult_habbits$vpadays[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$vpamin <- adult_habbits$vpamin[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$mpadays <- adult_habbits$mpadays[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$mpamin <- adult_habbits$mpamin[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]

# if vpa10=1 and vpadays=. then do vpadays=1; end;
lifestyle$vpadays[which(lifestyle$vpa10 == 1 & is.na(lifestyle$vpadays))] <- 1

# if vpa10=1 and vpamin=. then do vpamin=10; end;
lifestyle$vpamin[which(lifestyle$vpa10 == 1 & is.na(lifestyle$vpamin))] <- 10
cc <- lifestyle[which(lifestyle$vpa10 == 1 & is.na(lifestyle$vpamin)),]

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

# ## check few variables for sanity check between the two datasets
# table(original.adultdiet$VEGSRV == adultdiet$VEGSRV)
# table(original.adultdiet$WGRAINS == adultdiet$WGRAINS) 
# table(original.adultdiet$M_EGG == adultdiet$M_EGG)
# 
# original.adultdiet$sjlid[which(original.adultdiet$M_EGG != adultdiet$M_EGG)]
# 
# ## Variables for SJL5248101 seem to be inconsistent between the Siddhant's adultdiet data and SJLIFE data
# adultdiet$M_EGG[grepl("SJL5248101", adultdiet$sjlid)]
# original.adultdiet$M_EGG[grepl("SJL5248101", original.adultdiet$sjlid)]

table(original.adultdiet$KEY == adultdiet$KEY)
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


DIET <- {}
for (i in 1:length(samples.sjlife)){
  print(paste0("Doing ", i))
  dat <- adultdiet[adultdiet$sjlid == samples.sjlife[i],]
  if (max(dat$AGE_at_Visit, na.rm = T) >= 18){
    print("YES")
    dat2 <- dat[dat$AGE_at_Visit >= 18,]
    DIET.tmp <- dat2[which(dat2$AGE_at_Visit == min(dat2$AGE_at_Visit, na.rm = T)),] # Keep the earliest age after 18 years
  }
  DIET <- rbind.data.frame(DIET, DIET.tmp)
}

save.DIET <- DIET

DIET <- distinct(DIET)
sum(duplicated(DIET$sjlid))
# 0
DIET$sjlid[duplicated(DIET$sjlid)]
nrow(DIET)
# 3543
#######################
## Physical Activity ##
#######################

# Use variable 'wtlt': I do activities to increase muscle strength, such as lifting weights or aerobics, once a week or more 
lifestyle$PhysicalActivity_yn <- as.numeric(ifelse(lifestyle$wtlt == 1, 1, 0))

#############
## Obesity ##
#############

# DIET$Not_obese_yn <- as.numeric(ifelse(DIET$DIET < 30, 1, 0))

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
colnames(adultdiet)
# > colnames(adultdiet)
# [1] "mrn"                 "sjlid"               "DateVisitStart"      "visittype"           "AGE"                 "SEX"                 "DIET"                
# [8] "HEI2005_TOTAL_SCORE" "HEI2010_TOTAL_SCORE" "HEI2015_TOTAL_SCORE" "AHEI_VEGS"           "AHEI_FRUITS"         "AHEI_SUGBEVS"        "AHEI_NUTLEGS"       
# [15] "AHEI_RMEATS"         "AHEI_TRFATPCT"       "AHEI2010"            "VEGSRV"              "FRUITSRV"            "GRAINSRV"            "MEATSRV"            
# [22] "WGRAINS"             "DAIRYSRV"            "FATSRV"              "DT_SODI"             "DT_TFAT"             "DT_CARB"             "DT_PROT"            
# [29] "M_EGG"               "AV_TOT_S"            "AF_TOT_S"            "R_MEAT_S"            "A_NUT_S"             "A_BEAN_S"            "KEY"                
# [36] "NUTSFREQ"            "SOFTDRINKSFREQ"      "DOB"                 "AGE_at_Visit"   

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

 
#######################################
## Merge DIET and Lifestyle datasets ##
#######################################
ALL.LIFESTYLE <- merge(lifestyle, DIET, by.x = "SJLIFEID", by.y = "sjlid", all = T)

## Only Keep the ones that are needed
ALL.LIFESTYLE <- ALL.LIFESTYLE[c("agesurvey", "SJLIFEID", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE", "smoker_never_yn", "smoker_former_or_never_yn", "PhysicalActivity_yn", "NOT_RiskyHeavyDrink_yn", "HEALTHY_Diet_yn")]


##########################
## Merge BMI from Yadav ##
##########################
ALL.LIFESTYLE <- merge(ALL.LIFESTYLE, bmi_iid_dob_18_uniq, by.x = "SJLIFEID", by.y = "sjlid", all = T)


# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA")

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA")