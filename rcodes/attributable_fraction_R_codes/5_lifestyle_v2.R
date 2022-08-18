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

###############
## Adult BMI ##
###############
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
adultbmi$ANYFISHFREQ <- original.adultbmi$ANYFISHFREQ # 3 indicates 2 serving days last week
adultbmi$SOFTDRINKSFREQ <-  as.numeric(original.adultbmi$SOFTDRINKSFREQ) # Sugary beverage; 5 indicates 1.0 per Week; 1 Neve

adultbmi$AGE <- as.numeric(adultbmi$AGE) ## This age seems to be wrong; for example, SJL1287901 age


## Keep only those BMI data only for the samples in phenotype 
adultbmi <- adultbmi[adultbmi$sjlid %in% PHENO.ANY_SN$sjlid,]; dim(adultbmi)

## Add DOB
adultbmi$DOB <- PHENO.ANY_SN$dob[match(adultbmi$sjlid, PHENO.ANY_SN$sjlid)]



# adultbmi$AGE_at_Visit <- floor(time_length(interval(as.Date(adultbmi$DOB), as.Date(adultbmi$DateVisitStart)), "years"))
# table(adultbmi$AGE_at_Visit == adultbmi$AGE) # 125 $AGE seem to be wrong, so calculating the Age as AGE_at_Visit below
# WRONG.AGE <- adultbmi[which(adultbmi$AGE_at_Visit != adultbmi$AGE),]

adultbmi$AGE_at_Visit <- time_length(interval(as.Date(adultbmi$DOB), as.Date(adultbmi$DateVisitStart)), "years")

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

# Use variable 'nopa': I rarely or never do any physical activities 
lifestyle$PhysicalActivity_YN <- factor(ifelse(lifestyle$nopa, "Y", "N"))

#############
## Obesity ##
#############

BMI$Obesity_YN <- factor(ifelse(BMI$BMI < 30, "N", "Y"))

#############
## Alcohol ##
#############
# If heavy or risky drink any is 1, then Yes for RiskyHeavyDrink_YN
lifestyle$RiskyHeavyDrink_YN <- factor(ifelse (rowSums(lifestyle[c("heavydrink", "riskydrink")])>=1, "Y", "N"))


##########
## Diet ##
##########
colnames(adultbmi)
# VEGSRV"               
# [19] "FRUITSRV"              "GRAINSRV"              "MEATSRV"               "WGRAINS"               "DAIRYSRV"              "FATSRV"               
# [25] "DT_SODI"               "DT_TFAT"               "DT_CARB"               "DT_PROT"               "M_EGG"                 "AV_TOT_S"             
# [31] "AF_TOT_S"              "R_MEAT_S"              "A_NUT_S"               "A_BEAN_S"     "NUTSFREQ_YN"
adultbmi$FRUITSRV_YN <- factor(ifelse(adultbmi$FRUITSRV >= 3, 1, 0)) # fruits
adultbmi$NUTSFREQ_YN <- factor(ifelse(adultbmi$NUTSFREQ >= 1, 1, 0))
adultbmi$VEGSRV_YN <- factor(ifelse(adultbmi$VEGSRV >= 3, 1, 0)) # Veggie 
adultbmi$WGRAINS_YN <- factor(ifelse(adultbmi$WGRAINS >= 3, 1, 0)) # whole grains
adultbmi$ANYFISHFREQ <- factor(ifelse(adultbmi$ANYFISHFREQ >= 3, 1, 0)) # Any fish; 3 indicates 2 days last week; 1 Not eaten in the last week
adultbmi$DAIRYSRV_YN <- factor(ifelse(adultbmi$DAIRYSRV >= 2.5, 1, 0)) # Dairy
adultbmi$GRAINSRV_YN <- factor(ifelse(adultbmi$GRAINSRV <= 1.5, 1, 0)) # GRAINSRV
# No processed meats
adultbmi$R_MEAT_S_YN <- factor(ifelse(adultbmi$R_MEAT_S <= 1.5, 1, 0)) # Red meat
adultbmi$DT_TFAT_cohort_median_YN <- factor(ifelse(adultbmi$DT_TFAT <= median(adultbmi$DT_TFAT, na.rm = T), 1, 0)) # Fat serving less than or equal to cohort median
adultbmi$SOFTDRINKSFREQ_YN <- factor(ifelse(adultbmi$SOFTDRINKSFREQ <= 1.5, 1, 0)) # Sugary beverage; 5 indicates 1.0 per Week; 1 Never
adultbmi$DT_SODI_YN <- factor(ifelse(adultbmi$DT_SODI <= 2000, 1, 0))

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

save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA")
