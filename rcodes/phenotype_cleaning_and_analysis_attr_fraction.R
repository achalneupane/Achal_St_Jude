

###########################################################################################
## On 05/26/2022; we received Phenotype for Email subject: 'Attribution fraction for SN' ##     
###########################################################################################

##################
##################
## CLINICAL SET ##
##################
##################

#####################
## wgspop.sas7bdat ##
#####################
# WORKDIR:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE")

# ## Recode any values other than 0-5 as NAs
# data.df$grade[!grepl("0|1|2|3|4|5", data.df$grade)] <- NA
# 
# ## Omit NA grades
# # data.df <- data.df[!is.na(data.df$grade),]
# 
# ## Clean condition stings
# data.df$condition <- gsub("_$", "", (gsub("_+", "_",  gsub("[^A-Za-z0-9]", "_", data.df$condition))))
# dim(data.df)

wgspop <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/wgspop.sas7bdat")
head(wgspop)


demog <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/demog.sas7bdat")
head(demog)
demog <- demog[,c("MRN", "dob", "gender", "race", "ethnic", "agedx", "agelstcontact")]

## Add SJLIFE ID
clinical.dat <- cbind.data.frame(wgspop[match(demog$MRN, wgspop$MRN), c("sjlid")], demog)

wgsdiag <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/wgsdiag.sas7bdat")
head(wgsdiag)
# table(wgsdiag$diaggrp)
clinical.dat <- cbind.data.frame(clinical.dat, wgsdiag[match(clinical.dat$MRN, wgsdiag$MRN), c("diagdt", "diaggrp")])

###############
## Radiation ##
###############

## Get diagrp and diagdt
radiation <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/radiation.sas7bdat")
head(radiation)

## Add all from radiation
clinical.dat <- cbind.data.frame(clinical.dat, radiation[match(clinical.dat$MRN, radiation$MRN), ])


#############################
## Adult habits/ Lifestyle ##
#############################
## For each samples, get habits immediately after 18 years of age in agesurvey

# adlthabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.sas7bdat")
adlthabits <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.txt", header = T, sep ="\t")
head(adlthabits)
# remove duplicated rows
adlthabits <- distinct(adlthabits)
adlthabits$agesurvey <- floor(adlthabits$agesurvey_diff_of_dob_datecomp)

samples.sjlife <- unique(adlthabits$SJLIFEID)

lifestyle <- {}
for (i in 1:length(samples.sjlife)){
  print(paste0("Doing ", i))
  dat <- adlthabits[adlthabits$SJLIFEID == samples.sjlife[i],]
  if (max(dat$agesurvey) >= 18){
    print("YES")
    dat2 <- dat[dat$agesurvey >= 18,]
    lifestyle.tmp <- dat2[which(dat2$agesurvey == min(dat2$agesurvey)),]
    # } else {
    #   print("NO")
    #   lifestyle.tmp <-  dat[which(dat$agesurvey == max(dat$agesurvey)),]
    #   lifestyle.tmp[9:ncol(lifestyle.tmp)] <- NA
  }
  lifestyle <- rbind.data.frame(lifestyle, lifestyle.tmp)
}

sum(duplicated(lifestyle$SJLIFEID))
lifestyle$SJLIFEID[duplicated(lifestyle$SJLIFEID)]
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
lifestyle$smoker <- ifelse(lifestyle$smoker == "Current", "Y", "N")

## Recode to Y/N
lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 1 ] <- "Y"
lifestyle[grepl("nopa|ltpa", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa", colnames(lifestyle))] == 2 ] <- "N"
lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 0 ] <- "N"

# change the format of dates YYYY-MM-DD
lifestyle$datecomp <- gsub("\\/", "-", lifestyle$datecomp)
lifestyle$datecomp <- paste(sapply(strsplit(lifestyle$datecomp, "-"), `[`, 3), sapply(strsplit(lifestyle$datecomp, "-"), `[`, 1), sapply(strsplit(lifestyle$datecomp, "-"), `[`, 2), sep ="-")

lifestyle$dob <- gsub("\\/", "-", lifestyle$dob)
lifestyle$dob <- paste(sapply(strsplit(lifestyle$dob, "-"), `[`, 3), sapply(strsplit(lifestyle$dob, "-"), `[`, 1), sapply(strsplit(lifestyle$dob, "-"), `[`, 2), sep ="-")

#######################
## Adolescent habits ##
#######################

adolhabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adolhabits.sas7bdat")
head(adolhabits)

###############
## Adult BMI ##
###############

adultbmi <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.sas7bdat")
head(adultbmi)

## Add BMI and Nutrition to Lifestyle. Here, I am extracting the the BMI and Nutrition for the same age for each sample 
lifestyle$BMI_KEY <- paste(lifestyle$sjlid, lifestyle$agesurvey, sep = ":")

length(unique(adultbmi$sjlid))
# 3640
adultbmi$BMI_KEY <- paste(adultbmi$sjlid, adultbmi$AGE, sep = ":")
sum(lifestyle$BMI_KEY %in% adultbmi$BMI_KEY)
# 2964 
## samples that did not match by corresponding age
cc <- lifestyle[!lifestyle$BMI_KEY %in% adultbmi$BMI_KEY,]

lifestyle <- cbind.data.frame(lifestyle, adultbmi[match(lifestyle$BMI_KEY, adultbmi$BMI_KEY), c("BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")])
##########
## Drug ##
##########

drug <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/drug.sas7bdat")
head(drug)


## Add drug to clinical data
clinical.dat <- cbind.data.frame(clinical.dat, drug[match(clinical.dat$MRN, drug$MRN), grep("dose_any",colnames(drug))])

colnames(clinical.dat)[grepl("_yn|anyrt_|AnyRT", colnames(clinical.dat))]
# [1] "AnyRT"            "anyrt_prim"       "anyrt_5"          "anyrt_10"         "brainorheadrt_yn" "brainrt_yn"       "chestrt_yn"      
# [8] "neckrt_yn"        "pelvisrt_yn"      "abdomenrt_yn" 

# anyrt: 1Y, 0N;  brainorheadrt_yn : 1Y, 2N; brainrt_yn: 1Y, 2N; chestrt_yn: 1Y, 2N; neckrt_yn: 1Y, 2N; pelvisrt_yn: 1Y, 2N; abdomenrt_yn: 1Y, 2N

clinical.dat[grepl("_yn|anyrt_|AnyRT", colnames(clinical.dat))][clinical.dat[grepl("_yn|anyrt_|AnyRT", colnames(clinical.dat))] == 1 ] <- "Y"
clinical.dat[grepl("_yn", colnames(clinical.dat))][clinical.dat[grepl("_yn", colnames(clinical.dat))] == 2 ] <- "N"
clinical.dat[grepl("anyrt_|AnyRT", colnames(clinical.dat))][clinical.dat[grepl("anyrt_|AnyRT", colnames(clinical.dat))] == 0 ] <- "N"


#########################
## Subsequent Neoplasm ##
#########################

subneo <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/subneo.sas7bdat")
head(subneo)
table(subneo$diaggrp)
# add DOB
subneo$DOB <- demog$dob[match(subneo$MRN, demog$MRN)]

library(lubridate)
subneo$AGE.exact.ANY_SN <- time_length(interval(as.Date(subneo$DOB), as.Date(subneo$gradedt)), "years")
subneo$AGE.ANY_SN <- floor(subneo$AGE.exact.ANY_SN)
############
## Any SNs 
############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

# ## Add all samples
# ANY_SNs <- cbind.data.frame(wgspop[,c("MRN", "sjlid")], ANY_SNs[match(wgspop$MRN, ANY_SNs$MRN), ])
# ANY_SNs <- ANY_SNs[-c(3,4)]
# ###

PHENO.ANY_SN <- clinical.dat
PHENO.ANY_SN$ANY_SN <- ifelse(PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, "Y", "N")
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ANY_SNs[match(PHENO.ANY_SN$sjlid, ANY_SNs$sjli), c("gradedt", "AGE.ANY_SN")])

## Merge lifestyle
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, lifestyle[match(PHENO.ANY_SN$sjlid, lifestyle$SJLIFEID),c("datecomp", "agesurvey", "relation", "smoker", "nopa", "ltpa", "drk5", "bingedrink", "heavydrink", "riskydrink")])

# If ANY_SN is Yes and the survey date of lifestyle (datecomp) is later than the SN grade
# date (gradedt), then the lifestyle variables for those samples will be
# irrelevant
PHENO.ANY_SN[which(PHENO.ANY_SN[PHENO.ANY_SN$ANY_SN == "Y", "datecomp"] >  PHENO.ANY_SN[PHENO.ANY_SN$ANY_SN == "Y", "gradedt"]), c("smoker", "nopa", "ltpa", "drk5", "bingedrink", "heavydrink", "riskydrink")] <- NA

# Now adding BMI and nutrition; extracting the BMI values immediately before or on the date of SN gradedt
adultbmi.wanted <- {}
for (i in 1:length(PHENO.ANY_SN$sjlid)){
  tmp.adultbmi <- adultbmi[adultbmi$sjlid %in% PHENO.ANY_SN$sjlid[i], c("sjlid", "DateVisitStart", "BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")]
  tmp.adultbmi <- tmp.adultbmi [tmp.adultbmi$DateVisitStart < PHENO.ANY_SN$gradedt[i],]
# If there are multiple dates, I will take the closest date to (on or before) gradedt
  print(paste0(i, "--", nrow(tmp.adultbmi)))
  tmp.adultbmi <- tmp.adultbmi [which.closest(tmp.adultbmi$DateVisitStart, PHENO.ANY_SN$gradedt[i]),]
  adultbmi.wanted <- rbind.data.frame(adultbmi.wanted,tmp.adultbmi)
}
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, adultbmi.wanted[match(PHENO.ANY_SN$sjlid, adultbmi.wanted$sjlid), c("DateVisitStart", "BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")])
PHENO.ANY_SN$obesity <- ifelse(PHENO.ANY_SN$BMI >= 30, "Y", "N")


#############
## Any SMNs
#############
# This will include any SNs excluding NMSCs
table(ANY_SNs$diaggrp)
SMNs <- ANY_SNs[!grepl("basal cell|melanoma|squamous cell", ANY_SNs$diag, ignore.case = T),]
# Excluding this sample: SJL5352907
SMNs <- SMNs[!grepl("SJL5352907", SMNs$sjlid),]
nrow(SMNs)
# 413
table(SMNs$diaggrp)

##########
## NMSCs
##########
# This will include basal cell, squamous cell and melanoma
table(ANY_SNs$diaggrp)
# tt <- ANY_SNs[grepl("basal|melanoma|squamous", ANY_SNs$diag, ignore.case = T),]
NMSCs <- ANY_SNs[grepl("basal cell|melanoma|squamous cell", ANY_SNs$diag, ignore.case = T),]
tt$sjlid[!tt$sjlid %in% NMSCs$sjlid]
nrow(NMSCs)
# 222
table(NMSCs$diaggrp)

##################
## Breast cancer
##################
BREASTcancer <- ANY_SNs[grepl("breast", ANY_SNs$diaggrp, ignore.case = T),]
nrow(BREASTcancer)
table(BREASTcancer$diaggrp)

##################
## Thyroid cancer
##################
THYROIDcancer <- ANY_SNs[grepl("thyroid", ANY_SNs$diaggrp, ignore.case = T),]
nrow(THYROIDcancer)
table(THYROIDcancer$diaggrp)


###############
## Meningioma
###############
MENINGIOMA <- ANY_SNs[grepl("meningioma", ANY_SNs$diaggrp, ignore.case = T),]
nrow(MENINGIOMA)
table(MENINGIOMA$diaggrp)


############
## Sarcoma
############
SARCOMA <- ANY_SNs[grepl("sarcoma", ANY_SNs$diaggrp, ignore.case = T),]
nrow(SARCOMA)
table(SARCOMA$diaggrp)





# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/phenotype_cleaning_attr_fraction.RDATA")


