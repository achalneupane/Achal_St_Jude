
library(haven)
library(benchmarkme)
library(dplyr)
library(plyr)
library(data.table)
library (birk)
library(gtools)
# benchmarkme::get_ram()


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


# #############################
# ## Adult habits/ Lifestyle ##
# #############################
# ## For each samples, get habits immediately after 18 years of age in agesurvey
# 
# # adlthabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.sas7bdat")
# adlthabits <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.txt", header = T, sep ="\t")
# head(adlthabits)
# # remove duplicated rows
# adlthabits <- distinct(adlthabits)
# adlthabits$agesurvey <- floor(adlthabits$agesurvey_diff_of_dob_datecomp)
# 
# samples.sjlife <- unique(adlthabits$SJLIFEID)
# 
# lifestyle <- {}
# for (i in 1:length(samples.sjlife)){
#   print(paste0("Doing ", i))
#   dat <- adlthabits[adlthabits$SJLIFEID == samples.sjlife[i],]
#   if (max(dat$agesurvey) >= 18){
#     print("YES")
#     dat2 <- dat[dat$agesurvey >= 18,]
#     lifestyle.tmp <- dat2[which(dat2$agesurvey == min(dat2$agesurvey)),]
#     # } else {
#     #   print("NO")
#     #   lifestyle.tmp <-  dat[which(dat$agesurvey == max(dat$agesurvey)),]
#     #   lifestyle.tmp[9:ncol(lifestyle.tmp)] <- NA
#   }
#   lifestyle <- rbind.data.frame(lifestyle, lifestyle.tmp)
# }
# 
# sum(duplicated(lifestyle$SJLIFEID))
# lifestyle$SJLIFEID[duplicated(lifestyle$SJLIFEID)]
# ## Remove duplicate row
# lifestyle <- lifestyle[!duplicated(lifestyle$SJLIFEID),]
# 
# ## Add all samples
# # lifestyle <- cbind.data.frame(wgspop[,1:2], lifestyle[match(wgspop$MRN, lifestyle$mrn), ])
# # lifestyle <- lifestyle[-c(3,4)]
# # tt <- lifestyle
# ## Recode categorical variables
# lifestyle$relation[lifestyle$relation == 1] <- "Self"
# lifestyle$relation[lifestyle$relation == 2] <- "Parent"
# lifestyle$relation[lifestyle$relation == 3] <- "Other"
# 
# ## Recode smoker
# lifestyle$smoker[lifestyle$smoker == 1] <- "Past"
# lifestyle$smoker[lifestyle$smoker == 2] <- "Current"
# lifestyle$smoker[lifestyle$smoker == 3] <- "Never"
# lifestyle$smoker <- ifelse(lifestyle$smoker == "Current", 1, 0)
# 
# ## Recode to Y/N
# lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 1 ] <- 1
# lifestyle[grepl("nopa|ltpa", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa", colnames(lifestyle))] == 2 ] <- 0
# lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 0 ] <- 0
# 
# # change the format of dates YYYY-MM-DD
# lifestyle$datecomp <- gsub("\\/", "-", lifestyle$datecomp)
# lifestyle$datecomp <- paste(sapply(strsplit(lifestyle$datecomp, "-"), `[`, 3), sapply(strsplit(lifestyle$datecomp, "-"), `[`, 1), sapply(strsplit(lifestyle$datecomp, "-"), `[`, 2), sep ="-")
# 
# lifestyle$dob <- gsub("\\/", "-", lifestyle$dob)
# lifestyle$dob <- paste(sapply(strsplit(lifestyle$dob, "-"), `[`, 3), sapply(strsplit(lifestyle$dob, "-"), `[`, 1), sapply(strsplit(lifestyle$dob, "-"), `[`, 2), sep ="-")
# 
# #######################
# ## Adolescent habits ##
# #######################
# 
# adolhabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adolhabits.sas7bdat")
# head(adolhabits)
# 
# ###############
# ## Adult BMI ##
# ###############
# 
# adultbmi <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.sas7bdat")
# head(adultbmi)
# 
# ## Add BMI and Nutrition to Lifestyle. Here, I am extracting the the BMI and Nutrition for the same age for each sample 
# lifestyle$BMI_KEY <- paste(lifestyle$SJLIFEID, lifestyle$agesurvey, sep = ":")
# 
# length(unique(adultbmi$sjlid))
# # 3640
# adultbmi$BMI_KEY <- paste(adultbmi$sjlid, adultbmi$AGE, sep = ":")
# sum(lifestyle$BMI_KEY %in% adultbmi$BMI_KEY)
# # 2964 
# ## samples that did not match by corresponding age
# cc <- lifestyle[!lifestyle$BMI_KEY %in% adultbmi$BMI_KEY,]
# 
# lifestyle <- cbind.data.frame(lifestyle, adultbmi[match(lifestyle$BMI_KEY, adultbmi$BMI_KEY), c("BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")])
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

########################
## Merge Genetic data ##
########################
QIN_vars <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/Na_Qin_vars/Qin_et_al_all_vars_final_recodeA.raw", header = T)
QIN_vars <- QIN_vars[-grep("FID|PAT|MAT|SEX|PHENOTYPE", colnames(QIN_vars))]
QIN_vars$Qin_Non.Ref.Counts <- rowSums(QIN_vars[-1]) 
clinical.dat$Qin_Non.Ref.Counts <- QIN_vars$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, QIN_vars$IID)]
clinical.dat$Qin_carriers <- ifelse(clinical.dat$Qin_Non.Ref.Counts > 0, "Y", "N")

Zhaoming_vars <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/Zhaoming_Wang_vars/Zhaoming_et_al_all_vars_final_recodeA.raw", header = T)
Zhaoming_vars <- Zhaoming_vars[-grep("FID|PAT|MAT|SEX|PHENOTYPE", colnames(Zhaoming_vars))]
Zhaoming_vars$Zhaoming_Non.Ref.Counts <- rowSums(Zhaoming_vars[-1]) 
clinical.dat$Zhaoming_Non.Ref.Counts <- Zhaoming_vars$Zhaoming_Non.Ref.Counts [match(clinical.dat$sjlid, Zhaoming_vars$IID)]
clinical.dat$Zhaoming_carriers <- ifelse(clinical.dat$Zhaoming_Non.Ref.Counts > 0, "Y", "N")

# Add ethnicity from PCA
EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_EUR_Per_PCA.txt", stringsAsFactors = F, header = F)
AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_AFR_Per_PCA.txt", stringsAsFactors = F, header = F)
EAS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_EAS_Per_PCA.txt", stringsAsFactors = F, header = F)
sum((clinical.dat$sjlid %in% EUR$V1))

clinical.dat$PCA.ethnicity <- ifelse(clinical.dat$sjlid %in% EUR$V1, "EUR", "")
clinical.dat$PCA.ethnicity[!grepl("EUR", clinical.dat$PCA.ethnicity)]    <- ifelse(clinical.dat$sjlid[!grepl("EUR", clinical.dat$PCA.ethnicity)] %in% AFR$V1, "AFR", "")
clinical.dat$PCA.ethnicity[!grepl("EUR|AFR", clinical.dat$PCA.ethnicity)]    <- ifelse(clinical.dat$sjlid[!grepl("EUR|AFR", clinical.dat$PCA.ethnicity)] %in% EAS$V1, "EAS", "")

# test <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/test.txt", header =T, check.names = F,   colClasses = c("character"))
# test <- as.data.frame(t(test))
# test$row.NAMES <- as.character(rownames(test))
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

# Found two duplicate samples with exact same values, removing them
PHENO.ANY_SN <- PHENO.ANY_SN[!duplicated(PHENO.ANY_SN$sjlid),]


# ## Merge lifestyle
# PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, lifestyle[match(PHENO.ANY_SN$sjlid, lifestyle$SJLIFEID),c("datecomp", "agesurvey", "relation", "smoker", "nopa", "ltpa", "drk5", "bingedrink", "heavydrink", "riskydrink")])
# 
# # If ANY_SN is Yes and the survey date of lifestyle (datecomp) is later than the SN grade
# # date (gradedt), then the lifestyle variables for those samples will be
# # irrelevant
# PHENO.ANY_SN[which(PHENO.ANY_SN[PHENO.ANY_SN$ANY_SN == "Y, "datecomp"] >  PHENO.ANY_SN[PHENO.ANY_SN$ANY_SN == "Y", "gradedt"]), c("smoker", "nopa", "ltpa", "drk5", "bingedrink", "heavydrink", "riskydrink")] <- NA
# 
# # Now adding BMI and nutrition; extracting the BMI values immediately before or on the date of SN gradedt
# adultbmi.wanted <- {}
# for (i in 1:length(PHENO.ANY_SN$sjlid)){
#   tmp.adultbmi <- adultbmi[adultbmi$sjlid %in% PHENO.ANY_SN$sjlid[i], c("sjlid", "DateVisitStart", "BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")]
#   tmp.adultbmi <- tmp.adultbmi [tmp.adultbmi$DateVisitStart < PHENO.ANY_SN$gradedt[i],]
# # If there are multiple dates, I will take the closest date to (on or before) gradedt
#   print(paste0(i, "--", nrow(tmp.adultbmi)))
#   tmp.adultbmi <- tmp.adultbmi [which.closest(tmp.adultbmi$DateVisitStart, PHENO.ANY_SN$gradedt[i]),]
#   adultbmi.wanted <- rbind.data.frame(adultbmi.wanted,tmp.adultbmi)
# }
# PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, adultbmi.wanted[match(PHENO.ANY_SN$sjlid, adultbmi.wanted$sjlid), c("DateVisitStart", "BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")])
# PHENO.ANY_SN$obesity <- ifelse(PHENO.ANY_SN$BMI >= 30, "Y", "N")


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




#################################################################################################
####################################### ANALYSIS ################################################
#################################################################################################
################################
## Replicating Zhaoming et al ##
################################
colnames(PHENO.ANY_SN)

# PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "gender", "agedx", "AnyRT", "anyrt_prim", "anyrt_5", "anyrt_10", "brainorheadrt_yn",
#                                "brainrt_yn", "maxseg1dose", "maxseg2dose", "maxseg3dose", "maxseg4dose", "maxsegrtdose", "chestrt_yn",
#                                "maxchestrtdose", "neckrt_yn", "maxneckrtdose", "pelvisrt_yn", "maxpelvisrtdose", "abdomenrt_yn", "maxabdrtdose",
#                                "aa_class_dose_any", "aa_hvymtl_dose_any", "carbo_dose_any", "cisplat_dose_any", "cisplateq_dose_any",
#                                "epitxn_dose_any", "Qin_carriers", "Zhaoming_carriers",
#                                "PCA.ethnicity", "ANY_SN", "AGE.ANY_SN" )]


PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "gender", "agedx", "brainrt_yn", "chestrt_yn", "neckrt_yn", "pelvisrt_yn", "abdomenrt_yn", 
                               "aa_class_dose_any", "epitxn_dose_any", "Zhaoming_carriers", "PCA.ethnicity", "ANY_SN", "AGE.ANY_SN")]

PHENO.ANY_SN$ANY_SN <- factor(PHENO.ANY_SN$ANY_SN)

## Gene mutation
PHENO.ANY_SN$Zhaoming_carriers <- factor(PHENO.ANY_SN$Zhaoming_carriers, levels = c("N", "Y"))

## Age at diagnosis
PHENO.ANY_SN$agedx <- floor(PHENO.ANY_SN$agedx)
PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- PHENO.ANY_SN$agedx
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 0 & PHENO.ANY_SN$agedx <= 4 ] <- "0-4"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 5 & PHENO.ANY_SN$agedx <= 9 ] <- "5-9"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 10 & PHENO.ANY_SN$agedx <= 14 ] <- "10-14"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 15 ] <- "≥15"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- factor(PHENO.ANY_SN$AGE_AT_DIAGNOSIS, levels = c("0-4", "5-9", "10-14", "≥15")) # first level will be treated as reference
# PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- relevel(PHENO.ANY_SN$AGE_AT_DIAGNOSIS, ref = "0-4")

## Sex
PHENO.ANY_SN$gender <- factor(PHENO.ANY_SN$gender, levels = c("Male", "Female"))

## Radiation
PHENO.ANY_SN$brainrt_yn <- factor(PHENO.ANY_SN$brainrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$neckrt_yn <- factor(PHENO.ANY_SN$neckrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$chestrt_yn <- factor(PHENO.ANY_SN$chestrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$abdomenrt_yn <- factor(PHENO.ANY_SN$abdomenrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$pelvisrt_yn <- factor(PHENO.ANY_SN$pelvisrt_yn, levels = c("N", "Y"))

## Alkylating agents
PHENO.ANY_SN$Alkylating_agent [PHENO.ANY_SN$aa_class_dose_any == 0] <- "None"
PHENO.ANY_SN$Alkylating_agent [PHENO.ANY_SN$aa_class_dose_any > 0 & PHENO.ANY_SN$aa_class_dose_any <= 6617] <- "1st"
PHENO.ANY_SN$Alkylating_agent [PHENO.ANY_SN$aa_class_dose_any > 6617 & PHENO.ANY_SN$aa_class_dose_any <= 10500] <- "2nd"
PHENO.ANY_SN$Alkylating_agent [PHENO.ANY_SN$aa_class_dose_any > 10500] <- "3rd"
PHENO.ANY_SN$Alkylating_agent <- factor(PHENO.ANY_SN$Alkylating_agent, levels = c("None", "1st", "2nd", "3rd"))

## Anthracyclines


## Epidophyllotoxin
PHENO.ANY_SN$Epidophyllotoxin [PHENO.ANY_SN$epitxn_dose_any == 0] <- "None"
PHENO.ANY_SN$Epidophyllotoxin [PHENO.ANY_SN$epitxn_dose_any > 0 & PHENO.ANY_SN$epitxn_dose_any <= 1588] <- "1st"
PHENO.ANY_SN$Epidophyllotoxin [PHENO.ANY_SN$epitxn_dose_any > 1588 & PHENO.ANY_SN$epitxn_dose_any <= 6704] <- "2nd"
PHENO.ANY_SN$Epidophyllotoxin [PHENO.ANY_SN$epitxn_dose_any > 6704] <- "3rd"
PHENO.ANY_SN$Epidophyllotoxin <- factor(PHENO.ANY_SN$Epidophyllotoxin, levels = c("None", "1st", "2nd", "3rd"))



PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]



mod1 <- glm(ANY_SN ~ Qin_carriers + gender + agedx + AnyRT + maxseg1dose, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)

mod1 <- glm(ANY_SN ~ . , family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1)

mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_DIAGNOSIS + gender + brainrt_yn + neckrt_yn + chestrt_yn + abdomenrt_yn + pelvisrt_yn + Alkylating_agent + Epidophyllotoxin, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1)




# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/phenotype_cleaning_attr_fraction.RDATA")


