
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

# drug <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/drug.sas7bdat")
drug <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_drug.txt", sep = "\t", stringsAsFactors = F)
head(drug)
# drug$cisplat_dose_any

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
colnames(QIN_vars) <- str_split(gsub("\\.", ":", colnames(QIN_vars)), "_", simplify=T)[,1]
## Qin's Pathways
QIN.Pathways <- read.delim("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Na_Qin_Pathogenic Germline Mutations in DNA Repair Genes in Combination With Cancer Treatment Exposures/Qin_variant_pathways.txt", header = T)
QIN.Pathways <- QIN.Pathways[QIN.Pathways$MATCH_YN == "Y",]

# HR pathway
HR.pathways <- QIN.Pathways[grepl("HR", QIN.Pathways$DNA_Repair_Pathway),]
HR.pathways <- HR.pathways[!duplicated(HR.pathways$VarKEY_IN_SJLIFE),]
HR.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% HR.pathways$VarKEY_IN_SJLIFE))]

HR.pathways$Qin_Non.Ref.Counts <- rowSums(HR.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.HR.pathways <- HR.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, HR.pathways$IID)]
clinical.dat$Qin_carriers.HR.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.HR.pathways == 0, "N", "Y"))

# FA Pathway
FA.pathways <- QIN.Pathways[grepl("FA", QIN.Pathways$DNA_Repair_Pathway),]
FA.pathways <- FA.pathways[!duplicated(FA.pathways$VarKEY_IN_SJLIFE),]
FA.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% FA.pathways$VarKEY_IN_SJLIFE))]

FA.pathways$Qin_Non.Ref.Counts <- rowSums(FA.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.FA.pathways <- FA.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, FA.pathways$IID)]
clinical.dat$Qin_carriers.FA.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.FA.pathways == 0, "N", "Y"))

# MMR Pathway
MMR.pathways <- QIN.Pathways[grepl("MMR", QIN.Pathways$DNA_Repair_Pathway),]
MMR.pathways <- MMR.pathways[!duplicated(MMR.pathways$VarKEY_IN_SJLIFE),]
MMR.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% MMR.pathways$VarKEY_IN_SJLIFE))]

MMR.pathways$Qin_Non.Ref.Counts <- rowSums(MMR.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.MMR.pathways <- MMR.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, MMR.pathways$IID)]
clinical.dat$Qin_carriers.MMR.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.MMR.pathways == 0, "N", "Y"))

# BER Pathway
BER.pathways <- QIN.Pathways[grepl("BER", QIN.Pathways$DNA_Repair_Pathway),]
BER.pathways <- BER.pathways[!duplicated(BER.pathways$VarKEY_IN_SJLIFE),]
BER.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% BER.pathways$VarKEY_IN_SJLIFE))]

BER.pathways$Qin_Non.Ref.Counts <- rowSums(BER.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.BER.pathways <- BER.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, BER.pathways$IID)]
clinical.dat$Qin_carriers.BER.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.BER.pathways == 0, "N", "Y"))

# NER Pathway
NER.pathways <- QIN.Pathways[grepl("NER", QIN.Pathways$DNA_Repair_Pathway),]
NER.pathways <- NER.pathways[!duplicated(NER.pathways$VarKEY_IN_SJLIFE),]
NER.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% NER.pathways$VarKEY_IN_SJLIFE))]

NER.pathways$Qin_Non.Ref.Counts <- rowSums(NER.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.NER.pathways <- NER.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, NER.pathways$IID)]
clinical.dat$Qin_carriers.NER.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.NER.pathways == 0, "N", "Y"))

# NHEJ Pathway
NHEJ.pathways <- QIN.Pathways[grepl("NHEJ", QIN.Pathways$DNA_Repair_Pathway),]
NHEJ.pathways <- NHEJ.pathways[!duplicated(NHEJ.pathways$VarKEY_IN_SJLIFE),]
NHEJ.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% NHEJ.pathways$VarKEY_IN_SJLIFE))]

NHEJ.pathways$Qin_Non.Ref.Counts <- rowSums(NHEJ.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.NHEJ.pathways <- NHEJ.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, NHEJ.pathways$IID)]
clinical.dat$Qin_carriers.NHEJ.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.NHEJ.pathways == 0, "N", "Y"))

QIN_vars$Qin_Non.Ref.Counts <- rowSums(QIN_vars[-1]) 
clinical.dat$Qin_Non.Ref.Counts <- QIN_vars$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, QIN_vars$IID)]
clinical.dat$Qin_carriers <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts == 0, "N", "Y"))





Zhaoming_vars <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/Zhaoming_Wang_vars/Zhaoming_et_al_all_vars_final_recodeA.raw", header = T)
Zhaoming_vars <- Zhaoming_vars[-grep("FID|PAT|MAT|SEX|PHENOTYPE", colnames(Zhaoming_vars))]
# edit column names
colnames(Zhaoming_vars) <- str_split(gsub("\\.", ":", colnames(Zhaoming_vars)), "_", simplify=T)[,1]

## 
# Sys.setlocale("LC_ALL", "C")
zhaoming_tables <- read.delim("/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Zhaoming_Wang_Genetic_risk_for_subsequent_neoplasms_JCO.2018.77.8589/P-PL-sjlife-genetics-sn_SNV_INDELS_07_13_2022.txt", sep = "\t")
dim(zhaoming_tables)
sum(!duplicated(zhaoming_tables$FULLKEY))
search_list <- read.delim("/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Zhaoming_Wang_Genetic_risk_for_subsequent_neoplasms_JCO.2018.77.8589/search_list.txt", sep = "\t")
sum(unique(zhaoming_tables$FULLKEY) %in% search_list$KEY.varID)
# 299
tableA4 <- zhaoming_tables[grepl ("A4", zhaoming_tables$Table),]
dim(tableA4)
# 166
tableA10 <- zhaoming_tables[grepl ("A10", zhaoming_tables$Table),]
sum(unique(search_list$KEY.varID) %in% unique(tableA4$FULLKEY))
# 148
search_list <- search_list[search_list$KEY.varID %in% tableA4$FULLKEY,]


# sum(tableA4$Gene %in% tableA10$Gene)
# sum(tableA4$FULLKEY %in% tableA10$FULLKEY)


# Keeping variants in Table A4 only
Zhaoming_vars <- Zhaoming_vars[c(1,which(colnames(Zhaoming_vars) %in% search_list$ID))]
Zhaoming_vars$Zhaoming_Non.Ref.Counts <- rowSums(Zhaoming_vars[-1]) 

clinical.dat$Zhaoming_Non.Ref.Counts <- Zhaoming_vars$Zhaoming_Non.Ref.Counts [match(clinical.dat$sjlid, Zhaoming_vars$IID)]
clinical.dat$Zhaoming_carriers <- factor(ifelse(clinical.dat$Zhaoming_Non.Ref.Counts == 0, "N", "Y"))


# Add ethnicity from PCA
EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_EUR_Per_PCA.txt", stringsAsFactors = F, header = F)
AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_AFR_Per_PCA.txt", stringsAsFactors = F, header = F)
EAS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_EAS_Per_PCA.txt", stringsAsFactors = F, header = F)
sum((clinical.dat$sjlid %in% EUR$V1))
# 3358

clinical.dat$PCA.ethnicity <- ifelse(clinical.dat$sjlid %in% EUR$V1, "EUR", "")
clinical.dat$PCA.ethnicity[!grepl("EUR", clinical.dat$PCA.ethnicity)]    <- ifelse(clinical.dat$sjlid[!grepl("EUR", clinical.dat$PCA.ethnicity)] %in% AFR$V1, "AFR", "")
clinical.dat$PCA.ethnicity[!grepl("EUR|AFR", clinical.dat$PCA.ethnicity)]    <- ifelse(clinical.dat$sjlid[!grepl("EUR|AFR", clinical.dat$PCA.ethnicity)] %in% EAS$V1, "EAS", "")


dim(QIN_vars)
dim(clinical.dat)

sum(QIN_vars$IID %in% clinical.dat$sjlid)
# 4401
# test <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/test.txt", header =T, check.names = F,   colClasses = c("character"))
# test <- as.data.frame(t(test))
# test$row.NAMES <- as.character(rownames(test))
#########################
## Subsequent Neoplasm ##
#########################

subneo <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/subneo.sas7bdat")
head(subneo)
table(subneo$diaggrp)

# add diagnosis date 
subneo$diagdt <-  clinical.dat$diagdt [match(subneo$sjlid , clinical.dat$sjlid)]
subneo$agedx <-  clinical.dat$agedx [match(subneo$sjlid , clinical.dat$sjlid)]
# add DOB
subneo$DOB <- demog$dob[match(subneo$MRN, demog$MRN)]

subneo$AGE.ANY_SN <- time_length(interval(as.Date(subneo$DOB), as.Date(subneo$gradedt)), "years")
# subneo$AGE.exact.ANY_SN <- time_length(interval(as.Date(subneo$DOB), as.Date(subneo$gradedt)), "years")
# subneo$AGE.ANY_SN <- floor(subneo$AGE.exact.ANY_SN)

subneo$AGE.ANY_SN.after.childhood.cancer <- time_length(interval(as.Date(subneo$diagdt), as.Date(subneo$gradedt)), "years")
# subneo$AGE.ANY_SN.after.childhood.cancer <- floor(subneo$AGE.exact.ANY_SN.after.childhood.cancer)

subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


########################################
# How many SNs after 5 years
sjlife1.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife_1/sample_list.txt")
sjlife2.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife_2/sample_list.txt")

subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$sjlid))
# 619
subneo.after5.sjlife1 <- subneo.after5[subneo.after5$sjlid %in% sjlife1.samples$V1,]
length(unique(subneo.after5.sjlife1$sjlid))
# 562

subneo.after5.sjlife2 <- subneo.after5[subneo.after5$sjlid %in% sjlife2.samples$V1,]
length(unique(subneo.after5.sjlife2$sjlid))
# 45

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
dim(subneo.within5)
# 24
############
## Any SNs 
############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


# Removing samples with with SN within the 5 years of childhood cancer
ANY_SNs <- ANY_SNs[!ANY_SNs$sjlid %in% subneo.within5$sjlid,]
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
saved.PHENO.ANY_SN <- PHENO.ANY_SN

# PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "gender", "agedx", "agelstcontact", "AnyRT", "anyrt_prim", "anyrt_5", "anyrt_10", "brainorheadrt_yn",
#                                "brainrt_yn", "maxseg1dose", "maxseg2dose", "maxseg3dose", "maxseg4dose", "maxsegrtdose", "chestrt_yn",
#                                "maxchestrtdose", "neckrt_yn", "maxneckrtdose", "pelvisrt_yn", "maxpelvisrtdose", "abdomenrt_yn", "maxabdrtdose",
#                                "aa_class_dose_any", "aa_hvymtl_dose_any", "carbo_dose_any", "cisplat_dose_any", "cisplateq_dose_any",
#                                "epitxn_dose_any", "Qin_carriers", "Zhaoming_carriers",
#                                "PCA.ethnicity", "ANY_SN", "AGE.ANY_SN" )]


# PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "gender", "agedx", "agelstcontact", "brainrt_yn", "chestrt_yn", "neckrt_yn", "pelvisrt_yn", "abdomenrt_yn", 
#                                "aa_class_dose_any", "epitxn_dose_any", "cisplat_dose_any", "aa_hvymtl_dose_any", "Zhaoming_carriers", "Qin_carriers", "PCA.ethnicity", "ANY_SN", "AGE.ANY_SN")]

PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "gender", "agedx", "agelstcontact", "brainrt_yn", "maxsegrtdose", "chestrt_yn", "maxchestrtdose", "anthra_jco_dose_any","neckrt_yn", 
                               "maxneckrtdose", "pelvisrt_yn", "maxpelvisrtdose","abdomenrt_yn", "maxabdrtdose", "aa_class_dose_any", "epitxn_dose_any",
                               "cisplat_dose_any", "aa_hvymtl_dose_any", "Zhaoming_carriers", "Qin_carriers", "Qin_carriers.HR.pathways", "Qin_carriers.FA.pathways",
                               "Qin_carriers.MMR.pathways", "Qin_carriers.BER.pathways", "Qin_carriers.NER.pathways", "Qin_carriers.NHEJ.pathways", "PCA.ethnicity", "ANY_SN", "AGE.ANY_SN")]


PHENO.ANY_SN$ANY_SN <- factor(PHENO.ANY_SN$ANY_SN, levels = c("N", "Y"))

## Age at diagnosis
# PHENO.ANY_SN$agedx <- floor(PHENO.ANY_SN$agedx)
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 0 & PHENO.ANY_SN$agedx < 5 ] <- "0-4"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 5 & PHENO.ANY_SN$agedx < 10 ] <- "5-9"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 10 & PHENO.ANY_SN$agedx < 15 ] <- "10-14"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 15 ] <- ">=15"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- factor(PHENO.ANY_SN$AGE_AT_DIAGNOSIS, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference
# PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- relevel(PHENO.ANY_SN$AGE_AT_DIAGNOSIS, ref = "0-4")

# PHENO.ANY_SN$AGE_AT_DIAGNOSIS.1 <- cut(PHENO.ANY_SN$agedx, breaks = c(0, 4.9999, 9.9999, 14.9999, max(PHENO.ANY_SN$agedx, na.rm = T)),
#                                           labels = c("0-4", "5-9", "10-14", ">=15"),
#                                           include.lowest = TRUE)



## Age at last contact
PHENO.ANY_SN$AGE_AT_LAST_CONTACT[PHENO.ANY_SN$agelstcontact >= 0 & PHENO.ANY_SN$agelstcontact < 25 ] <- "0-24"
PHENO.ANY_SN$AGE_AT_LAST_CONTACT[PHENO.ANY_SN$agelstcontact >= 25 & PHENO.ANY_SN$agelstcontact < 35 ] <- "25-34"
PHENO.ANY_SN$AGE_AT_LAST_CONTACT[PHENO.ANY_SN$agelstcontact >= 35 & PHENO.ANY_SN$agelstcontact < 45 ] <- "35-44"
PHENO.ANY_SN$AGE_AT_LAST_CONTACT[PHENO.ANY_SN$agelstcontact >= 45 ] <- ">=45"
PHENO.ANY_SN$AGE_AT_LAST_CONTACT <- factor(PHENO.ANY_SN$AGE_AT_LAST_CONTACT, levels = c("0-24", "25-34", "35-44", ">=45")) # first level will be treated as reference


## Age at last contact (cubic spline)
source("https://raw.githubusercontent.com/achalneupane/Achal_St_Jude/main/rcodes/cubic_spline.r")

breaks = seq(5, 95, 22.5)

cp = quantile(PHENO.ANY_SN$agelstcontact, breaks/100, na.rm = T)

cs = cubic_spline(PHENO.ANY_SN$agelstcontact, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, cs)
# Merge cs to your original data.frame and adjust in the logistic regression




## Sex
PHENO.ANY_SN$gender <- factor(PHENO.ANY_SN$gender, levels = c("Male", "Female"))

## Radiation
PHENO.ANY_SN$brainrt_yn <- factor(PHENO.ANY_SN$brainrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$neckrt_yn <- factor(PHENO.ANY_SN$neckrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$chestrt_yn <- factor(PHENO.ANY_SN$chestrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$abdomenrt_yn <- factor(PHENO.ANY_SN$abdomenrt_yn, levels = c("N", "Y"))
PHENO.ANY_SN$pelvisrt_yn <- factor(PHENO.ANY_SN$pelvisrt_yn, levels = c("N", "Y"))

####################
## Radiation dose ##
####################
## Maxsegrtdose (Cranial RT)
PHENO.ANY_SN$maxsegrtdose.category <- cut(PHENO.ANY_SN$maxsegrtdose, breaks = c(0, 200, 1799, 2999, max(PHENO.ANY_SN$maxsegrtdose, na.rm = T)),
    labels = c("None", ">0-<18", ">=18-<30", ">=30"),
    include.lowest = TRUE)
levels(PHENO.ANY_SN$maxsegrtdose.category) <- c(levels(PHENO.ANY_SN$maxsegrtdose.category), "Unknown")
PHENO.ANY_SN$maxsegrtdose.category [is.na(PHENO.ANY_SN$maxsegrtdose)] <- "Unknown"
# ## Sanity check
# # tt <- PHENO.ANY_SN[c("maxsegrtdose", "maxsegrtdose.category")]
# # sum(is.na(PHENO.ANY_SN$maxsegrtdose.category ))
# # cc <- tt[tt$maxsegrtdose < 200,]
# # cc <- tt[tt$maxsegrtdose >= 2999 ,]
# # cc <- tt[is.na(tt$maxsegrtdose),]
# # is.na(PHENO.ANY_SN$maxsegrtdose.category)

## Abdomen RT
PHENO.ANY_SN$maxabdrtdose.category <- cut(PHENO.ANY_SN$maxabdrtdose, breaks = c(0, 200, 2999, max(PHENO.ANY_SN$maxabdrtdose, na.rm = T)),
                                          labels = c("None", ">0-<30", ">=30"),
                                          include.lowest = TRUE)
levels(PHENO.ANY_SN$maxabdrtdose.category) <- c(levels(PHENO.ANY_SN$maxabdrtdose.category), "Unknown")
PHENO.ANY_SN$maxabdrtdose.category [is.na(PHENO.ANY_SN$maxabdrtdose)] <- "Unknown"

## Pelvis RT
PHENO.ANY_SN$maxpelvisrtdose.category <- cut(PHENO.ANY_SN$maxpelvisrtdose, breaks = c(0, 200, 1999, max(PHENO.ANY_SN$maxpelvisrtdose, na.rm = T)),
                                          labels = c("None", ">0-<20", ">=20"),
                                          include.lowest = TRUE)
levels(PHENO.ANY_SN$maxpelvisrtdose.category) <- c(levels(PHENO.ANY_SN$maxpelvisrtdose.category), "Unknown")
PHENO.ANY_SN$maxpelvisrtdose.category [is.na(PHENO.ANY_SN$maxpelvisrtdose)] <- "Unknown"

## Chest RT
PHENO.ANY_SN$maxchestrtdose.category <- cut(PHENO.ANY_SN$maxchestrtdose, breaks = c(0, 200, 1999, max(PHENO.ANY_SN$maxchestrtdose, na.rm = T)),
                                             labels = c("None", ">0-<20", ">=20"),
                                             include.lowest = TRUE)
levels(PHENO.ANY_SN$maxchestrtdose.category) <- c(levels(PHENO.ANY_SN$maxchestrtdose.category), "Unknown")
PHENO.ANY_SN$maxchestrtdose.category [is.na(PHENO.ANY_SN$maxchestrtdose)] <- "Unknown"


## Neck RT
PHENO.ANY_SN$maxneckrtdose.category <- cut(PHENO.ANY_SN$maxneckrtdose, breaks = c(0, 200, 1099, 1999, 2999, max(PHENO.ANY_SN$maxneckrtdose, na.rm = T)),
                                          labels = c("None", ">0-<11", ">=11-<20", ">=20-<30", ">=30"),
                                          include.lowest = TRUE)
levels(PHENO.ANY_SN$maxneckrtdose.category) <- c(levels(PHENO.ANY_SN$maxneckrtdose.category), "Unknown")
PHENO.ANY_SN$maxneckrtdose.category [is.na(PHENO.ANY_SN$maxneckrtdose)] <- "Unknown"




## Alkylating agents (Y/N and Tertiles)
PHENO.ANY_SN$Alkylating_agent_yn <-  factor(ifelse(PHENO.ANY_SN$aa_class_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$aa_class_dose_any[PHENO.ANY_SN$aa_class_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$Alkylating_agent.category <- cut(PHENO.ANY_SN$aa_class_dose_any, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)


## Platinum agents (Y/N and Tertiles)
PHENO.ANY_SN$cisplat_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$cisplat_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$cisplat_dose_any[PHENO.ANY_SN$cisplat_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$cisplat_dose_any.category <- cut(PHENO.ANY_SN$cisplat_dose_any, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)


## Heavy metals (Y/N and Tertiles)
PHENO.ANY_SN$aa_hvymtl_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$aa_hvymtl_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$aa_hvymtl_dose_any[PHENO.ANY_SN$aa_hvymtl_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$aa_hvymtl_dose_any.category <- cut(PHENO.ANY_SN$aa_hvymtl_dose_any, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)


## Anthracyclines (Y/N and Tertiles)
PHENO.ANY_SN$anthra_jco_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$anthra_jco_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$anthra_jco_dose_any[PHENO.ANY_SN$anthra_jco_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$anthra_jco_dose_any.category <- cut(PHENO.ANY_SN$anthra_jco_dose_any, breaks = c(0, 0.001, TERT),
                                                labels = c("None", "1st", "2nd", "3rd"),
                                                include.lowest = TRUE)

## Epidophyllotoxin
PHENO.ANY_SN$epitxn_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$epitxn_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$epitxn_dose_any[PHENO.ANY_SN$epitxn_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$epitxn_dose_any.category <- cut(PHENO.ANY_SN$epitxn_dose_any, breaks = c(0, 0.001, TERT),
                                                labels = c("None", "1st", "2nd", "3rd"),
                                                include.lowest = TRUE)


PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]

########################################
## MODEL TEST for Zhaoming's variants ##
########################################
## SJLIFE (ALL)
# mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + Alkylating_agent_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN)
mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
# mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + Alkylating_agent_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


## Prevalence
prevalence.counts <- sum(Zhaoming_vars$Zhaoming_Non.Ref.Counts > 0)
prop.test(prevalence.counts, 4507)


###################################
## MODEL TEST for Qin's variants ##
###################################

###########################
## 1. Qin carriers (all) ##
###########################
## SJLIFE (ALL) 
# mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + Alkylating_agent_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN)
mod1 <- glm(ANY_SN ~ Qin_carriers + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
# mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + Alkylating_agent_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
mod1.EUR <- glm(ANY_SN ~ Qin_carriers + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)

## Prevalence
prevalence.counts <- sum(QIN_vars$Qin_Non.Ref.Counts > 0)
prop.test(prevalence.counts, 4507)
# SJLIFE (10.89% prevalence; 95% CI, 10% to 11.8%)
# QIN: (11.5% prevalence; 95% CI, 10.6% to 12.5%)

## Prevalence HR pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.HR.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4507)
# SJLIFE: (3.9% prevalence; 95% CI, 3.4% to 4.5%)
# QIN: (4.2% prevalence; 95% CI, 3.6% to 4.8%)

## Prevalence MMR pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.MMR.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4507)
# SJLIFE: (0.7% prevalence; 95% CI, 0.4% to 1.0%)
# QIN: (0.8% prevalence; 95% CI, 0.6% to 1.1%)

## Prevalence NER pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.NER.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4507)
# SJLIFE: (2.04% prevalence; 95% CI, 1.65% to 2.5%)
# QIN:  (2.2% prevalence; 95% CI, 1.8% to 2.7%)

## Prevalence FA pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.FA.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4507)
# SJLIFE: (2.9% prevalence; 95% CI, 2.48% to 3.49%)
# QIN:  ( 3.2% prevalence; 95% CI, 2.7% to 3.7%)

## Prevalence NHEJ pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.NHEJ.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4507)
# SJLIFE: (0.5% prevalence; 95% CI, 0.3% to 0.7%)
# QIN:  ( 0.6% prevalence; 95% CI, 0.4% to 0.9%)

## Prevalence BER pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.BER.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4507)
# SJLIFE: (2.4% prevalence; 95% CI, 2.0% to 2.9%)
# QIN:  ( 2.5% prevalence; 95% CI, 2.1% to 3.0%)


####################
## 2. HR Pathways ##
####################
## SJLIFE (ALL) 
# mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + Alkylating_agent_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN)
mod1 <- glm(ANY_SN ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
# mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + Alkylating_agent_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
mod1.EUR <- glm(ANY_SN ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)







#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
## Sjlife 1 sample list (used by Zhaoming); checking only in these samples
PHENO.ANY_SN.sjlife1 <- PHENO.ANY_SN[PHENO.ANY_SN$sjlid %in% sjlife1.samples$V1,]
PHENO.ANY_SN.sjlife1.EUR <- PHENO.ANY_SN.sjlife1[PHENO.ANY_SN.sjlife1$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN.sjlife1))]

mod.sjlife1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin, family = binomial(link = "logit"), data = PHENO.ANY_SN.sjlife1)
summary(mod.sjlife1)

mod.sjlife1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin, family = binomial(link = "logit"), data = PHENO.ANY_SN.sjlife1.EUR)
summary(mod.sjlife1.EUR)

#########################################
## Sjlife 2 sample list (used by Zhaoming); checking only in these samples
PHENO.ANY_SN.sjlife2 <- PHENO.ANY_SN[PHENO.ANY_SN$sjlid %in% sjlife1.samples$V1,]
PHENO.ANY_SN.sjlife2.EUR <- PHENO.ANY_SN.sjlife2[PHENO.ANY_SN.sjlife2$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN.sjlife2))]

mod.sjlife2 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin, family = binomial(link = "logit"), data = PHENO.ANY_SN.sjlife2)
summary(mod.sjlife2)

mod.sjlife2.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin, family = binomial(link = "logit"), data = PHENO.ANY_SN.sjlife2.EUR)
summary(mod.sjlife2.EUR)

#########################################

## Sjlife 1 sample list (used by Qin); checking only in these samples
PHENO.ANY_SN.sjlife1 <- PHENO.ANY_SN[PHENO.ANY_SN$sjlid %in% sjlife1.samples$V1,]
PHENO.ANY_SN.sjlife1.EUR <- PHENO.ANY_SN.sjlife1[PHENO.ANY_SN.sjlife1$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN.sjlife1))]

mod.sjlife1 <- glm(ANY_SN ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.sjlife1)
summary(mod.sjlife1)

mod.sjlife1.EUR <- glm(ANY_SN ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.sjlife1.EUR)
summary(mod.sjlife1.EUR)

########################################
## cross tab of categorical variables (Zhaoming)
library(expss)

CROSS_CASES.df <- PHENO.ANY_SN[c("ANY_SN", "Zhaoming_carriers" , "AGE_AT_LAST_CONTACT", "AGE_AT_DIAGNOSIS", "gender", "brainrt_yn", "chestrt_yn", "abdomenrt_yn", "Epidophyllotoxin")]
CROSS_CASES.df <- apply_labels(CROSS_CASES.df,
             ANY_SN = "ANY_SN", Zhaoming_carriers = "Zhaoming_carriers", AGE_AT_LAST_CONTACT = "AGE_AT_LAST_CONTACT",
             AGE_AT_DIAGNOSIS = "AGE_AT_DIAGNOSIS", gender = "gender", brainrt_yn  = "brainrt_yn", chestrt_yn = "chestrt_yn", abdomenrt_yn = "abdomenrt_yn", Epidophyllotoxin = "Epidophyllotoxin")

CROSS_CASES.df %>%
cross_cases(ANY_SN, list(Zhaoming_carriers , AGE_AT_LAST_CONTACT, AGE_AT_DIAGNOSIS, gender, brainrt_yn, chestrt_yn, abdomenrt_yn, Epidophyllotoxin))
# cross_cases(PHENO.ANY_SN, ANY_SN, list(Zhaoming_carriers , AGE_AT_LAST_CONTACT, AGE_AT_DIAGNOSIS, gender, brainrt_yn, chestrt_yn, abdomenrt_yn, Epidophyllotoxin))


## cross tab of categorical variables (Qin)

CROSS_CASES.df <- PHENO.ANY_SN[c("ANY_SN", "Qin_carriers", "AGE_AT_LAST_CONTACT", "AGE_AT_DIAGNOSIS", "gender", "maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_any.category")]
CROSS_CASES.df <- apply_labels(CROSS_CASES.df,
                               ANY_SN = "ANY_SN", Qin_carriers = "Qin_carriers", AGE_AT_LAST_CONTACT = "AGE_AT_LAST_CONTACT",
                               AGE_AT_DIAGNOSIS = "AGE_AT_DIAGNOSIS", gender = "gender", maxsegrtdose.category = "maxsegrtdose.category",
                               maxabdrtdose.category = "maxabdrtdose.category", maxchestrtdose.category = "maxchestrtdose.category",
                               epitxn_dose_any.category = "epitxn_dose_any.category")

CROSS_CASES.df %>%
  cross_cases(ANY_SN, list(Qin_carriers, AGE_AT_LAST_CONTACT, AGE_AT_DIAGNOSIS, gender, maxsegrtdose.category, maxabdrtdose.category, maxchestrtdose.category, epitxn_dose_any.category))
# cross_cases(PHENO.ANY_SN, ANY_SN, list(Qin_carriers , AGE_AT_LAST_CONTACT, AGE_AT_DIAGNOSIS, gender, brainrt_yn, chestrt_yn, abdomenrt_yn, Epidophyllotoxin))


########################################
## Checking prevalence
# also checking prevalence
Zhaoming_vars.sjlife1 <- Zhaoming_vars[Zhaoming_vars$IID %in% sjlife1.samples$V1,]
table(ifelse(Zhaoming_vars.sjlife1$Zhaoming_Non.Ref.Counts > 0, "Y", "N"))
# N    Y 
# 2664  322 

prop.test(375,4132,correct=FALSE)



# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/phenotype_cleaning_attr_fraction.RDATA")


