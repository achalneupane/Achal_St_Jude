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


##########
## Drug ##
##########

# drug <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/drug.sas7bdat")
drug <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_drug.txt", sep = "\t", stringsAsFactors = F)
head(drug)
# drug$cisplat_dose_any

## Add drug to clinical data
clinical.dat <- cbind.data.frame(clinical.dat, drug[match(clinical.dat$MRN, drug$MRN), grep("dose_any|_5",colnames(drug))])

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
QIN.Pathways <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Additional_files/Qin_variant_pathways.txt", header = T)
QIN.Pathways <- QIN.Pathways[QIN.Pathways$MATCH_YN == "Y",]

# HR pathway
HR.pathways <- QIN.Pathways[grepl("HR", QIN.Pathways$DNA_Repair_Pathway),]
HR.pathways <- HR.pathways[!duplicated(HR.pathways$VarKEY_IN_SJLIFE),]
HR.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% HR.pathways$VarKEY_IN_SJLIFE))]
# variants in SJLIFE
ncol(HR.pathways)-1
# 154
# variants in QIN
# 157

HR.pathways$Qin_Non.Ref.Counts <- rowSums(HR.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.HR.pathways <- HR.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, HR.pathways$IID)]
clinical.dat$Qin_carriers.HR.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.HR.pathways == 0, "N", "Y"))

# FA Pathway
FA.pathways <- QIN.Pathways[grepl("FA", QIN.Pathways$DNA_Repair_Pathway),]
FA.pathways <- FA.pathways[!duplicated(FA.pathways$VarKEY_IN_SJLIFE),]
FA.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% FA.pathways$VarKEY_IN_SJLIFE))]
# variants in SJLIFE
ncol(FA.pathways)-1
# 101
# variants in QIN
# 104

FA.pathways$Qin_Non.Ref.Counts <- rowSums(FA.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.FA.pathways <- FA.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, FA.pathways$IID)]
clinical.dat$Qin_carriers.FA.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.FA.pathways == 0, "N", "Y"))

# MMR Pathway
MMR.pathways <- QIN.Pathways[grepl("MMR", QIN.Pathways$DNA_Repair_Pathway),]
MMR.pathways <- MMR.pathways[!duplicated(MMR.pathways$VarKEY_IN_SJLIFE),]
MMR.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% MMR.pathways$VarKEY_IN_SJLIFE))]
# variants in SJLIFE
ncol(MMR.pathways)-1
# 28
# variants in QIN
# 30

MMR.pathways$Qin_Non.Ref.Counts <- rowSums(MMR.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.MMR.pathways <- MMR.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, MMR.pathways$IID)]
clinical.dat$Qin_carriers.MMR.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.MMR.pathways == 0, "N", "Y"))

# BER Pathway
BER.pathways <- QIN.Pathways[grepl("BER", QIN.Pathways$DNA_Repair_Pathway),]
BER.pathways <- BER.pathways[!duplicated(BER.pathways$VarKEY_IN_SJLIFE),]
BER.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% BER.pathways$VarKEY_IN_SJLIFE))]
# variants in SJLIFE
ncol(BER.pathways)-1
# 66
# variants in QIN
# 67


BER.pathways$Qin_Non.Ref.Counts <- rowSums(BER.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.BER.pathways <- BER.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, BER.pathways$IID)]
clinical.dat$Qin_carriers.BER.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.BER.pathways == 0, "N", "Y"))

# NER Pathway
NER.pathways <- QIN.Pathways[grepl("NER", QIN.Pathways$DNA_Repair_Pathway),]
NER.pathways <- NER.pathways[!duplicated(NER.pathways$VarKEY_IN_SJLIFE),]
NER.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% NER.pathways$VarKEY_IN_SJLIFE))]
# variants in SJLIFE
ncol(NER.pathways)-1
# 73
# variants in QIN
# 76

NER.pathways$Qin_Non.Ref.Counts <- rowSums(NER.pathways[-1]) 
clinical.dat$Qin_Non.Ref.Counts.NER.pathways <- NER.pathways$Qin_Non.Ref.Counts [match(clinical.dat$sjlid, NER.pathways$IID)]
clinical.dat$Qin_carriers.NER.pathways <- factor(ifelse(clinical.dat$Qin_Non.Ref.Counts.NER.pathways == 0, "N", "Y"))

# NHEJ Pathway
NHEJ.pathways <- QIN.Pathways[grepl("NHEJ", QIN.Pathways$DNA_Repair_Pathway),]
NHEJ.pathways <- NHEJ.pathways[!duplicated(NHEJ.pathways$VarKEY_IN_SJLIFE),]
NHEJ.pathways <- QIN_vars[c(1,which(colnames(QIN_vars) %in% NHEJ.pathways$VarKEY_IN_SJLIFE))]
# variants in SJLIFE
ncol(NHEJ.pathways)-1
# 18
# variants in QIN
# 18

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
Sys.setlocale("LC_ALL", "C")
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


PHENO.ANY_SN <- clinical.dat
# Found two duplicate samples with exact same values, removing them
PHENO.ANY_SN <- PHENO.ANY_SN[!duplicated(PHENO.ANY_SN$sjlid),]

## Source variable for subsequent neoplasm from R code: subsequent_neoplasm.r


#################################################################################################
####################################### ANALYSIS ################################################
#################################################################################################
## Here, I am taking survivors only for the analysis
wgs.pop.sjlife <- wgspop[grepl("SJLIFE", wgspop$wgs_cohort),]
PHENO.ANY_SN <- PHENO.ANY_SN[PHENO.ANY_SN$sjlid %in% wgs.pop.sjlife$sjlid,]

colnames(PHENO.ANY_SN)
saved.PHENO.ANY_SN <- PHENO.ANY_SN


###################################################

PHENO.ANY_SN <- PHENO.ANY_SN[c("sjlid", "MRN", "gender", "agedx", "diaggrp", "agelstcontact", "AnyRT", "anyrt_5", "brainrt_yn", "maxsegrtdose", "chestrt_yn", "maxchestrtdose", "anthra_jco_dose_any", "anthra_jco_dose_5","neckrt_yn", 
                               "maxneckrtdose", "pelvisrt_yn", "maxpelvisrtdose","abdomenrt_yn", "maxabdrtdose", "aa_class_dose_any", "aa_class_dose_5", "epitxn_dose_any", "epitxn_dose_5",
                               "cisplat_dose_any", "cisplateq_dose_5", "aa_hvymtl_dose_any", "aa_hvymtl_dose_5", "Zhaoming_carriers", "Qin_carriers", "Qin_carriers.HR.pathways", "Qin_carriers.FA.pathways",
                               "Qin_carriers.MMR.pathways", "Qin_carriers.BER.pathways", "Qin_carriers.NER.pathways", "Qin_carriers.NHEJ.pathways", "PCA.ethnicity")]


## Age at diagnosis
# PHENO.ANY_SN$agedx <- floor(PHENO.ANY_SN$agedx)
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 0 & PHENO.ANY_SN$agedx < 5 ] <- "0-4"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 5 & PHENO.ANY_SN$agedx < 10 ] <- "5-9"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 10 & PHENO.ANY_SN$agedx < 15 ] <- "10-14"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS[PHENO.ANY_SN$agedx >= 15 ] <- ">=15"
PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- factor(PHENO.ANY_SN$AGE_AT_DIAGNOSIS, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference


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
PHENO.ANY_SN$aa_class_dose_any_yn <-  factor(ifelse(PHENO.ANY_SN$aa_class_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$aa_class_dose_any[PHENO.ANY_SN$aa_class_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$aa_class_dose_any.category <- cut(PHENO.ANY_SN$aa_class_dose_any, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)
levels(PHENO.ANY_SN$aa_class_dose_any.category) <- c(levels(PHENO.ANY_SN$aa_class_dose_any.category), "Unknown")
PHENO.ANY_SN$aa_class_dose_any.category [is.na(PHENO.ANY_SN$aa_class_dose_any.category)] <- "Unknown"


PHENO.ANY_SN$aa_class_dose_5_yn <-  factor(ifelse(PHENO.ANY_SN$aa_class_dose_5 == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$aa_class_dose_5[PHENO.ANY_SN$aa_class_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$aa_class_dose_5.category <- cut(PHENO.ANY_SN$aa_class_dose_5, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)
levels(PHENO.ANY_SN$aa_class_dose_5.category) <- c(levels(PHENO.ANY_SN$aa_class_dose_5.category), "Unknown")
PHENO.ANY_SN$aa_class_dose_5.category [is.na(PHENO.ANY_SN$aa_class_dose_5.category)] <- "Unknown"


## Platinum agents (Y/N)
PHENO.ANY_SN$cisplat_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$cisplat_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$cisplat_dose_any[PHENO.ANY_SN$cisplat_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$cisplat_dose_any.category <- cut(PHENO.ANY_SN$cisplat_dose_any, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)
levels(PHENO.ANY_SN$cisplat_dose_any.category) <- c(levels(PHENO.ANY_SN$cisplat_dose_any.category), "Unknown")
PHENO.ANY_SN$cisplat_dose_any.category [is.na(PHENO.ANY_SN$cisplat_dose_any.category)] <- "Unknown"


PHENO.ANY_SN$cisplateq_dose_5_yn <- factor(ifelse(PHENO.ANY_SN$cisplateq_dose_5 == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$cisplateq_dose_5[PHENO.ANY_SN$cisplateq_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$cisplateq_dose_5.category <- cut(PHENO.ANY_SN$cisplateq_dose_5, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)
levels(PHENO.ANY_SN$cisplateq_dose_5.category) <- c(levels(PHENO.ANY_SN$cisplateq_dose_5.category), "Unknown")
PHENO.ANY_SN$cisplateq_dose_5.category [is.na(PHENO.ANY_SN$cisplateq_dose_5.category)] <- "Unknown"


## Heavy metals (Y/N and Tertiles)
PHENO.ANY_SN$aa_hvymtl_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$aa_hvymtl_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$aa_hvymtl_dose_any[PHENO.ANY_SN$aa_hvymtl_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$aa_hvymtl_dose_any.category <- cut(PHENO.ANY_SN$aa_hvymtl_dose_any, breaks = c(0, 0.001, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)
levels(PHENO.ANY_SN$aa_hvymtl_dose_any.category) <- c(levels(PHENO.ANY_SN$aa_hvymtl_dose_any.category), "Unknown")
PHENO.ANY_SN$aa_hvymtl_dose_any.category [is.na(PHENO.ANY_SN$aa_hvymtl_dose_any.category)] <- "Unknown"


PHENO.ANY_SN$aa_hvymtl_dose_5_yn <- factor(ifelse(PHENO.ANY_SN$aa_hvymtl_dose_5 == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$aa_hvymtl_dose_5[PHENO.ANY_SN$aa_hvymtl_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$aa_hvymtl_dose_5.category <- cut(PHENO.ANY_SN$aa_hvymtl_dose_5, breaks = c(0, 0.001, TERT),
                                                labels = c("None", "1st", "2nd", "3rd"),
                                                include.lowest = TRUE)
levels(PHENO.ANY_SN$aa_hvymtl_dose_5.category) <- c(levels(PHENO.ANY_SN$aa_hvymtl_dose_5.category), "Unknown")
PHENO.ANY_SN$aa_hvymtl_dose_5.category [is.na(PHENO.ANY_SN$aa_hvymtl_dose_5.category)] <- "Unknown"


## Anthracyclines (Y/N and Tertiles)
PHENO.ANY_SN$anthra_jco_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$anthra_jco_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$anthra_jco_dose_any[PHENO.ANY_SN$anthra_jco_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$anthra_jco_dose_any.category <- cut(PHENO.ANY_SN$anthra_jco_dose_any, breaks = c(0, 0.001, TERT),
                                                labels = c("None", "1st", "2nd", "3rd"),
                                                include.lowest = TRUE)
levels(PHENO.ANY_SN$anthra_jco_dose_any.category) <- c(levels(PHENO.ANY_SN$anthra_jco_dose_any.category), "Unknown")
PHENO.ANY_SN$anthra_jco_dose_any.category [is.na(PHENO.ANY_SN$anthra_jco_dose_any.category)] <- "Unknown"


PHENO.ANY_SN$anthra_jco_dose_5_yn <- factor(ifelse(PHENO.ANY_SN$anthra_jco_dose_5 == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$anthra_jco_dose_5[PHENO.ANY_SN$anthra_jco_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$anthra_jco_dose_5.category <- cut(PHENO.ANY_SN$anthra_jco_dose_5, breaks = c(0, 0.001, TERT),
                                                 labels = c("None", "1st", "2nd", "3rd"),
                                                 include.lowest = TRUE)
levels(PHENO.ANY_SN$anthra_jco_dose_5.category) <- c(levels(PHENO.ANY_SN$anthra_jco_dose_5.category), "Unknown")
PHENO.ANY_SN$anthra_jco_dose_5.category [is.na(PHENO.ANY_SN$anthra_jco_dose_5.category)] <- "Unknown"


## Epidophyllotoxin
PHENO.ANY_SN$epitxn_dose_any_yn <- factor(ifelse(PHENO.ANY_SN$epitxn_dose_any == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$epitxn_dose_any[PHENO.ANY_SN$epitxn_dose_any !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$epitxn_dose_any.category <- cut(PHENO.ANY_SN$epitxn_dose_any, breaks = c(0, 0.001, TERT),
                                                labels = c("None", "1st", "2nd", "3rd"),
                                                include.lowest = TRUE)
levels(PHENO.ANY_SN$epitxn_dose_any.category) <- c(levels(PHENO.ANY_SN$epitxn_dose_any.category), "Unknown")
PHENO.ANY_SN$epitxn_dose_any.category [is.na(PHENO.ANY_SN$epitxn_dose_any.category)] <- "Unknown"


PHENO.ANY_SN$epitxn_dose_5_yn <- factor(ifelse(PHENO.ANY_SN$epitxn_dose_5 == 0, "N", "Y"))

TERT = unname(quantile(PHENO.ANY_SN$epitxn_dose_5[PHENO.ANY_SN$epitxn_dose_5 !=0], c(1/3, 2/3, 1), na.rm = T))
PHENO.ANY_SN$epitxn_dose_5.category <- cut(PHENO.ANY_SN$epitxn_dose_5, breaks = c(0, 0.001, TERT),
                                             labels = c("None", "1st", "2nd", "3rd"),
                                             include.lowest = TRUE)
levels(PHENO.ANY_SN$epitxn_dose_5.category) <- c(levels(PHENO.ANY_SN$epitxn_dose_5.category), "Unknown")
PHENO.ANY_SN$epitxn_dose_5.category [is.na(PHENO.ANY_SN$epitxn_dose_5.category)] <- "Unknown"

save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/1_demographics.RDATA")


