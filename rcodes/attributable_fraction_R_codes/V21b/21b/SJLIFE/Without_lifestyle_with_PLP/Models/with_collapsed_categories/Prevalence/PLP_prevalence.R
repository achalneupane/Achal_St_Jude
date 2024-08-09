
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

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
wanted.samples <- PHENO.ANY_SN$sjlid
length(wanted.samples)

#####################
## wgspop.sas7bdat ##
#####################
# WORKDIR:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE")


########################
## Merge Genetic data ##
########################

QIN_vars <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/Na_Qin_vars/Qin_et_al_all_vars_final_recodeA.raw", header = T)
QIN_vars <- QIN_vars[QIN_vars$IID %in% wanted.samples,]
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


#######################
## Zhaoming variants ## 
#######################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")

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
# 166 variants listed in table A4
tableA10 <- zhaoming_tables[grepl ("A10", zhaoming_tables$Table),]
sum(unique(search_list$KEY.varID) %in% unique(tableA4$FULLKEY))
# 148 unique variants in table A4
search_list <- search_list[search_list$KEY.varID %in% tableA4$FULLKEY,]


# sum(tableA4$Gene %in% tableA10$Gene)
# sum(tableA4$FULLKEY %in% tableA10$FULLKEY)


# Keeping variants in Table A4 only
Zhaoming_vars <- Zhaoming_vars[c(1,which(colnames(Zhaoming_vars) %in% search_list$ID))]
Zhaoming_vars$Zhaoming_Non.Ref.Counts <- rowSums(Zhaoming_vars[-1]) 


prevalence.counts <- sum(Zhaoming_vars$Zhaoming_Non.Ref.Counts > 0)
table(ifelse(Zhaoming_vars$Zhaoming_Non.Ref.Counts > 0, "Y", "N"))
prop.test(prevalence.counts, 4401)
# 4.2%

###########################
## 2. Qin carriers (all) ##
###########################
#0n 05/17/2024
prevalence.counts <- sum(PHENO.ANY_SN$Qin_Non.Ref.Counts > 0)
prop.test(prevalence.counts, 4401)
# 11.11%


QIN_vars.wo.zhaoming <- QIN_vars[!colnames(QIN_vars) %in% colnames(Zhaoming_vars)]
QIN_vars.wo.zhaoming$Qin_without_Zhaoming_vars.Non.Ref.Counts <- rowSums(QIN_vars.wo.zhaoming[-1]) 

prevalence.counts <- sum(QIN_vars.wo.zhaoming$Qin_without_Zhaoming_vars.Non.Ref.Counts > 0)
table(ifelse(QIN_vars.wo.zhaoming$Qin_without_Zhaoming_vars.Non.Ref.Counts > 0, "Y", "N"))
prop.test(prevalence.counts, 4401)
# 10.18%



### OlD analysis
## SJLIFE (ALL) 
mod1 <- glm(ANY_SN ~ Qin_carriers + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
mod1.EUR <- glm(ANY_SN ~ Qin_carriers + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)

## Prevalence
prevalence.counts <- sum(PHENO.ANY_SN$Qin_Non.Ref.Counts > 0)
prop.test(prevalence.counts, 4401)
# SJLIFE (11.1% prevalence; 95% CI, 10% to 11.8%)
# QIN: (11.5% prevalence; 95% CI, 10.6% to 12.5%)

## Prevalence HR pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.HR.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4401)
# SJLIFE: (3.9% prevalence; 95% CI, 3.4% to 4.5%)
# QIN: (4.2% prevalence; 95% CI, 3.6% to 4.8%)

## Prevalence MMR pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.MMR.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4401)
# SJLIFE: (0.7% prevalence; 95% CI, 0.4% to 1.0%)
# QIN: (0.8% prevalence; 95% CI, 0.6% to 1.1%)

## Prevalence NER pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.NER.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4401)
# SJLIFE: (2.04% prevalence; 95% CI, 1.65% to 2.5%)
# QIN:  (2.2% prevalence; 95% CI, 1.8% to 2.7%)

## Prevalence FA pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.FA.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4401)
# SJLIFE: (2.9% prevalence; 95% CI, 2.48% to 3.49%)
# QIN:  ( 3.2% prevalence; 95% CI, 2.7% to 3.7%)

## Prevalence NHEJ pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.NHEJ.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4401)
# SJLIFE: (0.5% prevalence; 95% CI, 0.3% to 0.7%)
# QIN:  ( 0.6% prevalence; 95% CI, 0.4% to 0.9%)

## Prevalence BER pathways
prevalence.counts <- sum(PHENO.ANY_SN$Qin_carriers.BER.pathways == "Y", na.rm = T)
prop.test(prevalence.counts, 4401)
# SJLIFE: (2.4% prevalence; 95% CI, 2.0% to 2.9%)
# QIN:  ( 2.5% prevalence; 95% CI, 2.1% to 3.0%)


####################
## 3. HR Pathways ##
####################
## Checking with Qi's data
# sum(ANY_SNs$sjlid %in% qi.df.SN.filtered$sjlid)
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
ANY_SNs <- ANY_SNs[ANY_SNs$sjlid %in% qi.df.SN.filtered$sjlid,]
dim(ANY_SNs)
# 491
PHENO.ANY_SN$ANY_SN <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))

## SJLIFE (ALL) 
mod1 <- glm(ANY_SN ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)
estimates <- as.data.frame(summary(mod1)$coefficients)[c(1,2)]
estimates$Estimate <- as.numeric(estimates$Estimate)
estimates$OR <- round(exp(estimates$Estimate),2)
## CI
# exp(5.319e-01-(1.96*1.009e-01))
# exp(5.319e-01+(1.96*1.009e-01))
estimates$S.error <- as.numeric(as.character(estimates$`Std. Error`))
CI <-  paste0("(", paste0(round(exp(estimates$Estimate - (1.96*estimates$S.error)),1), " to ", round(exp(estimates$Estimate + (1.96*estimates$S.error)),1)),")")
estimates$OR <- paste(estimates$OR, CI, sep = " ")
estimates

# write.table(estimates, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Additional_files/estimates_Qin.txt", sep = "\t", col.names = T, row.names = T, quote = F)


## SJLIFE (EUR)
mod1.EUR <- glm(ANY_SN ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


#####################
## 4. MMR Pathways ##
#####################
mod1 <- glm(ANY_SN ~ Qin_carriers.MMR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

# tt <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$Qin_carriers.MMR.pathways),]
## SJL5450006 seems to be missing in genetic data

# tt <- tt[c("ANY_SN", "Qin_carriers.MMR.pathways", "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", "AGE_AT_DIAGNOSIS", "gender", "maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_any.category")]
# 
# mod1 <- glm(ANY_SN ~ Qin_carriers.MMR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = tt)
# summary(mod1)

# sum(is.na(tt$Qin_carriers.MMR.pathways))
# sum(is.na(tt$AGE_AT_LAST_CONTACT.cs1))
# sum(is.na(tt$AGE_AT_LAST_CONTACT.cs2))
# sum(is.na(tt$AGE_AT_LAST_CONTACT.cs3))
# sum(is.na(tt$AGE_AT_LAST_CONTACT.cs4))
# sum(is.na(tt$AGE_AT_DIAGNOSIS))
# sum(is.na(tt$gender))
# sum(is.na(tt$maxsegrtdose.category))
# sum(is.na(tt$maxabdrtdose.category))
# sum(is.na(tt$maxchestrtdose.category))
# sum(is.na(tt$epitxn_dose_any.category))



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
## Contingency table/cross tab of categorical variables (Zhaoming)
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





##################################
## Sanity check with Qin's data ##
##################################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/phenotype_cleaning_attr_fraction.RDATA")
pheno.crosscheck <- PHENO.ANY_SN
pheno.crosscheck <- pheno.crosscheck[mixedorder(pheno.crosscheck$sjlid),]
qi.df <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/toyadav.sas7bdat")
qi.df <- qi.df[mixedorder(qi.df$sjlid),]
head(qi.df)
qi.df.unique <- qi.df[!duplicated(qi.df$sjlid),]
# length(unique(qi.df.unique$sjlid))
# # 4402
# sum(pheno.crosscheck$sjlid %in% qi.df.unique$sjlid)
# # 4402
# qi.df <- qi.df[!duplicated(qi.df$sjlid),]
table(qi.df$MutationStatus)
table(qi.df$sn)
# 1263

qi.df.SN <- qi.df[grepl(1, qi.df$sn),]
dim(qi.df.SN)
# 1263
# Keeping unique rows
qi.df.SN <- distinct(qi.df.SN)
dim(qi.df.SN)



qi.df.not.SN <- qi.df[!grepl(1, qi.df$sn),]
dim(qi.df.not.SN)
# 3911
# 4402-3911 = 491

# Extracting SNs for the first time from Qi's data
qi.df.SN.filtered <- setDT(qi.df.SN)[,.SD[which.min(evaldt)],by=sjlid][order(evaldt, decreasing = FALSE)]
dim(qi.df.SN.filtered)

# Now merge back the data with SN and without SN 
qi.df.final <- rbind.data.frame(qi.df.not.SN, qi.df.SN.filtered)
dim(qi.df.final)
qi.df.final <- qi.df.final[mixedorder(qi.df.final$sjlid),]

paste0(table(qi.df.final$anyrt_5)[2], ", ", table(pheno.crosscheck$anyrt_5)[2]) # anyrt_5
# "2169, 2174"
paste0(table(qi.df.final$brainrt_yn)[2], ", ", table(pheno.crosscheck$brainrt_yn)[2]) # Brain
# "858, 1095"
paste0(table(qi.df.final$neckrt_yn)[2], ", ", table(pheno.crosscheck$neckrt_yn)[2]) # Neck
# "601, 797"
paste0(table(qi.df.final$chest_yn)[2], ", ", table(pheno.crosscheck$chestrt_yn)[2]) # Chest
# 628, 861
paste0(table(qi.df.final$abdomenrt_yn)[2], ", ", table(pheno.crosscheck$abdomenrt_yn)[2]) # Abdomen
# "567, 790"
paste0(table(qi.df.final$pelvis_yn)[2], ", ", table(pheno.crosscheck$pelvisrt_yn)[2]) # Pelvis
# "506, 688"
paste0(table(qi.df.final$aaclassic_5 ==1)[2], ", ", table(pheno.crosscheck$aa_class_dose_5_yn)[2]) # Alkylating_5
# "2480, 2473"


paste0(table(qi.df.final$anthracyclines_5 > 0)[2], ", ", table(pheno.crosscheck$anthra_jco_dose_5_yn)[2]) # Anthracyclines_5
# "2455, 2451"
paste0(table(qi.df.final$epipodophyllotoxins_5 > 0)[2], ", ", table(pheno.crosscheck$epitxn_dose_5_yn)[2]) # Epipodophyllotoxins_5
# "1496, 1488"

# Age at diagnosis
paste0(paste0(paste0( median(qi.df.final$agedx),  " (", paste0( quantile(qi.df.final$agedx, 1 / 4), "???",
        quantile(qi.df.final$agedx, 3 / 4)), ")" )), ", ", paste0(paste0(
        round(median(pheno.crosscheck$agedx), 2), " (", paste0(round(quantile( pheno.crosscheck$agedx, 1 / 4 ), 2), "???",
        round(quantile( pheno.crosscheck$agedx, 3 / 4), 2), ")")))) 
# "6.275 (2.79???12.525), 6.28 (2.79???12.52)"


table(qi.df.final$diaggrp == pheno.crosscheck$diaggrp)
# FALSE  TRUE 
# 1087  3315

table(qi.df.final$epitxn_dose)
sum(is.na(qi.df.final$epipodophyllotoxins_5))
# 10 
qi.df.final[is.na(qi.df.final$epipodophyllotoxins_5),c("sjlid", "epipodophyllotoxins_5")]

pheno.crosscheck.SN.491 <- pheno.crosscheck[pheno.crosscheck$sjlid %in%  qi.df.SN.filtered$sjlid,]
pheno.crosscheck.SN.491 <- pheno.crosscheck.SN.491[mixedorder(pheno.crosscheck.SN.491$sjlid),]
table(pheno.crosscheck.SN.491$diaggrp == qi.df.SN.filtered$diaggrp)
# FALSE  TRUE 
# 376   115 

pheno.crosscheck.not.SN <- pheno.crosscheck[pheno.crosscheck$sjlid %in%  qi.df.not.SN$sjlid,]
pheno.crosscheck.not.SN <- pheno.crosscheck.not.SN[mixedorder(pheno.crosscheck.not.SN$sjlid),]
table(pheno.crosscheck.not.SN$diaggrp == qi.df.not.SN$diaggrp)

# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/phenotype_cleaning_attr_fraction.RDATA")


