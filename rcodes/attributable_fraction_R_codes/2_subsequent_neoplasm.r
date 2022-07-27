#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/demographic.RDATA")
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
#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


# Removing samples with with SN within the 5 years of childhood cancer
ANY_SNs <- ANY_SNs[!ANY_SNs$sjlid %in% subneo.within5$sjlid,]
dim(ANY_SNs)
# 612
##############
## Any SMNs ##
##############
# This will include any SNs excluding NMSCs
SMNs <- subneo[!grepl("basal cell|squamous cell", subneo$diag, ignore.case = T),]
SMNs <- setDT(SMNs)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(SMNs)
# 482
table(SMNs$diaggrp)

# Removing samples with SNs within 5 years of childhood cancer
SMNs <- SMNs[!SMNs$sjlid %in% subneo.within5$sjlid,]
nrow(SMNs)
# 459
###########
## NMSCs ##
###########
# This will include basal cell, squamous cell and melanoma
NMSCs <- subneo[grepl("basal cell|squamous cell", subneo$diag, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(NMSCs)
# 255
table(NMSCs$diaggrp)

# Removing samples with SNs within 5 years of childhood cancer
NMSCs <- NMSCs[!NMSCs$sjlid %in% subneo.within5$sjlid,]
nrow(NMSCs)
# 252

###################
## Breast cancer ##
###################
BREASTcancer <- subneo[grepl("breast", subneo$diag, ignore.case = T),]
BREASTcancer <- setDT(BREASTcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(THYROIDcancer)
# 86
# Removing samples with SNs within 5 years of childhood cancer
BREASTcancer <- BREASTcancer[!BREASTcancer$sjlid %in% subneo.within5$sjlid,]
nrow(BREASTcancer)
# 77

##################
## Thyroid cancer
##################
THYROIDcancer <- subneo[grepl("thyroid", subneo$diag, ignore.case = T),]
THYROIDcancer <- setDT(THYROIDcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(THYROIDcancer)
# Removing samples with SNs within 5 years of childhood cancer
THYROIDcancer <- THYROIDcancer[!THYROIDcancer$sjlid %in% subneo.within5$sjlid,]
nrow(THYROIDcancer)


###############
## Meningioma
###############
MENINGIOMA <- subneo[grepl("meningioma", subneo$diag, ignore.case = T),]
MENINGIOMA <- setDT(MENINGIOMA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(MENINGIOMA)
# 150
# Removing samples with SNs within 5 years of childhood cancer
MENINGIOMA <- MENINGIOMA[!MENINGIOMA$sjlid %in% subneo.within5$sjlid,]
nrow(MENINGIOMA)

#############
## Sarcoma ##
#############
SARCOMA <- subneo[grepl("sarcoma", subneo$diag, ignore.case = T),]
SARCOMA <- setDT(SARCOMA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]
nrow(SARCOMA)
# 36
# Removing samples with SNs within 5 years of childhood cancer
SARCOMA <- SARCOMA[!SARCOMA$sjlid %in% subneo.within5$sjlid,]
nrow(SARCOMA)
# 33

#####################
## Add SN variable ##
#####################
# 1. Any SN

PHENO.ANY_SN$ANY_SN <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ANY_SNs[match(PHENO.ANY_SN$sjlid, ANY_SNs$sjli), c("gradedt", "AGE.ANY_SN")])



########################################
## MODEL TEST for Zhaoming's variants ##
########################################
## SJLIFE (ALL)
# mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + aa_class_dose_any_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN)
mod1 <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
# mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + Epidophyllotoxin + aa_class_dose_any_yn + cisplat_dose_any_yn + aa_hvymtl_dose_any_yn, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
mod1.EUR <- glm(ANY_SN ~ Zhaoming_carriers + AGE_AT_LAST_CONTACT + AGE_AT_DIAGNOSIS + gender + brainrt_yn + chestrt_yn + abdomenrt_yn + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


## Prevalence
prevalence.counts <- sum(Zhaoming_vars$Zhaoming_Non.Ref.Counts > 0)
prop.test(prevalence.counts, 4507)


###################################
## MODEL TEST for Qin's variants ##
###################################

###########################
## 1. Qin baseline model ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(ANY_SN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)
estimates <- as.data.frame(summary(mod1)$coefficients)[c(1,2)]
estimates$Estimate <- as.numeric(estimates$Estimate)
estimates$RR <- round(exp(estimates$Estimate),2)
## CI
# exp(5.319e-01-(1.96*1.009e-01))
# exp(5.319e-01+(1.96*1.009e-01))
estimates$S.error <- as.numeric(as.character(estimates$`Std. Error`))
CI <-  paste0("(", paste0(round(exp(estimates$Estimate - (1.96*estimates$S.error)),1), " to ", round(exp(estimates$Estimate + (1.96*estimates$S.error)),1)),")")

estimates$RR <- paste(estimates$RR, CI, sep = " ")
write.table(estimates, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Additional_files/estimates_Qin.txt", sep = "\t", col.names = T, row.names = T, quote = F)

## SJLIFE (EUR)
mod1.EUR <- glm(ANY_SN ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


###########################
## 2. Qin carriers (all) ##
###########################
## SJLIFE (ALL) 
mod1 <- glm(ANY_SN ~ Qin_carriers + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
summary(mod1)

## SJLIFE (EUR)
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
## 3. HR Pathways ##
####################
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


## SJLIFE (EUR)
mod1.EUR <- glm(ANY_SN ~ Qin_carriers.HR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN.EUR)
summary(mod1.EUR)


#####################
## 4. MMR Pathways ##
#####################
mod1 <- glm(ANY_SN ~ Qin_carriers.MMR.pathways + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_any.category, family = binomial(link = "logit"), data = PHENO.ANY_SN)
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





