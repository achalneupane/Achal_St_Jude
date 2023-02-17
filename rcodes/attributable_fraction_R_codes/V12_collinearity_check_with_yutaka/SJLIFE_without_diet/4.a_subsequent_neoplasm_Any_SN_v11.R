#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)
#########################
## Subsequent Neoplasm ##
#########################

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


subneo <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/subneo.sas7bdat")
head(subneo)
table(subneo$diaggrp)
dim(subneo)
# 1731 9

subneo <- subneo[subneo$sjlid %in% PHENO.ANY_SN$sjlid ,]
dim(subneo)
# 1717
# add diagnosis date 
subneo$diagdt <-  PHENO.ANY_SN$diagdt [match(subneo$sjlid , PHENO.ANY_SN$sjlid)]
subneo$agedx <-  PHENO.ANY_SN$agedx [match(subneo$sjlid , PHENO.ANY_SN$sjlid)]
# add DOB
subneo$DOB <- PHENO.ANY_SN$dob[match(subneo$sjlid, PHENO.ANY_SN$sjlid)]

subneo$AGE.ANY_SN <- time_length(interval(as.Date(subneo$DOB), as.Date(subneo$gradedt)), "years")

## These two dates should be the (almost) same
subneo$AGE.ANY_SN.after.childhood.cancer <- time_length(interval(as.Date(subneo$diagdt), as.Date(subneo$gradedt)), "years")
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$sjlid))
# 612

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 22
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
# 605

PHENO.ANY_SN$ANY_SN <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))

# write.table(PHENO.ANY_SN$sjlid, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/4401_attributable_fraction_ids.txt", col.names = F, row.names = F, quote = F)

#############################
## Add Lifestyle variables ##
#############################
# Define CA/CO status in lifestyle
ALL.LIFESTYLE$CACO <- factor(ifelse(!ALL.LIFESTYLE$SJLIFEID %in% ANY_SNs$sjlid, 0, 1))



## Get date (gradedt) and age at diagnosis of SN
ALL.LIFESTYLE$ANY.SN_gradedate <- ANY_SNs$gradedt[match(ALL.LIFESTYLE$SJLIFEID, ANY_SNs$sjlid)]
ALL.LIFESTYLE$AGE.ANY_SN <- ANY_SNs$AGE.ANY_SN[match(ALL.LIFESTYLE$SJLIFEID, ANY_SNs$sjlid)]

## In CASES, if age survey is greater than age at diagnosis; NULLIFY the favorable_lifestyle.category. That information is not useful
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("Current_smoker_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$PhysicalActivity_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("PhysicalActivity_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("RiskyHeavyDrink_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEALTHY_Diet_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEALTHY_Diet_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$Not_obese_yn_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("Obese_yn")] <- NA
ALL.LIFESTYLE[which(ALL.LIFESTYLE$CACO == 1 & ALL.LIFESTYLE$HEI2015_TOTAL_SCORE_agesurvey > ALL.LIFESTYLE$AGE.ANY_SN), c("HEI2015_TOTAL_SCORE")] <- NA


#############################
## Addd lifestyle to Pheno ##
#############################
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ALL.LIFESTYLE[match(PHENO.ANY_SN$sjlid, ALL.LIFESTYLE$SJLIFEID),c("HEI2015_TOTAL_SCORE", "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Obese_yn")])

# Count missing
PHENO.ANY_SN$missing.lifestyles <- rowSums(is.na(PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", "Obese_yn")]))
table(PHENO.ANY_SN$missing.lifestyles)

## Relevel 6 lifestyle variables
PHENO.ANY_SN$Current_smoker_yn[is.na(PHENO.ANY_SN$Current_smoker_yn)] <- "Unknown"
PHENO.ANY_SN$Current_smoker_yn <- factor(PHENO.ANY_SN$Current_smoker_yn, level = c("No", "Yes", "Unknown")) 

PHENO.ANY_SN$PhysicalActivity_yn[is.na(PHENO.ANY_SN$PhysicalActivity_yn)] <- "Unknown"
PHENO.ANY_SN$PhysicalActivity_yn <- factor(PHENO.ANY_SN$PhysicalActivity_yn, level = c("Yes", "No", "Unknown")) 

PHENO.ANY_SN$RiskyHeavyDrink_yn[is.na(PHENO.ANY_SN$RiskyHeavyDrink_yn)] <- "Unknown"
PHENO.ANY_SN$RiskyHeavyDrink_yn <- factor(PHENO.ANY_SN$RiskyHeavyDrink_yn, level = c("No", "Yes", "Unknown")) 

PHENO.ANY_SN$HEALTHY_Diet_yn[is.na(PHENO.ANY_SN$HEALTHY_Diet_yn)] <- "Unknown"
PHENO.ANY_SN$HEALTHY_Diet_yn <- factor(PHENO.ANY_SN$HEALTHY_Diet_yn, level = c("Yes", "No", "Unknown")) 

PHENO.ANY_SN$Obese_yn[is.na(PHENO.ANY_SN$Obese_yn)] <- "Unknown";
PHENO.ANY_SN$Obese_yn <- factor(PHENO.ANY_SN$Obese_yn, level = c("No", "Yes", "Unknown")) 

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])



####################
## Count Unknowns ##
####################
columns <- c("Current_smoker_yn", "PhysicalActivity_yn",
             "RiskyHeavyDrink_yn", "Obese_yn")
# Get the first letter of each column name
col_names <- substr(columns, 1, 1)

## Loop through each row of the dataframe
# PHENO.ANY_SN <- PHENO.ANY_SN[columns]
PHENO.ANY_SN$UNKNOWNS <- apply(PHENO.ANY_SN[columns], 1, function(x) {
  ifelse(sum(x == "Unknown") > 0, paste(col_names[x == "Unknown"], collapse = ""), NA)
})

table(PHENO.ANY_SN$UNKNOWNS)


library(purrr)
library(stringr)
library(gtools)
library(data.table)

## Loop over nchar>= 2
# life.vars <- unique(PHENO.ANY_SN$UNKNOWNS[which(nchar(PHENO.ANY_SN$UNKNOWNS) >= 2)])
elements <- c("C","O","P","R")
combinations <- lapply(2:length(elements), function(i) combn(elements, i, simplify = FALSE))
combinations <- unlist(combinations, recursive = FALSE)
combinations <- unique(combinations)
life.vars <- sapply(combinations, function(x) paste0(x, collapse = ""))

life.vars <- sapply(base::strsplit(life.vars, ""), function(x) paste(sort(x), collapse = ""))
life.vars
# "CO"   "CP"   "CR"   "OP"   "OR"   "PR"   "COP"  "COR"  "CPR"  "OPR"  "COPR"

cont.table <- {}
table(PHENO.ANY_SN$UNKNOWNS)

for (i in 1:length(life.vars)){
  PHENO.ANY_SN$UNKNOWNS.tmp <- gsub(paste0("[^",life.vars[i],"]"), "", PHENO.ANY_SN$UNKNOWNS)
  ## Sort the missing codes by letters
  PHENO.ANY_SN$UNKNOWNS.tmp <- sapply(base::strsplit(PHENO.ANY_SN$UNKNOWNS.tmp, ""), function(x) paste(sort(x), collapse = ""))
  table(PHENO.ANY_SN$UNKNOWNS.tmp)
  new_col_name <- paste0(life.vars[i], "_missing")
  PHENO.ANY_SN[[new_col_name]] <- as.factor(ifelse(grepl(life.vars[i], PHENO.ANY_SN$UNKNOWNS.tmp), "Yes", "No"))
  
  ## Get contingency table
  # CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]
  CROSS_CASES.df <- PHENO.ANY_SN[c("ANY_SN", new_col_name)]
  new_labels <- gsub("C", "Smoking", new_col_name)
  new_labels <- gsub("P", "Physical", new_labels)
  new_labels <- gsub("R", "Heavydrink", new_labels)
  new_labels <- gsub("O", "Obese", new_labels)
  new_labels <- gsub("(?<=[a-z])(?=[A-Z])", "_", new_labels)
  colnames(CROSS_CASES.df)[2] <- new_labels
  cont.table[[i]] <- table(CROSS_CASES.df[1:2])
}


# ###########################################
# ## Check data in each category/cross tab ##
# ###########################################
# library(expss)
# 
# # Getting counts for non-missing data only; 6 samples do not have admixture ancestry
# CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]
# 
# CROSS_CASES.df <- CROSS_CASES.df[c("ANY_SN", "Current_smoker_yn", "PhysicalActivity_yn",
#                                    "RiskyHeavyDrink_yn", "HEALTHY_Diet_yn", Obese_yn = "Obese_yn")]
# 
# CROSS_CASES.df <- apply_labels(CROSS_CASES.df, ANY_SN = "ANY_SN", 
#                                Current_smoker_yn = "Current_smoker_yn", PhysicalActivity_yn = "PhysicalActivity_yn",
#                                RiskyHeavyDrink_yn = "RiskyHeavyDrink_yn", HEALTHY_Diet_yn = "HEALTHY_Diet_yn", Obese_yn = "Obese_yn")
# 
# as.data.frame(t(CROSS_CASES.df %>%
#                         cross_cases(ANY_SN, list(Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Obese_yn))))
# 
# 
# cc <- as.data.frame(t(CROSS_CASES.df %>%
#                   cross_cases(ANY_SN, list(Current_smoker_yn, PhysicalActivity_yn, RiskyHeavyDrink_yn, HEALTHY_Diet_yn, Obese_yn))))
# 
# rownames(cc) <- NULL 
# # View(cc)

############################################################
## Drop Unknown level from the lifestyle factor variables ##
############################################################

## Recode lifestyle variables to fit the model for missingness
PHENO.ANY_SN$Current_smoker_yn[PHENO.ANY_SN$Current_smoker_yn == "Unknown"] <- "No"
PHENO.ANY_SN$Current_smoker_yn <- droplevels(PHENO.ANY_SN$Current_smoker_yn)

PHENO.ANY_SN$PhysicalActivity_yn[PHENO.ANY_SN$PhysicalActivity_yn == "Unknown"] <- "Yes"
PHENO.ANY_SN$PhysicalActivity_yn <- droplevels(PHENO.ANY_SN$PhysicalActivity_yn)

PHENO.ANY_SN$RiskyHeavyDrink_yn[PHENO.ANY_SN$RiskyHeavyDrink_yn == "Unknown"] <- "No"
PHENO.ANY_SN$RiskyHeavyDrink_yn <- droplevels(PHENO.ANY_SN$RiskyHeavyDrink_yn)

PHENO.ANY_SN$Obese_yn[PHENO.ANY_SN$Obese_yn == "Unknown"] <- "No"
PHENO.ANY_SN$Obese_yn <- droplevels(PHENO.ANY_SN$Obese_yn)




######################################
## Attributable fraction of Any SNs ##
######################################
# live.vars: "CPRO" "CPR"  "RO"   "CR"   "PR"
colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))] <- gsub("C", "Smoking", colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))])
colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))] <- gsub("P", "Physical", colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))])
colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))] <- gsub("R", "Heavydrink", colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))])
colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))] <- gsub("O", "Obese", colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))])
colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))] <- gsub("(?<=[a-z])(?=[A-Z])", "_", colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))])
colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))] 

# [1] "Smoking_Obese_missing"                     "Smoking_Physical_missing"                  "Smoking_Heavydrink_missing"               
# [4] "Obese_Physical_missing"                    "Obese_Heavydrink_missing"                  "Physical_Heavydrink_missing"              
# [7] "Smoking_Obese_Physical_missing"            "Smoking_Obese_Heavydrink_missing"          "Smoking_Physical_Heavydrink_missing"      
# [10] "Obese_Physical_Heavydrink_missing"         "Smoking_Obese_Physical_Heavydrink_missing"


models.types <- colnames(PHENO.ANY_SN)[grepl("_missing", colnames(PHENO.ANY_SN))] 

sn.model <- {}
for (i in 1:length(models.types)){
formula = paste0("ANY_SN ~ Pleiotropy_PRSWEB_PRS.tertile.category + AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category + Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn + EAS + AFR + ", models.types[i])

dat_all = PHENO.ANY_SN
fit_all = glm(formula = formula,
              family = binomial,
              data = dat_all)


summary(fit_all)

(output <- summary(fit_all)$coefficients)
as.data.frame(apply(output, 2, formatC, format="f", digits=4))
# options(scipen=999)
estimate <- format(round(output[,1],3), nsmall = 3)
std.error <- format(round(output[,2],3), nsmall = 3)
# P.val <- formatC(output[,4], format="G", digits=3)
P.val <- output[,4]
P.val[P.val < 0.001] <- "<0.001"
P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
sn.model[[i]] <- (setNames(cbind.data.frame(estimate, std.error, P.val
), c("Estimate", "Std.error", "P")))
sn.model[[i]] <- sn.model[[i]][!grepl("AGE_AT_LAST_CONTACT", row.names(sn.model[[i]])),]
}


ChangeNames <- function(x) {
  row.names(x)[grepl("PRS.tertile.category2nd", row.names(x))] <- "PRS.tertile.category2nd"
  row.names(x)[grepl("PRS.tertile.category3rd", row.names(x))] <- "PRS.tertile.category3rd"
  return(x)
}

## change PRS labels
df_list <- sn.model
df_list <- lapply(df_list, ChangeNames)

names(df_list) <- models.types
list2env(df_list ,.GlobalEnv)

View(Smoking_Obese_missing)
View(Smoking_Physical_missing)
View(Smoking_Heavydrink_missing) 
View(Obese_Physical_missing)
View(Obese_Heavydrink_missing)
View(Physical_Heavydrink_missing)
View(Smoking_Obese_Physical_missing)
View(Smoking_Obese_Heavydrink_missing)
View(Smoking_Physical_Heavydrink_missing)
View(Obese_Physical_Heavydrink_missing)
View(Smoking_Obese_Physical_Heavydrink_missing)


##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

###############
## Treatment ##
###############

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all

dat_tx$maxsegrtdose.category =
dat_tx$maxabdrtdose.category =
dat_tx$maxchestrtdose.category =
dat_tx$epitxn_dose_5.category = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_plp.prs$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$Current_smoker_yn = "No"
dat_lifestyle$PhysicalActivity_yn = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn = "No"
# dat_lifestyle$HEALTHY_Diet_yn = "Yes"
dat_lifestyle$Obese_yn = "No"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)

#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxsegrtdose.category =
  dat_tx.plp.prs.lifestyle$maxabdrtdose.category =
  dat_tx.plp.prs.lifestyle$maxchestrtdose.category =
  dat_tx.plp.prs.lifestyle$epitxn_dose_5.category = "None"

## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Pleiotropy_PRSWEB_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$Current_smoker_yn = "No"
dat_tx.plp.prs.lifestyle$PhysicalActivity_yn = "Yes"
dat_tx.plp.prs.lifestyle$RiskyHeavyDrink_yn = "No"
# dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn = "Yes"
dat_tx.plp.prs.lifestyle$Obese_yn = "No"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)

SN.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
SN.res
