# Yutakas email on 03/08/2023: I think this result looks good to me.  The only
# thing is that large (positive or negative) coefficients are clustered in
# "Unknown" treatment categories.  This is presumably due to the same reason as
# the lifestyle unknowns where "unknowns" of different items are overlapping.
# In the case of treatment, if the survivor's medical record abstraction is
# missing, all treatment variables will be unknown.  Thus, I suggest you do the
# same for the treatment "Unknowns" as you did lifestyle "Unknowns", i.e.,
# create an indicator variable "Any tx Unknown" and put everyone with "Unknown"
# at any of the treatment variables there.  Then, you will remove "Unknown" from
# each of the treatment variables and we will not have the large coefficients (I
# think).
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
# Get THYROIDcancer for the first time and Age at First THYROIDcancer.
# For this, I will first sort the table by date
library(data.table)
THYROIDcancer <- subneo[grepl("thyroid", subneo$diag, ignore.case = T),]
THYROIDcancer <- setDT(THYROIDcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
THYROIDcancer <- THYROIDcancer[!THYROIDcancer$sjlid %in% subneo.within5$sjlid,]
nrow(THYROIDcancer)
# 86

## Remove THYROIDcancer if younger than 18
PHENO.ANY_SN$AGE.ANY_SN <- THYROIDcancer$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, THYROIDcancer$sjlid)]
if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
}

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]

PHENO.ANY_SN$THYROIDcancer <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% THYROIDcancer$sjlid, 0, 1))

table(PHENO.ANY_SN$THYROIDcancer)
# 0    1 
# 4293 81

########################################
## Do the same for missing treatments ##
########################################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxneckrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


############################################################
## Drop Unknown level from the lifestyle factor variables ##
############################################################
## Missing tx
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxsegrtdose.category <- droplevels(PHENO.ANY_SN$maxsegrtdose.category)

PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxneckrtdose.category <- droplevels(PHENO.ANY_SN$maxneckrtdose.category)

PHENO.ANY_SN$maxabdrtdose.category[PHENO.ANY_SN$maxabdrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxabdrtdose.category <- droplevels(PHENO.ANY_SN$maxabdrtdose.category)

PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxchestrtdose.category <- droplevels(PHENO.ANY_SN$maxchestrtdose.category)

PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxpelvisrtdose.category <- droplevels(PHENO.ANY_SN$maxpelvisrtdose.category)

PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$epitxn_dose_5.category <- droplevels(PHENO.ANY_SN$epitxn_dose_5.category)

PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$anthra_jco_dose_5.category <- droplevels(PHENO.ANY_SN$anthra_jco_dose_5.category)

PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$aa_class_dose_5.category <- droplevels(PHENO.ANY_SN$aa_class_dose_5.category)


################
## Cross tabs ##
################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- PHENO.ANY_SN

CROSS_CASES.df <- CROSS_CASES.df[,c("THYROIDcancer", "maxneckrtdose.category", "epitxn_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, THYROIDcancer = "THYROIDcancer", 
                               maxneckrtdose.category = "maxneckrtdose.category", epitxn_dose_5.category = "epitxn_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(THYROIDcancer, list(maxneckrtdose.category, epitxn_dose_5.category))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(THYROIDcancer, list(maxneckrtdose.category, epitxn_dose_5.category))))

rownames(cc) <- NULL 
cc


# Create a cross-tabulation table between maxsegrtdose and maxchedtrtdose for cases
cases_table <- table(Max_SegmentedRT_Dose = PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$THYROIDcancer == 1],
                     Max_ChestRT_Dose = PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$THYROIDcancer == 1])

# Create a cross-tabulation table between maxsegrtdose and maxchedtrtdose for controls
control_table <- table(Max_SegmentedRT_Dose = PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$THYROIDcancer == 0],
                     Max_ChestRT_Dose = PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$THYROIDcancer == 0])


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.THYROIDcancer.V14-4-3.Rdata")
