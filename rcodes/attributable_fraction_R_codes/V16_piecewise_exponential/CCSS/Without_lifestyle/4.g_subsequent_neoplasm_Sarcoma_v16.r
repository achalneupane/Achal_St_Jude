rm()
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")

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
library(survival)

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$ccssid))
# 26


#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
SARCOMA <- subneo[grepl("sarcoma", subneo$groupdx3, ignore.case = T),]
SARCOMA <- setDT(SARCOMA)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(SARCOMA)
# 109

# based on Yadav's email on 03/09/2023, I am removing all benign diagnoses from the list of 52 survivors
to.remove <- as.character(c(1004973, 1005736, 3012171, 4073492, 5097496, 5146972, 8217873, 9059523, 9203577,
                                           10085746, 11108731, 12083337, 13054941, 13231652, 16041746, 16045012, 17050333,
                                           18080902, 18141511, 20024771, 20027745, 20032881, 20033541, 21228953, 22091488,
                                           22155815, 22156111, 22200376, 25017727, 26016681, 26018907, 26020735, 26056273,
                                           1262696, 2511092, 2518314, 5362062, 6302298, 8356277, 15283414, 19295502, 22434302,
                                           26403512, 27117943))
                                           


## remove samples
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% to.remove,]

## Remove SNs if younger than 18 **
dim(PHENO.ANY_SN)
# 7943   50

PHENO.ANY_SN$AGE.ANY_SN <- SARCOMA$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, SARCOMA$ccssid)]
# if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
#   PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
# }


# "a_dx"  : Primary cancer diagnosis age
# "a_end" : age at last contact
# "d_candx" : Date when second cancer was diagnosed                                    
# "a_candx": Age at second cancer diagnosis


# dat[,c("ccssid","strokedt","event","dob","agelstcontact","agedx")]
SARCOMA$gradeage <- SARCOMA$gradedt
SARCOMA$gradedt <- as.Date(SARCOMA$d_candx, format = "%d%b%Y")
## Calculate DOB
SARCOMA$dob <- SARCOMA$gradedt - as.numeric(SARCOMA$gradeage) * 365.2422
PHENO.ANY_SN$dob <- SARCOMA$dob[match(PHENO.ANY_SN$ccssid, SARCOMA$ccssid)] ## 2009-02-12
PHENO.ANY_SN$gradedt <- SARCOMA$gradedt[match(PHENO.ANY_SN$ccssid, SARCOMA$ccssid)] ## 2009-02-12


dim(PHENO.ANY_SN)
## 7928 52 ** END

# Removing samples with SN within the 5 years of childhood cancer **
sum(PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid,]
dim(PHENO.ANY_SN)
# 7906 ** END

## CA CO status
PHENO.ANY_SN$SARCOMA <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% SARCOMA$ccssid, 0, 1))
table(PHENO.ANY_SN$SARCOMA)
# 0    1 
# 7814   92


######################### **

#################
## Missingness ##
#################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("aa_class_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_tx_missing)
# No  Yes 
# 7255  651
PHENO.ANY_SN$any_chemo_missing <- apply(PHENO.ANY_SN[c("aa_class_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_chemo_missing  <- factor(ifelse(PHENO.ANY_SN$any_chemo_missing == FALSE, "No", "Yes"))
#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
ethnicity.admixture$INDIVIDUAL <- sapply(strsplit(ethnicity.admixture$INDIVIDUAL,"_"), `[`, 1)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


############################################################
## Drop Unknown level from the lifestyle factor variables ##
############################################################
## Recode tx variables to fit the model for missingness
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

# ###############################################
# ## Find out benign Sarcoma's and remove them ##
# ###############################################
# # # based on Yadav's email on 03/09/2023, I am removing all benign diagnoses from the list of 52 survivors
# ## This file is from Kyla
# KIRI.ccss <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Kyla/combinedsn_final_02_17_2023.csv", header = T, sep = ",", stringsAsFactors = F)
# dim(KIRI.ccss)
# ## Keep non-missing candxo3
# KIRI.ccss <- KIRI.ccss[!is.na(KIRI.ccss$candxo3),]
# # KIRI.ccss <- KIRI.ccss[KIRI.ccss$candxo3 !="",]
# KIRI.ccss <- KIRI.ccss[KIRI.ccss$d_candx !="",]
# dim(KIRI.ccss)
# KIRI.ccss$SN_diagnosis_date <- as.Date(KIRI.ccss$d_candx, format = "%d%b%Y")
# KIRI.ccss$SN_diagnosis_date <- format(KIRI.ccss$SN_diagnosis_date, "%m-%d-%Y") # 06-30-2008
# KIRI.ccss$KEY <- paste(KIRI.ccss$ccssid, KIRI.ccss$SN_diagnosis_date, sep = ":")
# 
# PHENO.ANY_SN$SN_diagnosis_date <- as.Date(PHENO.ANY_SN$d_candx, format = "%d%b%Y")
# PHENO.ANY_SN$SN_diagnosis_date  <- format(PHENO.ANY_SN$SN_diagnosis_date, "%m-%d-%Y") # "08-23-2016"
# PHENO.ANY_SN$ccssid <- gsub("_.*","",PHENO.ANY_SN$ccssid)
# PHENO.ANY_SN$KEY <- paste(PHENO.ANY_SN$ccssid, PHENO.ANY_SN$SN_diagnosis_date, sep = ":")
# 
# 
# table(PHENO.ANY_SN$KEY %in% KIRI.ccss$KEY)
# # FALSE  TRUE 
# # 6326  1580 



# SARCOMA <- PHENO.ANY_SN[PHENO.ANY_SN$SARCOMA == 1,]
# KIRI.ccss <- KIRI.ccss[KIRI.ccss$KEY %in% SARCOMA$KEY,]
# KIRI.save <- KIRI.ccss
# KIRI.ccss <- KIRI.ccss[grepl("\\/0|\\/1", KIRI.ccss$candxo3),]
# 
# PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% KIRI.ccss$ccssid,]
table(PHENO.ANY_SN$SARCOMA)
# 7814   53
##########################################

################
## Cross tabs ##
################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- PHENO.ANY_SN

CROSS_CASES.df <- CROSS_CASES.df[,c("SARCOMA", "aa_class_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, SARCOMA = "SARCOMA", 
                               aa_class_dose_5.category = "aa_class_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(SARCOMA, list(aa_class_dose_5.category))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(SARCOMA, list(aa_class_dose_5.category))))

rownames(cc) <- NULL 
cc


########################################
## Prepare data accoding to Qi's code ## 
########################################
data <- PHENO.ANY_SN

## Age at last contact for cases is SN diagnosis data
data$agelstcontact[!is.na(data$AGE.ANY_SN)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

data$event <- ifelse(!is.na(data$gradedt), 1, 0)

data$first <- ave(data$agelstcontact, data$ccssid, FUN = seq_along)
M <- max(data$first, na.rm = T)  ### maximum number of events
data[data$first==M,]$ccssid  #the id for the person with the maximum number of rows.
data[data$ccssid==data[data$first==M,]$ccssid,]

### Take the last row so we know the maximum number of events per person
event.number <- do.call(rbind, lapply(split(data, data$ccssid), tail, 1))[,c("ccssid","first")]
colnames(event.number) <- c("ccssid","maxE")

alldata <- merge(data,event.number,by.x="ccssid",by.y="ccssid")

alldata$start <- NULL
alldata$end <- NULL
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
### Achal: Since the analysis start from 5 years post DX, the above line has been revised to: alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0]+5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx and end is first event time
alldata$start[alldata$event==1 & alldata$first==1] <- alldata$agedx[alldata$event==1 & alldata$first==1] + 5
### Achal: also added +5 in the above line
alldata$end[alldata$event==1 & alldata$first==1] <- as.numeric(difftime(alldata$gradedt[alldata$event==1 & alldata$first==1],alldata$dob[alldata$event==1 & alldata$first==1], units = 'days')/365.25)

#### For events that are not the first, segments are from the previous event date to this event date
alldata$previous <- as.Date(c(NA,alldata$gradedt[1:length(alldata$gradedt)-1]),origin = "1970-01-01")
alldata[1:10,]
### if first>1 (i.e, 2 or more events, previous event time remained, others are missing)
alldata$previous[alldata$first==1] <- NA
alldata$start[alldata$first>1] <- as.numeric(difftime(alldata$previous[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)
alldata$end[alldata$first>1] <- as.numeric(difftime(alldata$gradedt[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)


### If one person has only 1 event, we need to have segment from agedx to event date (handled above), and then from event date to end of FU
### If one person has multiple events, we need to add segments from the last event to end of Fu.
row_add <- alldata[alldata$event==1 & alldata$first==alldata$maxE,] ## this includes (1) people with only 1 event (2) The last row for people with multiple events
row_add$start <- as.numeric(difftime(row_add$gradedt, row_add$dob, units = 'days')/365.25)
row_add$end <- row_add$agelstcontact;
#### remember to make these rows event as 0, as we are adding the time after the last event date.
row_add$event <- 0
### If we do first event analysis, we would not need these. I need to add an indicator here, in case I only need to get the first event analysis segments.
row_add$evt1 <- 0;

alldata$evt1 <- 0;
alldata$evt1[alldata$first==1] <- 1
adata <- rbind(alldata,row_add)
#### order by person and start date
adata <- adata[order(adata$ccssid, adata$start, decreasing = FALSE),]
table(adata$event)## Double check event numebr is correct

###any stop time <=start time
adata$end[adata$end<=adata$start] <- adata$end[adata$end<=adata$start] + 1/365
any <- adata[adata$end<=adata$start,] # 
### 2 people had the stroke date on the Fu date. In this case, either get rid of the 2 lines [i.e, no time is follow-up after the last event date], or add 1 day on end date of these 2 segments, assuming there were followed up 1 more day. Will not make much difference. 1 day out of 365 days is 0.0027.
diff=any$start-any$end ###Qi: These are people who had SN after the last contact date. Just wonder why this could happen. While it may not make the results differ, I wonder is there any reason to keep the events but change last contact date to be SN+1day? Depends on why there are SN after last contact date.
dim(adata)
final <- adata[adata$end>adata$start,]

minimum  <-  min(final$start, na.rm = TRUE) - 1
if (minimum < 0) minimum <- 0
maximum  <-  max(final$end, na.rm = TRUE) + 1
SNs_py <-  survSplit(final, cut=seq(minimum, maximum, 1), end="end",start="start",event="event") 
table(SNs_py$event)   ## Double check event numebr is correct
#### If you need the rows for first event analysis, take evt1=1
length(unique(SNs_py$ccssid[SNs_py$event==1]))
length(unique(SNs_py$ccssid[SNs_py$event==1 & SNs_py$evt1==1 ]))

SNs_py$PY <- SNs_py$end-SNs_py$start

###############
## Model fit ##
###############
PHENO.ANY_SN <- SNs_py[c("ccssid", "event", "Pleiotropy_PRSWEB_PRS.tertile.category", "BASALcell_PRS.tertile.category", 
                         "Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category", "Thyroid_PRS.tertile.category",
                         "Meningioma_PRS.tertile.category", "Sarcoma_Machiela_PRS.tertile.category",
                         "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4",
                         "AGE_AT_DIAGNOSIS", "gender", 
                         "maxsegrtdose.category", "maxneckrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category",
                         "maxpelvisrtdose.category", "epitxn_dose_5.category", "anthra_jco_dose_5.category", "aa_class_dose_5.category",
                         "EAS", "AFR", 
                         "any_tx_missing", "any_chemo_missing",
                         "PY","evt1")]


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.SARCOMA_without_lifestyle.V16.Rdata")

