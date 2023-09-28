# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
rm(list=ls())
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
library(survival)

# subneo <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/subneo.sas7bdat")
# subneo <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_07_22_2023_subneo.txt", header = T, stringsAsFactors = F)
subneo <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/subneoplasms.sas7bdat')

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

subneo$gradedt <- as.Date(subneo$gradedt, "%m/%d/%Y") ## **
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
# Get NMSCs for the first time and Age at First NMSCs.
# For this, I will first sort the table by date
library(data.table)
## Remove SNs as cases that are within 5 years of primary diagnosis
subneo <- subneo[!subneo$sjlid %in% subneo.within5$sjlid,]

NMSCs <- subneo[grepl("basal cell|squamous cell", subneo$diag, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


nrow(NMSCs)
# 249

## Remove NMSCs if younger than 18
PHENO.ANY_SN$AGE.ANY_SN <- NMSCs$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, NMSCs$sjlid)]
PHENO.ANY_SN$gradedt <- NMSCs$gradedt[match(PHENO.ANY_SN$sjlid, NMSCs$sjlid)]

# if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
# PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
# }

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]

PHENO.ANY_SN$NMSCs <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% NMSCs$sjlid, 0, 1))
table(PHENO.ANY_SN$NMSCs)
# 0    1 
# 4130  249

# #################
# ## missingness ##
# #################
# PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxpelvisrtdose.category")], 1, function(x) any("Unknown" %in% x))
# PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))
# 
# PHENO.ANY_SN$any_rt_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxpelvisrtdose.category")], 1, function(x) any("Unknown" %in% x))
# PHENO.ANY_SN$any_rt_missing  <- factor(ifelse(PHENO.ANY_SN$any_rt_missing == FALSE, "No", "Yes"))

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


# ############################################################
# ## Drop Unknown level from the lifestyle factor variables ##
# ############################################################
# ## Missing tx
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxsegrtdose.category <- droplevels(PHENO.ANY_SN$maxsegrtdose.category)
# 
# PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxneckrtdose.category <- droplevels(PHENO.ANY_SN$maxneckrtdose.category)
# 
# PHENO.ANY_SN$maxabdrtdose.category[PHENO.ANY_SN$maxabdrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxabdrtdose.category <- droplevels(PHENO.ANY_SN$maxabdrtdose.category)
# 
# PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxchestrtdose.category <- droplevels(PHENO.ANY_SN$maxchestrtdose.category)
# 
# PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxpelvisrtdose.category <- droplevels(PHENO.ANY_SN$maxpelvisrtdose.category)
# 
# PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "Unknown"] <- "None"
# PHENO.ANY_SN$epitxn_dose_5.category <- droplevels(PHENO.ANY_SN$epitxn_dose_5.category)
# 
# PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "Unknown"] <- "None"
# PHENO.ANY_SN$anthra_jco_dose_5.category <- droplevels(PHENO.ANY_SN$anthra_jco_dose_5.category)
# 
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "Unknown"] <- "None"
# PHENO.ANY_SN$aa_class_dose_5.category <- droplevels(PHENO.ANY_SN$aa_class_dose_5.category)




################
## Cross tabs ##
################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- PHENO.ANY_SN

CROSS_CASES.df <- CROSS_CASES.df[,c("NMSCs", "maxsegrtdose.category", "maxabdrtdose.category", "maxpelvisrtdose.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, NMSCs = "NMSCs", 
                               maxsegrtdose.category = "maxsegrtdose.category", maxabdrtdose.category="maxabdrtdose.category", maxpelvisrtdose.category = "maxpelvisrtdose.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(NMSCs, list(maxsegrtdose.category, maxabdrtdose.category, maxpelvisrtdose.category))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(NMSCs, list(maxsegrtdose.category, maxabdrtdose.category, maxpelvisrtdose.category))))

rownames(cc) <- NULL 
cc


########################################
## Prepare data accoding to Qi's code ## 
########################################
data <- PHENO.ANY_SN[!grepl("MRN", colnames(PHENO.ANY_SN))]

## Age at last contact for cases is SN diagnosis data
data$agelstcontact[!is.na(data$AGE.ANY_SN)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

data$event <- ifelse(!is.na(data$gradedt), 1, 0)

data$first <- ave(data$agelstcontact, data$sjlid, FUN = seq_along)
M <- max(data$first, na.rm = T)  ### maximum number of events
data[data$first==M,]$sjlid  #the id for the person with the maximum number of rows.
data[data$sjlid==data[data$first==M,]$sjlid,]

### Take the last row so we know the maximum number of events per person
event.number <- do.call(rbind, lapply(split(data, data$sjlid), tail, 1))[,c("sjlid","first")]
colnames(event.number) <- c("sjlid","maxE")

alldata <- merge(data,event.number,by.x="sjlid",by.y="sjlid")

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
adata <- adata[order(adata$sjlid, adata$start, decreasing = FALSE),]
table(adata$event)## Double check event numebr is correct

###any stop time <=start time
adata$end[adata$end<=adata$start] <- adata$end[adata$end<=adata$start] + 1/365
any <- adata[adata$end<=adata$start,] # 
diff=any$start-any$end 
dim(adata)
final <- adata[adata$end>adata$start,]

minimum  <-  min(final$start, na.rm = TRUE) - 1
if (minimum < 0) minimum <- 0
maximum  <-  max(final$end, na.rm = TRUE) + 1
SNs_py <-  survSplit(final, cut=seq(minimum, maximum, 1), end="end",start="start",event="event") 
table(SNs_py$event)   ## Double check event numebr is correct
#### If you need the rows for first event analysis, take evt1=1
length(unique(SNs_py$sjlid[SNs_py$event==1]))
length(unique(SNs_py$sjlid[SNs_py$event==1 & SNs_py$evt1==1 ]))

SNs_py$PY <- SNs_py$end-SNs_py$start

###############
## Model fit ##
###############
PHENO.ANY_SN <- SNs_py[c("sjlid", "event", "Pleiotropy_PRSWEB_PRS.tertile.category", "BASALcell_PRS.tertile.category", 
                         "Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category", "Thyroid_PRS.tertile.category",
                         "Meningioma_PRS.tertile.category", "Sarcoma_Machiela_PRS.tertile.category",
                         "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4",
                         "AGE_AT_DIAGNOSIS", "gender", 
                         "maxsegrtdose.category", "maxneckrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category",
                         "maxpelvisrtdose.category", "epitxn_dose_5.category", "anthra_jco_dose_5.category", "aa_class_dose_5.category",
                         "EAS", "AFR", 
                         "PY","evt1", "end")]

## Age attained age (cubic spline)
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/cubic_spline.r")
breaks = seq(5, 95, 22.5)
cp = quantile(PHENO.ANY_SN$end, breaks/100, na.rm = T)
cs = cubic_spline(PHENO.ANY_SN$end, knots = cp)
PHENO.ANY_SN <- PHENO.ANY_SN[!colnames(PHENO.ANY_SN) %in% c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")]
colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, cs)


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.NMSCs.V18.Rdata")
