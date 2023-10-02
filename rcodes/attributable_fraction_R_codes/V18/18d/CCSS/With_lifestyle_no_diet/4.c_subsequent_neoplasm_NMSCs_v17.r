# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression. **
obj_keep <- c("miss_Any_SN", "miss_SMN", "miss_NMSC", "miss_BREAST", "miss_THYROID", "miss_MENINGIOMA", "miss_SARCOMA")
rm(list = setdiff(ls(), obj_keep))
# rm(list=ls())
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

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

# ## Edit lifestyle variables
# source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
# PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)


data1 = read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/sns2022.sas7bdat")
data1=as.data.frame(data1)
# data1$ccssid <- paste0(data1$ccssid, "_", data1$ccssid)
data1$KEY <- paste0(data1$ccssid,":",data1$d_candx)

#####################
#####################
subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]

#########################
#########################
subneo$d_candx <- as.Date(subneo$d_candx, format = "%d%b%Y")
subneo$KEY <- paste0(subneo$ccssid,":",subneo$d_candx)
table(subneo$KEY %in% data1$KEY)
# FALSE  TRUE 
# 6307  3440 
table(data1$KEY %in% subneo$KEY)
# FALSE  TRUE 
# 4629  4434 
subneo$nmsc <- data1$nmsc[match(subneo$KEY, data1$KEY)]

# cc <- cbind.data.frame(subneo$KEY, subneo$nmsc, subneo$AGE.ANY_SN, subneo$groupdx3)

# Now get age of SN after first cancer
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
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)

# This will include basal cell, squamous cell and melanoma
NMSCs <- subneo[which((subneo$nmsc ==1| (subneo$nmsc == 2 & subneo$groupdx3 == "Skin"))),]
# cc <- cbind.data.frame(NMSC$KEY, NMSC$nmsc, NMSC$AGE.ANY_SN, NMSC$groupdx3)

# NMSCs <- subneo[grepl("skin", subneo$groupdx3, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(NMSCs)
# 775
dim(PHENO.ANY_SN)
# 7943   50

PHENO.ANY_SN$AGE.ANY_SN <- NMSCs$AGE.ANY_SN[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)]
if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
  PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
}

dim(PHENO.ANY_SN)

# Removing samples with SN within the 5 years of childhood cancer **
sum(PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid)
# 7

# "a_dx"  : Primary cancer diagnosis age
# "a_end" : age at last contact
# "d_candx" : Date when second cancer was diagnosed                                    
# "a_candx": Age at second cancer diagnosis


# dat[,c("ccssid","strokedt","event","dob","agelstcontact","agedx")]
NMSCs$gradeage <- NMSCs$AGE.ANY_SN
NMSCs$gradedt <- as.Date(NMSCs$d_candx, format = "%d%b%Y")
NMSCs$dob <- NMSCs$gradedt - as.numeric(NMSCs$gradeage) * 365.2422
PHENO.ANY_SN$dob <- NMSCs$dob[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)] ## 2009-02-12
PHENO.ANY_SN$gradedt <- NMSCs$gradedt[match(PHENO.ANY_SN$ccssid, NMSCs$ccssid)] ## 2009-02-12
###################

PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid,]
dim(PHENO.ANY_SN)
# 7912   51 ** END

PHENO.ANY_SN$NMSCs <- factor(ifelse(!PHENO.ANY_SN$ccssid %in% NMSCs$ccssid, 0, 1))
table(PHENO.ANY_SN$NMSCs)
# 0    1 
# 7149  763 


######################### **
############################################################################################
############################################################################################
#################################### Work for V11-4-v2 #####################################
############################################################################################
############################################################################################
# Yutaka's email 02/23/2023: The number of SN should decrease
# in the new analysis because those who developed SN before the first adult
# survey should be out of the analysis.  Also, we should use the "any missing in
# the lifestyle variables" rather than using the individual missing (with
# missing combined to the reference in each lifestyle variable).

PHENO.ANY_SN <- PHENO.ANY_SN[!(PHENO.ANY_SN$Current_smoker_yn == "Unknown" &
                                 PHENO.ANY_SN$PhysicalActivity_yn == "Unknown" &
                                 PHENO.ANY_SN$RiskyHeavyDrink_yn == "Unknown" &
                                 PHENO.ANY_SN$Obese_yn == "Unknown" ),]

dim(PHENO.ANY_SN)
# 7828   52

sum((PHENO.ANY_SN$Current_smoker_yn_agesurvey >= 18|
       PHENO.ANY_SN$PhysicalActivity_yn_agesurvey >= 18|
       PHENO.ANY_SN$RiskyHeavyDrink_yn_agesurvey >= 18|
       PHENO.ANY_SN$Obese_yn_agesurvey >= 18), na.rm = T)
# 7780


PHENO.ANY_SN <- PHENO.ANY_SN[which(PHENO.ANY_SN$Current_smoker_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$PhysicalActivity_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$RiskyHeavyDrink_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$Obese_yn_agesurvey >= 18),]

cols <- c(
  "Current_smoker_yn_agesurvey",
  "PhysicalActivity_yn_agesurvey",
  "RiskyHeavyDrink_yn_agesurvey",
  "Obese_yn_agesurvey"
)

# saved.cc <- ALL.LIFESTYLE[, c("SJLIFEID", cols)]
# ALL.LIFESTYLE[, cols] <- apply(PHENO.ANY_SN[, cols], 2, round)
library(matrixStats)
PHENO.ANY_SN$survey_min <- rowMins(as.matrix(PHENO.ANY_SN[, cols]), na.rm = TRUE)

test.1 <- PHENO.ANY_SN[c("ccssid", "survey_min", "Current_smoker_yn_agesurvey", "RiskyHeavyDrink_yn_agesurvey", "Obese_yn_agesurvey", "PhysicalActivity_yn_agesurvey","Current_smoker_yn", "RiskyHeavyDrink_yn", "Obese_yn", "PhysicalActivity_yn")]

# cc <- PHENO.ANY_SN[, c("SJLIFEID", cols, "survey_min")]
# PHENO.ANY_SN$Current_smoker_yn [which(PHENO.ANY_SN$Current_smoker_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"
# PHENO.ANY_SN$PhysicalActivity_yn [which(PHENO.ANY_SN$PhysicalActivity_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"
# PHENO.ANY_SN$RiskyHeavyDrink_yn [which(PHENO.ANY_SN$RiskyHeavyDrink_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"
# PHENO.ANY_SN$Obese_yn [which(PHENO.ANY_SN$Obese_yn_agesurvey != PHENO.ANY_SN$survey_min)] <- "Unknown"


PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$survey_min > PHENO.ANY_SN$AGE.ANY_SN),]
dim(PHENO.ANY_SN)
# 7775   53
######################### ** END


# PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn")]
# PHENO.ANY_SN$any_lifestyle_missing <- apply(PHENO.ANY_SN[c("Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn")], 1, function(x) any("Unknown" %in% x))
# PHENO.ANY_SN$any_lifestyle_missing  <- factor(ifelse(PHENO.ANY_SN$any_lifestyle_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_lifestyle_missing)
# No  Yes 
# 72 7703
########################################
########################################
# PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxpelvisrtdose.category")], 1, function(x) any("Unknown" %in% x))
# PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

table(PHENO.ANY_SN$any_tx_missing)
# No  Yes 
# 7131  514  
# PHENO.ANY_SN$any_rt_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxpelvisrtdose.category")], 1, function(x) any("Unknown" %in% x))
# PHENO.ANY_SN$any_rt_missing  <- factor(ifelse(PHENO.ANY_SN$any_rt_missing == FALSE, "No", "Yes"))

# PHENO.ANY_SN$any_tx_missing <- factor(PHENO.ANY_SN$any_tx_missing, levels = c("No", "Yes"))
# PHENO.ANY_SN$any_rt_missing <- factor(PHENO.ANY_SN$any_rt_missing, levels = c("No", "Yes"))
#########################
#########################
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
ethnicity.admixture$INDIVIDUAL <- sapply(strsplit(ethnicity.admixture$INDIVIDUAL,"_"), `[`, 1)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$ccssid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])

##############################
##############################
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/get_missing_combination_V17.R")
# Columns to check
columns_to_check <- c("PhysicalActivity_yn", "Current_smoker_yn", "RiskyHeavyDrink_yn", "Obese_yn")
# columns_to_check <- c("maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_5.category")
miss_NMSC <- get_missing_combinations(PHENO.ANY_SN, columns_to_check)

miss_NMSC <- calculate_missing_counts(PHENO.ANY_SN)
miss_NMSC <- calculate_missing_percentages(PHENO.ANY_SN)
############################################################
############################################################
# PHENO.ANY_SN$Current_smoker_yn[PHENO.ANY_SN$Current_smoker_yn == "Unknown"] <- "No"
# PHENO.ANY_SN$Current_smoker_yn <- droplevels(PHENO.ANY_SN$Current_smoker_yn)

# PHENO.ANY_SN$PhysicalActivity_yn[PHENO.ANY_SN$PhysicalActivity_yn == "Unknown"] <- "Yes"
# PHENO.ANY_SN$PhysicalActivity_yn <- droplevels(PHENO.ANY_SN$PhysicalActivity_yn)

# PHENO.ANY_SN$RiskyHeavyDrink_yn[PHENO.ANY_SN$RiskyHeavyDrink_yn == "Unknown"] <- "No"
# PHENO.ANY_SN$RiskyHeavyDrink_yn <- droplevels(PHENO.ANY_SN$RiskyHeavyDrink_yn)

# PHENO.ANY_SN$Obese_yn[PHENO.ANY_SN$Obese_yn == "Unknown"] <- "No"
# PHENO.ANY_SN$Obese_yn <- droplevels(PHENO.ANY_SN$Obese_yn)


# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxsegrtdose.category <- droplevels(PHENO.ANY_SN$maxsegrtdose.category)

# PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxneckrtdose.category <- droplevels(PHENO.ANY_SN$maxneckrtdose.category)

# PHENO.ANY_SN$maxabdrtdose.category[PHENO.ANY_SN$maxabdrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxabdrtdose.category <- droplevels(PHENO.ANY_SN$maxabdrtdose.category)

# PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxchestrtdose.category <- droplevels(PHENO.ANY_SN$maxchestrtdose.category)

# PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == "Unknown"] <- "None"
# PHENO.ANY_SN$maxpelvisrtdose.category <- droplevels(PHENO.ANY_SN$maxpelvisrtdose.category)

# PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "Unknown"] <- "None"
# PHENO.ANY_SN$epitxn_dose_5.category <- droplevels(PHENO.ANY_SN$epitxn_dose_5.category)

# PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "Unknown"] <- "None"
# PHENO.ANY_SN$anthra_jco_dose_5.category <- droplevels(PHENO.ANY_SN$anthra_jco_dose_5.category)

# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "Unknown"] <- "None"
# PHENO.ANY_SN$aa_class_dose_5.category <- droplevels(PHENO.ANY_SN$aa_class_dose_5.category)


################
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
########################################  ## change agedx to survey_min which is the age at first adult survey
data <- PHENO.ANY_SN

data$agelstcontact[!is.na(data$AGE.ANY_SN)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]
# If age at adult survey is greater than age at last contact for Controls, Use age at survey as age at last contact **
BOOL <- is.na(data$AGE.ANY_SN) & data$agelstcontact < data$survey_min
data$agelstcontact[BOOL] <- data$survey_min[BOOL]
# data$ccssid[BOOL]


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
### For those without event, start is first age at adult survey and end is Fu date === for SN, analysis starts from age at first adult survey
alldata$start[alldata$event==0] <- alldata$survey_min[alldata$event==0]
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is survey_min and end is first event time
alldata$start[alldata$event==1 & alldata$first==1] <- alldata$survey_min[alldata$event==1 & alldata$first==1] 
alldata$end[alldata$event==1 & alldata$first==1] <- as.numeric(difftime(alldata$gradedt[alldata$event==1 & alldata$first==1],alldata$dob[alldata$event==1 & alldata$first==1], units = 'days')/365.25)

#### For events that are not the first, segments are from the previous event date to this event date
alldata$previous <- as.Date(c(NA,alldata$gradedt[1:length(alldata$gradedt)-1]),origin = "1970-01-01")
alldata[1:10,]
### if first>1 (i.e, 2 or more events, previous event time remained, others are missing)
alldata$previous[alldata$first==1] <- NA
alldata$start[alldata$first>1] <- as.numeric(difftime(alldata$previous[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)
alldata$end[alldata$first>1] <- as.numeric(difftime(alldata$gradedt[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)


### If one person has only 1 event, we need to have segment from survey_min to event date (handled above), and then from event date to end of FU
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

final$ccssid[duplicated(final$ccssid)]


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
###############
PHENO.ANY_SN <- SNs_py[c("ccssid", "event", "Pleiotropy_PRSWEB_PRS.tertile.category", "BASALcell_PRS.tertile.category", 
                         "Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category", "Thyroid_PRS.tertile.category",
                         "Meningioma_PRS.tertile.category", "Sarcoma_Machiela_PRS.tertile.category",
                         "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4",
                         "AGE_AT_DIAGNOSIS", "gender", 
                         "maxsegrtdose.category", "maxneckrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category",
                         "maxpelvisrtdose.category", "epitxn_dose_5.category", "anthra_jco_dose_5.category", "aa_class_dose_5.category",
                         "EAS", "AFR", 
                         "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn", 
                         "PY","evt1", "end")]



source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/cubic_spline.r")
breaks = seq(5, 95, 22.5)
cp = quantile(PHENO.ANY_SN$end, breaks/100, na.rm = T)
cs = cubic_spline(PHENO.ANY_SN$end, knots = cp)
PHENO.ANY_SN <- PHENO.ANY_SN[!colnames(PHENO.ANY_SN) %in% c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")]
colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, cs)


rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.NMSCs.V18d_without_diet.Rdata")

