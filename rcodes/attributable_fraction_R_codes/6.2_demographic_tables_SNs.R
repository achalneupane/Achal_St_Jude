#################################
## 1....................SJLIFE ##
#################################
#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")

demo <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat')
demo <- demo[demo$sjlid %in% PHENO.ANY_SN$sjlid,]
PHENO.ANY_SN$ethnic <- demo$hispanic[match(PHENO.ANY_SN$sjlid, demo$sjlid)]


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
# Get SNs for the first time and Age at First SN. Note: I am only keeping first event SN 5 years post diagnosis
# For this, I will first sort the table by date
library(data.table)

## Remove SNs as cases that are within 5 years of primary diagnosis
subneo <- subneo[!subneo$sjlid %in% subneo.within5$sjlid,]
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


PHENO.ANY_SN$gradedt <- ANY_SNs$gradedt[match(PHENO.ANY_SN$sjlid, ANY_SNs$sjlid)]
PHENO.ANY_SN$AGE.ANY_SN <- ANY_SNs$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, ANY_SNs$sjlid)]

# ## Remove SNs if younger than 18 **
# if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
# PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
# }

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]


PHENO.ANY_SN$ANY_SNs <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))

table(PHENO.ANY_SN$ANY_SNs)

get.SNs <- PHENO.ANY_SN[PHENO.ANY_SN$ANY_SNs == 1,]

## SMN
# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)
#########################
## Subsequent Neoplasm ##
#########################

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

subneo$AGE.ANY_SN.after.childhood.cancer <- time_length(interval(as.Date(subneo$diagdt), as.Date(subneo$gradedt)), "years")
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$sjlid))
# 612

#############
## Any SNs ##
#############
# Get SMNs for the first time and Age at First SMNs.
# For this, I will first sort the table by date
library(data.table)

## Keeping only malignant based on ICDO3 behaviour
subneo <- subneo[grepl("^Malignant", subneo$icdo3behavior, ignore.case = T),]

SMNs <- subneo

# SN within 5 years
subneo.within5 <- SMNs[SMNs$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 18

# SMNs <- subneo[!grepl("basal cell|Basosquamous|squamous cell", subneo$diag, ignore.case = T),]
SMNs <- setDT(SMNs)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
SMNs <- SMNs[!SMNs$sjlid %in% subneo.within5$sjlid,]
nrow(SMNs)

table(SMNs$icdo3behavior)

PHENO.ANY_SN$gradedt <- SMNs$gradedt[match(PHENO.ANY_SN$sjlid, SMNs$sjlid)]
PHENO.ANY_SN$AGE.ANY_SN <- SMNs$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, SMNs$sjlid)]


## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 18
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]

PHENO.ANY_SN$SMNs <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% SMNs$sjlid, 0, 1))


table(PHENO.ANY_SN$SMNs)
# 0    1 
# 3920  463



# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
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

#############
## Any SNs ##
#############
# Get NMSCs for the first time and Age at First NMSCs.
# For this, I will first sort the table by date
library(data.table)
## Remove SNs as cases that are within 5 years of primary diagnosis
NMSCs <- subneo[grepl("basal cell|squamous cell", subneo$diag, ignore.case = T),]

# SN within 5 years
subneo.within5 <- NMSCs[NMSCs$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 1
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


nrow(NMSCs)



load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)
#########################
## Subsequent Neoplasm ##
#########################

library(haven)

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

#############
## Any SNs ##
#############
# Get BREASTcancer for the first time and Age at First BREASTcancer.
# For this, I will first sort the table by date
library(data.table)

BREASTcancer <- subneo[grepl("breast", subneo$diag, ignore.case = T),]

## SN within 5 years
subneo.within5 <- BREASTcancer[BREASTcancer$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

BREASTcancer <- setDT(BREASTcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
BREASTcancer <- BREASTcancer[!BREASTcancer$sjlid %in% subneo.within5$sjlid,]
nrow(BREASTcancer)
# 78


## Remove BREASTcancer if younger than 18
PHENO.ANY_SN$gradedt <- BREASTcancer$gradedt[match(PHENO.ANY_SN$sjlid, BREASTcancer$sjlid)]
PHENO.ANY_SN$AGE.ANY_SN <- BREASTcancer$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, BREASTcancer$sjlid)]
# if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
# PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
# }

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 0
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]


# ## Keep females only
# BREASTcancer$gender <-  PHENO.ANY_SN$gender[match(BREASTcancer$sjlid , PHENO.ANY_SN$sjlid)]
# BREASTcancer <- BREASTcancer[BREASTcancer$gender == "Female",]

# ## Remove those that are not breast cancer
# subneo.not.breast.cancer <- unique(subneo$sjlid[!subneo$sjlid %in% unique(BREASTcancer$sjlid)])
# PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.not.breast.cancer,]


PHENO.ANY_SN$BREASTcancer <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% BREASTcancer$sjlid, 0, 1))
PHENO.ANY_SN <- PHENO.ANY_SN[PHENO.ANY_SN$gender == "Female",]

table(PHENO.ANY_SN$BREASTcancer)

BREASTcancer <- BREASTcancer[BREASTcancer$sjlid %in% PHENO.ANY_SN$sjlid,]

#########################
## Load Phenotype data ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)
#########################
## Subsequent Neoplasm ##
#########################


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

subneo$AGE.ANY_SN.after.childhood.cancer <- time_length(interval(as.Date(subneo$diagdt), as.Date(subneo$gradedt)), "years")
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$sjlid))
# 612

#############
## Any SNs ##
#############
# Get THYROIDcancer for the first time and Age at First THYROIDcancer.
# For this, I will first sort the table by date
library(data.table)

THYROIDcancer <- subneo[grepl("thyroid", subneo$diag, ignore.case = T),]

## SN within 5 years
subneo.within5 <- THYROIDcancer[THYROIDcancer$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 1

THYROIDcancer <- setDT(THYROIDcancer)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
THYROIDcancer <- THYROIDcancer[!THYROIDcancer$sjlid %in% subneo.within5$sjlid,]
nrow(THYROIDcancer)
# 87


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

################
## Meningioma ##
################
# Get MENINGIOMA for the first time and Age at First MENINGIOMA.
# For this, I will first sort the table by date
library(data.table)

MENINGIOMA <- subneo[grepl("meningioma", subneo$diag, ignore.case = T),]

## SN within 5 years
subneo.within5 <- MENINGIOMA[MENINGIOMA$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

MENINGIOMA <- setDT(MENINGIOMA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
MENINGIOMA <- MENINGIOMA[!MENINGIOMA$sjlid %in% subneo.within5$sjlid,]
nrow(MENINGIOMA)
# 149


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

################
## Meningioma ##
################
# Get SARCOMA for the first time and Age at First SARCOMA.
# For this, I will first sort the table by date
library(data.table)

SARCOMA <- subneo[grepl("sarcoma", subneo$diag, ignore.case = T),]

## SN within 5 years
subneo.within5 <- SARCOMA[SARCOMA$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))

SARCOMA <- setDT(SARCOMA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

## Remove SNs as cases that are within 5 years of primary diagnosis
SARCOMA <- SARCOMA[!SARCOMA$sjlid %in% subneo.within5$sjlid,]
nrow(SARCOMA)
# 33




# all.malignants <- unique(c(as.character(SMNs$sjlid), as.character(NMSCs$sjlid), as.character(BREASTcancer$sjlid),
#                     as.character(THYROIDcancer$sjlid), as.character(MENINGIOMA$sjlid),
#                     as.character(SARCOMA$sjlid)))

all.malignants <- unique(c(as.character(NMSCs$sjlid), as.character(BREASTcancer$sjlid), 
                           as.character(THYROIDcancer$sjlid),
                           as.character(SARCOMA$sjlid)))

all.malignants <- unique(c(as.character(NMSCs$sjlid), as.character(BREASTcancer$sjlid), 
                           as.character(THYROIDcancer$sjlid),as.character(MENINGIOMA$sjlid),
                           as.character(SARCOMA$sjlid)))


length(all.malignants)
## 591

unique <- SMNs$sjlid[(!SMNs$sjlid %in% all.malignants)]

# ANY_SNs$sjlid %in% all.malignants
# 
# table(SMNs$sjlid %in% NMSCs$sjlid)

table(SMNs$sjlid %in% NMSCs$sjlid)

# table(NMSCs$sjlid %in% SMNs$sjlid)
# # FALSE  TRUE 
# # 14   238
# 
# table(BREASTcancer$sjlid %in% SMNs$sjlid)
# # FALSE  TRUE 
# # 14    62
# 
# table(THYROIDcancer$sjlid %in% SMNs$sjlid)
# # FALSE  TRUE 
# # 2    85
# THYROIDcancer$sjlid[!THYROIDcancer$sjlid %in% SMNs$sjlid]
# # [1] "SJL5060217" "SJL5204304"

# table(ANY_SNs$sjlid %in% all.malignants)
# # 19   586
# other <- ANY_SNs[!ANY_SNs$sjlid %in% all.malignants,]
# table(other$diag)
# # SJL5039417

other <- SMNs[!SMNs$sjlid %in% all.malignants,]
table(other$diag)




SMNs$KEY <- paste0(SMNs$sjlid, ":", SMNs$sncount)
NMSCs$KEY <- paste0(NMSCs$sjlid, ":", NMSCs$sncount)
BREASTcancer$KEY <- paste0(BREASTcancer$sjlid, ":", BREASTcancer$sncount)
THYROIDcancer$KEY <- paste0(THYROIDcancer$sjlid, ":", THYROIDcancer$sncount)
MENINGIOMA$KEY <- paste0(MENINGIOMA$sjlid, ":", MENINGIOMA$sncount)
SARCOMA$KEY <- paste0(SARCOMA$sjlid, ":", SARCOMA$sncount)

SMNs$smnTYPE <- "Unique"
SMNs$smnTYPE [SMNs$KEY %in% NMSCs$KEY] <- "NMSCs"
SMNs$smnTYPE [SMNs$KEY %in% BREASTcancer$KEY] <- "BREASTcancer"
SMNs$smnTYPE [SMNs$KEY %in% THYROIDcancer$KEY] <- "THYROIDcancer"
SMNs$smnTYPE [SMNs$KEY %in% MENINGIOMA$KEY] <- "MENINGIOMA"
SMNs$smnTYPE [SMNs$KEY %in% SARCOMA$KEY] <- "SARCOMA"
other <- SMNs[SMNs$smnTYPE == "Unique",]

table(SMNs$smnTYPE)

table(NMSCs$KEY %in% SMNs$KEY)
# FALSE  TRUE 
# 225   238
table(BREASTcancer$sjlid %in% SMNs$sjlid)
# FALSE  TRUE 
# 14    62 
table(THYROIDcancer$sjlid %in% SMNs$sjlid)
# FALSE  TRUE 
# 2    85 
table(MENINGIOMA$sjlid %in% SMNs$sjlid)
# FALSE  TRUE 
# 103    46 
table(SARCOMA$sjlid %in% SMNs$sjlid)
# FALSE  TRUE 
# 1    32

table(all.malignants %in% SMNs$sjlid)

##################################################
##################################################
##################################################
###########################
## Add demographic table ##
###########################

#################################
## 1....................SJLIFE ##
#################################

# male.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$gender == "Male",]

library(dplyr)
library(tibble)
# Fully self-contained function
create_formatted_table <- function(male_pheno_ids) {
  # Define all the ID lists inside the function
  ANY_SNs <- ANY_SNs$sjlid
  SMNs <- SMNs$sjlid
  MENINGIOMA <- MENINGIOMA$sjlid
  NMSCs <- NMSCs$sjlid
  BREASTcancer <- BREASTcancer$sjlid
  THYROIDcancer <- THYROIDcancer$sjlid
  SARCOMA <- SARCOMA$sjlid
  
  # List of conditions
  id_lists <- list(
    ANY_SNs,
    SMNs,
    MENINGIOMA,
    NMSCs,
    BREASTcancer,
    THYROIDcancer,
    SARCOMA
  )
  
  # Calculate counts for each list
  counts <- sapply(id_lists, function(ids) sum(male_pheno_ids %in% ids))
  
  # Create the formatted string
  formatted_table <- paste(counts, collapse = "/")
  
  return(formatted_table)
}

# # Example usage
# formatted_table <- create_formatted_table(male.pheno$sjlid)
# print(formatted_table)



## Diagnosis variable mapping based on SJLIFE <-> CCSS DXs mapping; Yadav's email: on 05/02/2023


library(haven)
## Load PHenotype
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)

demo <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat')
demo <- demo[demo$sjlid %in% PHENO.ANY_SN$sjlid,]
PHENO.ANY_SN$ethnic <- demo$hispanic[match(PHENO.ANY_SN$sjlid, demo$sjlid)]


wgsdiag <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/wgsdiag.sas7bdat")
head(wgsdiag)

PHENO.ANY_SN$diag <- wgsdiag$diag[match(PHENO.ANY_SN$MRN, wgsdiag$MRN)]

dim(PHENO.ANY_SN)
# View(PHENO.ANY_SN)

## Sex
table(PHENO.ANY_SN$gender)
Male = sum(PHENO.ANY_SN$gender == "Male")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$gender == "Male",]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table

Female = sum(PHENO.ANY_SN$gender != "Male")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$gender != "Male",]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table

# ## PCA ethnicity
# table(PHENO.ANY_SN$PCA.ethnicity)
# Black = sum(PHENO.ANY_SN$PCA.ethnicity == "AFR")
# White = sum(PHENO.ANY_SN$PCA.ethnicity == "EUR")
# # Asian = sum(PHENO.ANY_SN$PCA.ethnicity == "EAS")
# Other = 4401 - sum(Black, White)

## Race (Self-reported)
White = sum(PHENO.ANY_SN$race == "White")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$race == "White",]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "560/434/139/250/68/82/29"

Black = sum(PHENO.ANY_SN$race == "Black")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$race == "Black",]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "42/26/10/2/7/4/4"

# Unknown = sum(PHENO.ANY_SN$race == "Unknown"| PHENO.ANY_SN$race=="N/A")
# sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$race == "Unknown"| PHENO.ANY_SN$race=="N/A",]
# formatted_table <- create_formatted_table(sub.pheno$sjlid)
# formatted_table

Other = 4401 - sum(Black, White)
other.sjlids <- c(PHENO.ANY_SN[PHENO.ANY_SN$race == "White",]$sjlid, 
                  PHENO.ANY_SN[PHENO.ANY_SN$race == "Black",]$sjlid)
sub.pheno <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% other.sjlids,]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 3/3/0/0/1/1/0

## Ethnicity (Self-reported)
table(PHENO.ANY_SN$ethnic)
non_hispanic = sum(grepl("Non Hispanic", PHENO.ANY_SN$ethnic))

sub.pheno <- PHENO.ANY_SN[grepl("Non Hispanic", PHENO.ANY_SN$ethnic),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "595/455/148/249/72/83/33"

hispanic = sum(grepl("^Hispanic/Latino", PHENO.ANY_SN$ethnic)) # (Uknown + caribbean)
sub.pheno = PHENO.ANY_SN[grepl("^Hispanic/Latino", PHENO.ANY_SN$ethnic),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "8/7/0/2/3/4/0"


unknown = sum(grepl("Unknown", PHENO.ANY_SN$ethnic)) # (Uknown + caribbean)
sub.pheno = PHENO.ANY_SN[grepl("Unknown", PHENO.ANY_SN$ethnic),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 2/1/1/1/1/0/0

table(PHENO.ANY_SN$ethnic)
# Caribbean                              Cuban                    Mexican/Chicano 
# 5                                  2                                 27 
# Non Spanish speaking, Non Hispanic        NOS Spanish,Hispanic,Latino                       Puerto Rican 
# 4225                                 60                                 10 
# South or Central Amercian                          Southwest                            Unknown 
# 39                                  1                                 32 

as.data.frame(table(PHENO.ANY_SN$diaggrp))
# as.data.frame(table(PHENO.ANY_SN$diag))

###############
## DIAGNOSIS ##
###############

# ## Acute lymphoblastic leukemia
# ALL <- sum(grepl("lymphoblastic", PHENO.ANY_SN$diaggrp, ignore.case = T))
# AML <- sum(grepl("Acute myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
# other.leukemia <- sum(grepl("Other leukemia|Chronic myeloid leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
# 
# ## CNS
# CNS.tumors <- PHENO.ANY_SN[PHENO.ANY_SN$diaggrp == "Central nervous system (CNS)",]
# CNS.counts <- nrow(CNS.tumors)
# astrocytoma_glioma <- sum(grepl("astrocytoma|glioma", CNS.tumors$diag, ignore.case = T))
# medulloblastoma_PNET <- sum(grepl("Medulloblastoma|PNET", CNS.tumors$diag, ignore.case = T))
# ependymoma <- sum(grepl("Ependymoma", CNS.tumors$diag, ignore.case = T))
# other_CNS.tumor <- CNS.counts - (astrocytoma_glioma + medulloblastoma_PNET + ependymoma)
# 
# ## Lymphoma
# hodgkin.lymphoma <- sum(grepl("^Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# non.hodgkin.lymphoma <- sum(grepl("non-Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T)) 
# 
# ## Sarcoma
# sarcoma <- sum(grepl("sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# ewing.sarcoma <- sum(grepl("ewing", PHENO.ANY_SN$diaggrp, ignore.case = T))
# osteosarcoma <- sum(grepl("osteosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# rhabdomyosarcoma <- sum(grepl("^Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# nonrhabdomyosarcoma <- sum(grepl("Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# 
# 
# ## Embryonal
# wilms <- sum(grepl("wilms", PHENO.ANY_SN$diaggrp, ignore.case = T))
# neuroblastoma <- sum(grepl("neuroblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# germ.cell.tumor <- sum(grepl("Germ cell tumor", PHENO.ANY_SN$diaggrp, ignore.case = T))
# retinoblastoma <- sum(grepl("Retinoblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# hepatoblastoma <- sum(grepl("Liver malignancies", PHENO.ANY_SN$diaggrp, ignore.case = T))
# 
# ## Other
# melanoma <- sum(grepl("melanoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# carcinoma <- sum(grepl("carcinoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
# # Others

## Diagnosis as in CCSS (with common disease types)
leukemia <- sum(grepl("leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T))
# leukemia  <- sum(grepl("Acute lymphoblastic leukemia|Acute myeloid leukemia|Other leukemia|MDS/Acute myeloid", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("leukemia", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 233/156/112/109/11/28/5

CNS <- sum(grepl("Central nervous system", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("Central nervous system", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "52/32/23/14/1/6/6"

HD  <- sum(grepl("^Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("^Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "143/132/3/57/45/37/3"

NHL  <- sum(grepl("Non-Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("Non-Hodgkin lymphoma", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "36/31/4/22/5/2/0"

wilms  <- sum(grepl("Wilms tumor", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("Wilms tumor", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 25/22/0/12/2/2/1

bone.cancer <- sum(grepl("Osteosarcoma|Ewing sarcoma family of tumors", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("Osteosarcoma|Ewing sarcoma family of tumors", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 32/22/0/13/5/2/1

neuroblastoma  <- sum(grepl("Neuroblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("Neuroblastoma", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 20/18/1/6/1/2/2

soft.tissue.sarcoma  <- sum(grepl("Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))
sub.pheno = PHENO.ANY_SN[grepl("Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 14/11/1/3/3/2/1

Rhabdomyosarcoma  <- sum(grepl("Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))

other <- PHENO.ANY_SN[!grepl("leukemia|Central nervous system|Hodgkin|Wilms tumor|Osteosarcoma|Ewing sarcoma family of tumors|Neuroblastoma|Soft tissue sarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T),]
sub.pheno = other
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 50/39/5/16/3/6/14
##################
## Radiotherapy ##
##################
# brainRT <- sum(PHENO.ANY_SN$brainorheadrt_yn == "Y", na.rm = T)
# neckRT <- sum(PHENO.ANY_SN$neckrt_yn == "Y", na.rm = T)
# chestRT <- sum(PHENO.ANY_SN$chestrt_yn == "Y", na.rm = T)
# abdomenRT <- sum(PHENO.ANY_SN$abdomenrt_yn == "Y", na.rm = T)
# pelvisRT <- sum(PHENO.ANY_SN$pelvisrt_yn == "Y", na.rm = T)

##
brainRT <- sum(PHENO.ANY_SN$maxsegrtdose > 200, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$maxsegrtdose > 200),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 286/185/140/126/8/35/13

neckRT <- sum(PHENO.ANY_SN$maxneckrtdose > 200, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$maxneckrtdose > 200),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 238/200/41/99/45/60/6

chestRT <- sum(PHENO.ANY_SN$maxchestrtdose > 200, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$maxchestrtdose > 200),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 233/198/39/100/49/57/9

abdomenRT <- sum(PHENO.ANY_SN$maxabdrtdose > 200, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$maxabdrtdose > 200),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 214/178/39/94/36/39/7

pelvisRT <- sum(PHENO.ANY_SN$maxpelvisrtdose > 200, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$maxpelvisrtdose > 200),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# "175/137/38/71/19/35/8"
##################
## Chemotherapy ##
##################
# alkylating <- sum(PHENO.ANY_SN$aa_class_dose_5_yn == "Y", na.rm = T)
alkylating <- sum(PHENO.ANY_SN$aa_class_dose_5 > 0, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$aa_class_dose_5 > 0),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 383/297/93/158/47/57/24

# anthracyclines <- sum(PHENO.ANY_SN$anthra_jco_dose_5_yn == "Y", na.rm = T)
anthracyclines <- sum(PHENO.ANY_SN$anthra_jco_dose_5 > 0, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$anthra_jco_dose_5 > 0),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table
# 314/238/73/122/42/49/13

# epipodophyllotoxins <- sum(PHENO.ANY_SN$epitxn_dose_5_yn == "Y", na.rm = T)
epipodophyllotoxins <- sum(PHENO.ANY_SN$epitxn_dose_5 > 0, na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$epitxn_dose_5 > 0),]
formatted_table <- create_formatted_table(sub.pheno$sjlid)
formatted_table

## age at diagnosis
median.agedx <- round(median(PHENO.ANY_SN$agedx), 1)
agedx.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agedx, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agedx.IQR <- gsub(" ", "-", agedx.IQR)
agedx <- paste0(median.agedx, " (", agedx.IQR, ")")

## Age at follow up
median.age.followup <- round(median(PHENO.ANY_SN$agelstcontact), 1)
age.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
age.followup.IQR <- gsub(" ", "-", age.followup.IQR)
age.at.followup <- paste0(median.age.followup, " (", age.followup.IQR, ")")

# ## Length of follow up
# PHENO.ANY_SN$length.followup <- PHENO.ANY_SN$agelstcontact -PHENO.ANY_SN$agedx
# median.length.followup <- round(median(PHENO.ANY_SN$length.followup), 1)
# length.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$length.followup, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
# length.followup.IQR <- gsub(" ", "-", length.followup.IQR)
# lenght.followup <- paste0(median.length.followup, " (", length.followup.IQR, ")")


## rbind all
# df <- as.data.frame(t(cbind.data.frame(Male, Female, White, Black, Other, ALL, AML, other.leukemia, astrocytoma_glioma, medulloblastoma_PNET, ependymoma, other_CNS.tumor, hodgkin.lymphoma, non.hodgkin.lymphoma, ewing.sarcoma, osteosarcoma, rhabdomyosarcoma, nonrhabdomyosarcoma, wilms, neuroblastoma, germ.cell.tumor, retinoblastoma, hepatoblastoma, melanoma, carcinoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup, lenght.followup)))


# df.sjlife <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup, lenght.followup)))
df.sjlife <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, Rhabdomyosarcoma, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup)))
df.sjlife$percent <- c(round((as.numeric(df.sjlife$V1[!grepl("agedx|age.at", rownames(df.sjlife))])/4401)*100, 1), NA, NA)





###############################
## 2....................CCSS ##
###############################

## **
# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)


common.cols <- colnames(PHENO.ANY_SN)[colnames(PHENO.ANY_SN) %in% colnames(subneo)]
common.cols <- common.cols[-1]
subneo2 <- subneo[,!colnames(subneo) %in% common.cols]
merged_data <- merge(PHENO.ANY_SN, subneo2, by = c("ccssid"), all.x = TRUE)
merged_data <- merged_data[,!grepl("AGE_AT_LAST_CONTACT", colnames(merged_data))]
## Send this to QI
# saveRDS(merged_data, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_complete_data.rds")

subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]

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
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
ANY_SNs <- ANY_SNs[!ANY_SNs$ccssid %in% subneo.within5$ccssid ,]

##########
## SMNs ##
##########

# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

#########################
## Keep malignant only ##
#########################
subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

#############
## Any SNs ##
#############
table(subneo$seersmn)
# No  Yes 
# 2206  895 

SMNs <- subneo[grepl("Yes", subneo$seersmn),]

# SN within 5 years
subneo.within5 <- SMNs[SMNs$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

SMNs <- setDT(SMNs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
SMNs <- SMNs[!SMNs$ccssid %in% subneo.within5$ccssid ,]
nrow(SMNs)

########### 
## NMSCs ##
###########
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)


## Read NMSC data from Qi
data1 = read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/sns2022.sas7bdat")
data1=as.data.frame(data1)
# data1$ccssid <- paste0(data1$ccssid, "_", data1$ccssid)
data1$KEY <- paste0(data1$ccssid,":",data1$d_candx)


#########################
## Subsequent neoplasm ##
#########################
## ADD NMSC from Qi
subneo$d_candx <- as.Date(subneo$d_candx, format = "%d%b%Y")
subneo$KEY <- paste0(subneo$ccssid,":",subneo$d_candx)
table(subneo$KEY %in% data1$KEY)
# FALSE  TRUE 
# 6307  3440 
table(data1$KEY %in% subneo$KEY)
# FALSE  TRUE 
# 4629  4434 
subneo$nmsc <- data1$nmsc[match(subneo$KEY, data1$KEY)]
subneo$candxo3 <- data1$candxo3[match(subneo$KEY, data1$KEY)]
# cc <- cbind.data.frame(subneo$KEY, subneo$nmsc, subneo$AGE.ANY_SN, subneo$groupdx3)

# Now get age of SN after first cancer
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

#########################
## Keep malignant only ##
#########################
subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]

########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)

# This will include basal cell and squamous cell
# NMSCs <- subneo[which((subneo$nmsc ==1| (subneo$nmsc == 2 & subneo$groupdx3 == "Skin"))),]
NMSCs <- subneo[which(subneo$nmsc ==1),]

# cc <- cbind.data.frame(NMSC$KEY, NMSC$nmsc, NMSC$AGE.ANY_SN, NMSC$groupdx3)

subneo.within5 <- NMSCs[NMSCs$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

# NMSCs <- subneo[grepl("skin", subneo$groupdx3, ignore.case = T),]
NMSCs <- setDT(NMSCs)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(NMSCs)
# 729



# NMSCs.saved <- NMSCs
# NMSCs <- NMSCs[,-c(49,52:53)]
# save(NMSCs, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/NMSCs.RData")


## Based on Qi's email on 9/26/2024
# # Filter for NMSC cases
sngroups <- NMSCs
nmsc <- sngroups %>%
  filter(nmsc == 1)
# Frequency table for candxo3 in the nmsc dataset
table(nmsc$candxo3)
# Yadav said: In the final analysis, Basosquamous carcinoma was excluded, per Smita. (but we will include)
# Filter for Basosquamous carcinoma
basosquamous <- sngroups %>%
  filter(candxo3 == 8094.3)  # All 8094.3 are NMSC
# Frequency table for nmsc in the basosquamous dataset
table(basosquamous$nmsc)
# Create sngroup data frame with conditions for BCC and BCC_exclude
sngroup <- sngroups %>%
  mutate(BCC = ifelse(candxo3 > 8090 & candxo3 != 8094.3 & nmsc == 1, 1, 0),
         BCC_exclude = ifelse(candxo3 == 8094.3, 1, 0))
# Optionally, if you want to see the resulting dataset:
head(sngroup)


# exclude basosquamous
table(sngroup$BCC_exclude != 1)
# FALSE  TRUE 
# 1   728 
# sngroup <- sngroup[sngroup$BCC_exclude != 1]
NMSCs <- NMSCs[NMSCs$ccssid %in% sngroup$ccssid,]
NMSCs$nmsc_type <- ifelse(sngroup$BCC==1, "BCC", "SCC")
dim(NMSCs)
# 729
NMSCs$nmsc_type[NMSCs$ccssid %in% sngroup$ccssid[sngroup$BCC_exclude == 1]] <- "Basosquamous"
NMSCs <- NMSCs[!NMSCs$ccssid %in% subneo.within5$ccssid,]

###################
## Breast cancer ##
###################

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

#########################
## Keep malignant only ##
#########################
subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]
########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
BREASTcancer <- subneo[grepl("breast", subneo$groupdx3, ignore.case = T),]

## SN within 5 years
subneo.within5 <- BREASTcancer[BREASTcancer$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

BREASTcancer <- setDT(BREASTcancer)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(BREASTcancer)
# 295

# "a_dx"  : Primary cancer diagnosis age
# "a_end" : age at last contact
# "d_candx" : Date when second cancer was diagnosed                                    
# "a_candx": Age at second cancer diagnosis


# dat[,c("ccssid","strokedt","event","dob","agelstcontact","agedx")]
BREASTcancer$gradeage <- BREASTcancer$gradedt
BREASTcancer$gradedt <- as.Date(BREASTcancer$d_candx, format = "%d%b%Y")
## Calculate DOB
BREASTcancer$dob <- BREASTcancer$gradedt - as.numeric(BREASTcancer$gradeage) * 365.2422


BREASTcancer <- BREASTcancer[BREASTcancer$SEX == "Female",]

####################
## Thyroid cancer ##
####################

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

#########################
## Keep malignant only ##
#########################
subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]
########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
THYROIDcancer <- subneo[grepl("thyroid", subneo$groupdx3, ignore.case = T),]

## SN within 5 years
subneo.within5 <- THYROIDcancer[THYROIDcancer$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

THYROIDcancer <- setDT(THYROIDcancer)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(THYROIDcancer)
# 167

# "a_dx"  : Primary cancer diagnosis age
# "a_end" : age at last contact
# "d_candx" : Date when second cancer was diagnosed                                    
# "a_candx": Age at second cancer diagnosis


# dat[,c("ccssid","strokedt","event","dob","agelstcontact","agedx")]
THYROIDcancer$gradeage <- THYROIDcancer$gradedt
THYROIDcancer$gradedt <- as.Date(THYROIDcancer$d_candx, format = "%d%b%Y")
## Calculate DOB
THYROIDcancer$dob <- THYROIDcancer$gradedt - as.numeric(THYROIDcancer$gradeage) * 365.2422

# Removing samples with SN within the 5 years of childhood cancer **
THYROIDcancer <- THYROIDcancer[!(THYROIDcancer$ccssid %in% subneo.within5$ccssid),]
# 4

################
## Meningioma ##
################

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

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


subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]

#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
MENINGIOMA <- subneo[grepl("meningioma", subneo$groupdx3, ignore.case = T),]

## SN within 5 years
subneo.within5 <- MENINGIOMA[MENINGIOMA$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 0

MENINGIOMA <- setDT(MENINGIOMA)[,.SD[which.min(gradedt)],by=ccssid][order(gradedt, decreasing = FALSE)]
nrow(MENINGIOMA)
# 256

# "a_dx"  : Primary cancer diagnosis age
# "a_end" : age at last contact
# "d_candx" : Date when second cancer was diagnosed                                    
# "a_candx": Age at second cancer diagnosis


# dat[,c("ccssid","strokedt","event","dob","agelstcontact","agedx")]
MENINGIOMA$gradeage <- MENINGIOMA$gradedt
MENINGIOMA$gradedt <- as.Date(MENINGIOMA$d_candx, format = "%d%b%Y")
## Calculate DOB
MENINGIOMA$dob <- MENINGIOMA$gradedt - as.numeric(MENINGIOMA$gradeage) * 365.2422


# Removing samples with SN within the 5 years of childhood cancer **
sum(PHENO.ANY_SN$ccssid %in% subneo.within5$ccssid)
# 0
MENINGIOMA <- MENINGIOMA[!MENINGIOMA$ccssid %in% subneo.within5$ccssid,]
dim(MENINGIOMA)
# 256

#############
## Sarcoma ##
#############

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_Genetic_data_P_LP_v17.Rdata")

## Edit lifestyle variables
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

#########################
## Subsequent neoplasm ##
#########################
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


#########################
## Keep malignant only ##
#########################
subneo$malKey <- paste(subneo$ccssid, subneo$groupdx3, subneo$a_candx, subneo$count, sep = ":")
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus <- malignantStatus[malignantStatus$a_candx !=".",]
malignantStatus$malKey <- paste(malignantStatus$ccssid, malignantStatus$groupdx3, malignantStatus$a_candx, malignantStatus$count, sep = ":")
# malignantStatus$dupli <- duplicated(malignantStatus$Key)

## Add malignant status
subneo$seersmn <- malignantStatus$seersmn[match(subneo$malKey, malignantStatus$malKey)]
########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$ccssid))
# 1619

#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN.
# For this, I will first sort the table by date
library(data.table)
SARCOMA <- subneo[grepl("sarcoma", subneo$groupdx3, ignore.case = T),]

## SN within 5 years
subneo.within5 <- SARCOMA[SARCOMA$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))

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


SARCOMA <- SARCOMA[!SARCOMA$ccssid %in% to.remove,]

# "a_dx"  : Primary cancer diagnosis age
# "a_end" : age at last contact
# "d_candx" : Date when second cancer was diagnosed                                    
# "a_candx": Age at second cancer diagnosis


# dat[,c("ccssid","strokedt","event","dob","agelstcontact","agedx")]
SARCOMA$gradeage <- SARCOMA$gradedt
SARCOMA$gradedt <- as.Date(SARCOMA$d_candx, format = "%d%b%Y")
## Calculate DOB
SARCOMA$dob <- SARCOMA$gradedt - as.numeric(SARCOMA$gradeage) * 365.2422

# Removing samples with SN within the 5 years of childhood cancer **
sum(SARCOMA$ccssid %in% subneo.within5$ccssid)

SARCOMA <- SARCOMA[!(SARCOMA$ccssid %in% subneo.within5$ccssid),]

#######################
## Demographic table ##
#######################

# male.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$gender == "Male",]

library(dplyr)
library(tibble)
# Fully self-contained function
create_formatted_table <- function(male_pheno_ids) {
  # Define all the ID lists inside the function
  ANY_SNs <- ANY_SNs$ccssid
  SMNs <- SMNs$ccssid
  MENINGIOMA <- MENINGIOMA$ccssid
  NMSCs <- NMSCs$ccssid
  BREASTcancer <- BREASTcancer$ccssid
  THYROIDcancer <- THYROIDcancer$ccssid
  SARCOMA <- SARCOMA$ccssid
  
  # List of conditions
  id_lists <- list(
    ANY_SNs,
    SMNs,
    MENINGIOMA,
    NMSCs,
    BREASTcancer,
    THYROIDcancer,
    SARCOMA
  )
  
  # Calculate counts for each list
  counts <- sapply(id_lists, function(ids) sum(male_pheno_ids %in% ids))
  
  # Create the formatted string
  formatted_table <- paste(counts, collapse = "/")
  
  return(formatted_table)
}

# # Example usage
# formatted_table <- create_formatted_table(male.pheno$ccssid)
# print(formatted_table)

## Load PHenotype
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

dim(PHENO.ANY_SN)
# View(PHENO.ANY_SN)

## Sex
table(PHENO.ANY_SN$gender)
Male = sum(PHENO.ANY_SN$gender == "Male")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$gender == "Male",]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# "645/264/107/328/0/56/29"

Female = sum(PHENO.ANY_SN$gender != "Male")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$gender != "Male",]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 966/498/149/400/290/107/32

## Add race and Ethnicity
ccss.new <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/1_22_2025/add2var_race.txt", sep = "\t", header = T)
PHENO.ANY_SN$race <- ccss.new$racegroup[match(PHENO.ANY_SN$ccssid, ccss.new$ccssid)]
PHENO.ANY_SN$ethnic <- ccss.new$hispgroup[match(PHENO.ANY_SN$ccssid, ccss.new$ccssid)]

## Race (Self-reported)
White = sum(PHENO.ANY_SN$race == "White")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$race == "White",]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 1544/724/243/709/278/155/59

Black = sum(PHENO.ANY_SN$race == "Black")
sub.pheno <- PHENO.ANY_SN[PHENO.ANY_SN$race == "Black",]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 19/13/3/2/5/2/2




Other = 7943 - sum(Black, White)
other.ccssids <- c(PHENO.ANY_SN[PHENO.ANY_SN$race == "White",]$ccssid, 
                   PHENO.ANY_SN[PHENO.ANY_SN$race == "Black",]$ccssid)
sub.pheno <- PHENO.ANY_SN[!PHENO.ANY_SN$ccssid %in% other.ccssids,]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# "48/25/10/17/7/6/0"

## Ethnicity (Self-reported)
table(PHENO.ANY_SN$ethnic)
non_hispanic = sum(grepl("No", PHENO.ANY_SN$ethnic))

sub.pheno <- PHENO.ANY_SN[grepl("No", PHENO.ANY_SN$ethnic),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# "1490/706/229/672/272/144/57"

hispanic = sum(grepl("Yes", PHENO.ANY_SN$ethnic)) # (Uknown + caribbean)
sub.pheno = PHENO.ANY_SN[grepl("Yes", PHENO.ANY_SN$ethnic),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# "66/30/16/26/10/8/2"


unknown = sum(grepl("Unknown", PHENO.ANY_SN$ethnic)) # (Uknown + caribbean)
sub.pheno = PHENO.ANY_SN[grepl("Unknown", PHENO.ANY_SN$ethnic),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 55/26/11/30/8/11/2

as.data.frame(table(PHENO.ANY_SN$diagnose))
# as.data.frame(table(PHENO.ANY_SN$diag))


###############
## DIAGNOSIS ##
###############
leukemia <- sum(PHENO.ANY_SN$diagnose == "Leukemia", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "Leukemia"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 389/120/122/211/25/33/7

CNS <- sum(PHENO.ANY_SN$diagnose == "CNS", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "CNS"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 268/95/102/94/6/29/5

HD <- sum(PHENO.ANY_SN$diagnose == "HD", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "HD"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 420/235/6/244/153/47/14

NHL <- sum(PHENO.ANY_SN$diagnose == "NHL", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "NHL"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 100/58/9/39/14/8/2

wilms <- sum(PHENO.ANY_SN$diagnose == "Kidney (Wilms)", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "Kidney (Wilms)"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 83/46/1/33/15/8/4

bone.cancer <- sum(PHENO.ANY_SN$diagnose == "Bone cancer", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "Bone cancer"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 156/98/3/47/49/21/10
neuroblastoma <- sum(PHENO.ANY_SN$diagnose == "Neuroblastoma", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "Neuroblastoma"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 65/39/3/15/5/8/4

soft.tissue.sarcoma <- sum(PHENO.ANY_SN$diagnose == "Soft tissue sarcoma", na.rm = T)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$diagnose == "Soft tissue sarcoma"),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 130/71/10/45/23/9/15
# Rhabdomyosarcoma  <- sum(grepl("Rhabdomyosarcoma", PHENO.ANY_SN$diaggrp, ignore.case = T))

##################
## Radiotherapy ##
##################
PHENO.ANY_SN$maxsegrtdose <- as.numeric(PHENO.ANY_SN$maxsegrtdose)
PHENO.ANY_SN$neckmaxrtdose <- as.numeric(PHENO.ANY_SN$neckmaxrtdose)
PHENO.ANY_SN$chestmaxrtdose <- as.numeric(PHENO.ANY_SN$chestmaxrtdose)
PHENO.ANY_SN$abdmaxrtdose <- as.numeric(PHENO.ANY_SN$abdmaxrtdose)
PHENO.ANY_SN$pelvismaxrtdose <- as.numeric(PHENO.ANY_SN$pelvismaxrtdose)

brainRT <- sum(PHENO.ANY_SN$maxsegrtdose > 0, na.rm = T)
round((brainRT/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$maxsegrtdose > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 587/182/222/282/29/53/10

neckRT <- sum(PHENO.ANY_SN$neckmaxrtdose > 0, na.rm = T)
round((neckRT/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$neckmaxrtdose > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 631/314/88/347/154/88/19

chestRT <- sum(PHENO.ANY_SN$chestmaxrtdose > 0, na.rm = T)
round((chestRT/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$chestmaxrtdose > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 675/349/88/357/175/95/21

abdomenRT <- sum(PHENO.ANY_SN$abdmaxrtdose > 0, na.rm = T)
round((abdomenRT/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$abdmaxrtdose > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 559/285/83/284/115/68/27

pelvisRT <- sum(PHENO.ANY_SN$pelvismaxrtdose > 0, na.rm = T)
round((pelvisRT/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$pelvismaxrtdose > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 448/220/80/217/64/59/23
##################
## Chemotherapy ##
##################
alkylating <- sum(PHENO.ANY_SN$alk_CED5 > 0, na.rm = T)
round((alkylating/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$alk_CED5 > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 697/356/89/293/124/82/35

anthracyclines <- sum(PHENO.ANY_SN$anth_DED5 > 0, na.rm = T)
round((anthracyclines/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$anth_DED5 > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table
# 506/264/59/199/107/64/21

epipodophyllotoxins <- sum(PHENO.ANY_SN$epipdose5 > 0, na.rm = T)
round((epipodophyllotoxins/7943)*100,1)
sub.pheno = PHENO.ANY_SN[which(PHENO.ANY_SN$epipdose5 > 0),]
formatted_table <- create_formatted_table(sub.pheno$ccssid)
formatted_table

## age at diagnosis
median.agedx <- round(median(PHENO.ANY_SN$agedx), 1)
agedx.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agedx, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agedx.IQR <- gsub(" ", "-", agedx.IQR)
agedx <- paste0(median.agedx, " (", agedx.IQR, ")")

# ## Age at follow up
# median.age.followup <- round(median(PHENO.ANY_SN$agelstcontact), 1)
# age.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
# age.followup.IQR <- gsub(" ", "-", age.followup.IQR)
# age.at.followup <- paste0(median.age.followup, " (", age.followup.IQR, ")")

## Length of follow up
# PHENO.ANY_SN$length.followup <- PHENO.ANY_SN$agelstcontact -PHENO.ANY_SN$agedx
# median.length.followup <- round(median(PHENO.ANY_SN$length.followup), 1)
# length.followup.IQR <- paste0(unname(round((quantile(PHENO.ANY_SN$length.followup, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
# length.followup.IQR <- gsub(" ", "-", length.followup.IQR)
# lenght.followup <- paste0(median.length.followup, " (", length.followup.IQR, ")")


## Common primary dx in CCSS and SJLIFE
# CCSS: SJLIFE
# Bone cancer: Osteosarcoma|Ewing sarcoma family of tumors|
#   CNS: Central nervous system (CNS)
# Leukemia: Acute lymphoblastic leukemia|Acute myeloid leukemia|Other leukemia|MDS/Acute myeloid leukemia|Chronic myeloid leukemia
# Neuroblastoma: Neuroblastoma
# Kidney (Wilms): Wilms tumor
# HD: Hodgkin lymphoma
# NHL: Non-Hodgkin lymphoma
# Soft tissue sarcoma: Soft tissue sarcoma


# df.ccss <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup, lenght.followup)))
df.ccss <- as.data.frame(t(cbind.data.frame(Male, Female, bone.cancer, Rhabdomyosarcoma, CNS, HD, wilms, leukemia, neuroblastoma, NHL, soft.tissue.sarcoma, brainRT, neckRT, chestRT, abdomenRT, pelvisRT, alkylating, anthracyclines, epipodophyllotoxins, agedx, age.at.followup)))
df.ccss$percent <- c(round((as.numeric(df.ccss$V1[!grepl("agedx|age.at", rownames(df.ccss))])/7943)*100, 1), NA, NA)

colnames(df.sjlife) <- c("SJLIFE (n)", "sjlife%")
colnames(df.ccss) <- c("CCSS (n)", "CCSS%")

df <- cbind(df.sjlife, df.ccss)

#########################
## Age at last contact ##
#########################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/")
ccss <- readRDS("CCSS_complete_data.rds")

agelstcontact <- round(median(ccss$agelstcontact, na.rm= T),1)
agelstcontact.IQR <- paste0(unname(round((quantile(ccss$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agelstcontact.IQR <- gsub(" ", "-", agelstcontact.IQR)
agelstcontact <- paste0(agelstcontact, " (", agelstcontact.IQR, ")")
# "38 (30.9-46.3)"

## Attained age
ccss$a_candx <- as.numeric(ccss$a_candx)
ccss$agelstcontact[!is.na(ccss$a_candx)] <- ccss$a_candx[!is.na(ccss$a_candx)] 
agelstcontact <- round(median(ccss$agelstcontact, na.rm= T),1)
agelstcontact.IQR <- paste0(unname(round((quantile(ccss$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agelstcontact.IQR <- gsub(" ", "-", agelstcontact.IQR)
agelstcontact <- paste0(agelstcontact, " (", agelstcontact.IQR, ")")
# "36.1 (29.7-43.5)"

## SJLIFE
## Age at last contact
sj <- readRDS("SJLIFE_complete_data.rds")
agelstcontact <- round(median(sj$agelstcontact, na.rm= T),1)
agelstcontact.IQR <- paste0(unname(round((quantile(sj$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agelstcontact.IQR <- gsub(" ", "-", agelstcontact.IQR)
agelstcontact <- paste0(agelstcontact, " (", agelstcontact.IQR, ")")
agelstcontact
# 36.2 (26.2-46.9)

## Attained age
sj$dob <- as.Date(sj$dob)
sj$diagdt <- as.Date(sj$diagdt)

# Calculate age at diagnosis
sj$a_candx <- as.numeric(difftime(sj$gradedt, sj$dob, units = "days")) / 365.25

sj$agelstcontact[!is.na(sj$a_candx)] <- sj$a_candx[!is.na(sj$a_candx)] 
agelstcontact <- round(median(sj$agelstcontact, na.rm= T),1)
agelstcontact.IQR <- paste0(unname(round((quantile(sj$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
agelstcontact.IQR <- gsub(" ", "-", agelstcontact.IQR)
agelstcontact <- paste0(agelstcontact, " (", agelstcontact.IQR, ")")



















######
SMNs$KEY <- paste0(SMNs$ccssid, ":", SMNs$count)
NMSCs$KEY <- paste0(NMSCs$ccssid, ":", NMSCs$count)
BREASTcancer$KEY <- paste0(BREASTcancer$ccssid, ":", BREASTcancer$count)
THYROIDcancer$KEY <- paste0(THYROIDcancer$ccssid, ":", THYROIDcancer$count)
MENINGIOMA$KEY <- paste0(MENINGIOMA$ccssid, ":", MENINGIOMA$count)
SARCOMA$KEY <- paste0(SARCOMA$ccssid, ":", SARCOMA$count)

SMNs$smnTYPE <- "Unique"
SMNs$smnTYPE [SMNs$KEY %in% NMSCs$KEY] <- "NMSCs"
SMNs$smnTYPE [SMNs$KEY %in% BREASTcancer$KEY] <- "BREASTcancer"
SMNs$smnTYPE [SMNs$KEY %in% THYROIDcancer$KEY] <- "THYROIDcancer"
SMNs$smnTYPE [SMNs$KEY %in% MENINGIOMA$KEY] <- "MENINGIOMA"
SMNs$smnTYPE [SMNs$KEY %in% SARCOMA$KEY] <- "SARCOMA"
other <- SMNs[SMNs$smnTYPE == "Unique",]
######

table(SMNs$smnTYPE)
table(other$groupdx3)



