#################################
## 1....................SJLIFE ##
#################################
#########################
## Load Phenotype data ##
#########################
library(haven)
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


## Create SN and SMN tables


SMNs$KEY <- paste0(SMNs$sjlid, ":", SMNs$sncount)
NMSCs$KEY <- paste0(NMSCs$sjlid, ":", NMSCs$sncount)
BREASTcancer$KEY <- paste0(BREASTcancer$sjlid, ":", BREASTcancer$sncount)
THYROIDcancer$KEY <- paste0(THYROIDcancer$sjlid, ":", THYROIDcancer$sncount)
MENINGIOMA$KEY <- paste0(MENINGIOMA$sjlid, ":", MENINGIOMA$sncount)
SARCOMA$KEY <- paste0(SARCOMA$sjlid, ":", SARCOMA$sncount)

ANY_SNs.sjlife <- ANY_SNs
SMNs.sjlife <- SMNs

##########
## SMNs ##
##########
SMNs.sjlife$new.diagnosis <- NA

SMNs.sjlife$new.diagnosis[grepl("basal cell|squamous cell|squamous", SMNs.sjlife$diag, ignore.case = T)] <- "NMSCs"
SMNs.sjlife$new.diagnosis[grepl("breast", SMNs.sjlife$diaggrp, ignore.case = T)] <- "Breast"
SMNs.sjlife$new.diagnosis[grepl("thyroid", SMNs.sjlife$diag, ignore.case = T)] <- "Thyroid"
SMNs.sjlife$new.diagnosis[grepl("meningioma", SMNs.sjlife$diag, ignore.case = T)] <- "Meningioma"
SMNs.sjlife$new.diagnosis[grepl("Sarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Sarcoma"
SMNs.sjlife$new.diagnosis[is.na(SMNs.sjlife$new.diagnosis)] <- "Other"

other.smn.sjlife <- SMNs.sjlife[SMNs.sjlife$new.diagnosis=="Other",]
malignant.sjlife <- SMNs.sjlife[SMNs.sjlife$new.diagnosis!="Other",]

SMNs.sjlife$new.diagnosis.cleaned <- NA
SMNs.sjlife$new.diagnosis.cleaned[grepl("basal cell", SMNs.sjlife$diag, ignore.case = T)] <- "Basal cell carcinoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("squamous", SMNs.sjlife$diag, ignore.case = T)] <- "Squamous cell carcinoma"

SMNs.sjlife$new.diagnosis.cleaned[grepl("ductal|Mammary Carcinoma", SMNs.sjlife$diag, ignore.case = T)] <- "Infiltrating ductal carcinoma of breast"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Secretory", SMNs.sjlife$diag, ignore.case = T)] <- "Secretory carcinoma of breast"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Metaplastic, Breast", SMNs.sjlife$diag, ignore.case = T)] <- "Metaplastic carcinoma of breast"
cc <- SMNs.sjlife[SMNs.sjlife$new.diagnosis == "Breast",]


SMNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Papillary|Carcinoma, Hurthle Cell and Papillary|Papillary carcinoma", SMNs.sjlife$diag, ignore.case = T)] <- "Carcinoma, Papillary, Thyroid"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Adenocarcinoma, Papillary", SMNs.sjlife$diag, ignore.case = T)] <- "Adenocarcinoma, Papillary, Thyroid"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Follicular, Thyroid", SMNs.sjlife$diag, ignore.case = T)] <- "Adenocarcinoma, Follicular, Thyroid"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Microcarcinoma, Papillary", SMNs.sjlife$diag, ignore.case = T)] <- "Microcarcinoma, Papillary, Thyroid"

SMNs.sjlife$new.diagnosis.cleaned[grepl("meningioma", SMNs.sjlife$diag, ignore.case = T)] <- "Meningioma"

SMNs.sjlife$new.diagnosis.cleaned[grepl("Osteosarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Osteosarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Chondrosarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Chondrosarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Dermatofibrosarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Dermatofibrosarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Leiomyosarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Leiomyosarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Ewing", SMNs.sjlife$diag, ignore.case = T)] <- "Ewing's Sarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Histiocytic Sarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Histiocytic Sarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Epithelioid", SMNs.sjlife$diag, ignore.case = T)] <- "Epithelioid Sarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Liposarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Liposarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Sarcoma, High Grade Spindle Cell", SMNs.sjlife$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Sarcoma, High Grade, Adrenal Gland", SMNs.sjlife$diag, ignore.case = T)] <- "Sarcoma, High Grade Adrenal Gland"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Sarcoma, Soft Tissue", SMNs.sjlife$diag, ignore.case = T)] <- "Leiomyosarcoma"
SMNs.sjlife$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", SMNs.sjlife$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"

# SMNs.sjlife$diag[grepl("sarcoma" , SMNs.sjlife$diag, ignore.case = T)][!SMNs.sjlife$diag[grepl("sarcoma" , SMNs.sjlife$diag, ignore.case = T)] %in% SMNs.sjlife$diag[grepl("sarcoma" , SMNs.sjlife$new.diagnosis.cleaned, ignore.case = T)]]


SMNs.sjlife$new.diagnosis.cleaned[grepl("Sarcoma", SMNs.sjlife$diaggrp, ignore.case = T)& grepl("^Malignant|Dermatofibro|Ganglioma|papillary|firbromatosis|fibromyxoid|Giant cell|Neurilemoma|lipoma|Haemangioblastoma|comedocarcinoma", SMNs.sjlife$diag, ignore.case = T)] <- "Soft tissue tumor"
SMNs.sjlife$new.diagnosis.cleaned[is.na(SMNs.sjlife$new.diagnosis.cleaned)] <- "Other"

# Malignant Peripheral Nerve Sheath Tumor, Neck, Left

other.SMN.sjlife <- SMNs.sjlife[(SMNs.sjlife$new.diagnosis.cleaned=="Other"),]
table(ifelse(grepl("^Malignant", other.SMN.sjlife$icdo3behavior), "M", "B"))
# M 
# 96
SMNs.sjlife$new.diagnosis.cleaned2 <- SMNs.sjlife$new.diagnosis.cleaned
SMNs.sjlife$new.diagnosis.cleaned2[grepl("Soft tissue tumor", SMNs.sjlife$new.diagnosis.cleaned, ignore.case = T)] <- "Other"
other.SMN.sjlife <- SMNs.sjlife[(SMNs.sjlife$new.diagnosis.cleaned2=="Other"),]

cc <- SMNs.sjlife[grepl("Sarcoma, Soft Tissue", SMNs.sjlife$new.diagnosis.cleaned),]

View(table(SMNs.sjlife$new.diagnosis.cleaned2))

#############
## Any SNs ##
#############
ANY_SNs.sjlife$KEY <- paste0(ANY_SNs.sjlife$sjlid, ":", ANY_SNs.sjlife$sncount)

ANY_SNs.sjlife$new.diagnosis <- NA

ANY_SNs.sjlife$new.diagnosis[grepl("basal cell|squamous cell|squamous", ANY_SNs.sjlife$diag, ignore.case = T)] <- "NMSCs"
ANY_SNs.sjlife$new.diagnosis[grepl("breast", ANY_SNs.sjlife$diaggrp, ignore.case = T)] <- "Breast"
ANY_SNs.sjlife$new.diagnosis[grepl("thyroid", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Thyroid"
ANY_SNs.sjlife$new.diagnosis[grepl("meningioma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Meningioma"
ANY_SNs.sjlife$new.diagnosis[grepl("Sarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Sarcoma"
ANY_SNs.sjlife$new.diagnosis[is.na(ANY_SNs.sjlife$new.diagnosis)] <- "Other"

other.sjlife <- ANY_SNs.sjlife[ANY_SNs.sjlife$new.diagnosis=="Other",]
checkanySN.sjlife <- ANY_SNs.sjlife[ANY_SNs.sjlife$new.diagnosis!="Other",]

# Breast               Meningioma Non Melanoma Skin Cancer       Other Solid Tumors 
# 54                      125                      184                       22 
# Sarcoma                  Thyroid 
# 26                       68 

ANY_SNs.sjlife$new.diagnosis.cleaned <- NA
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("basal cell", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Basal cell carcinoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("squamous", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Squamous cell carcinoma"

ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("ductal|Mammary Carcinoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Infiltrating ductal carcinoma of breast"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Secretory", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Secretory carcinoma of breast"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Metaplastic, Breast", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Metaplastic carcinoma of breast"
cc <- ANY_SNs.sjlife[ANY_SNs.sjlife$new.diagnosis == "Breast",]


ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Papillary|Carcinoma, Hurthle Cell and Papillary|Papillary carcinoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Carcinoma, Papillary, Thyroid"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Adenocarcinoma, Papillary", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Adenocarcinoma, Papillary, Thyroid"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Carcinoma, Follicular, Thyroid", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Adenocarcinoma, Follicular, Thyroid"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Microcarcinoma, Papillary", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Microcarcinoma, Papillary, Thyroid"

ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("meningioma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Meningioma"

ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Osteosarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Osteosarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Chondrosarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Chondrosarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Dermatofibrosarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Sarcoma, Soft Tissue"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Ewing", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Ewing's Sarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Histiocytic Sarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Histiocytic Sarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Epithelioid", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Epithelioid Sarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Myxoid Liposarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Liposarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Sarcoma, High Grade Spindle Cell", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Sarcoma, High Grade, Adrenal Gland", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Sarcoma, High Grade Adrenal Gland"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Sarcoma, Soft Tissue", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Sarcoma, Soft Tissue"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
ANY_SNs.sjlife$new.diagnosis.cleaned[grepl("Leiomyosarcoma", ANY_SNs.sjlife$diag, ignore.case = T)] <- "Leiomyosarcoma"

ANY_SNs.sjlife$new.diagnosis.cleaned2 <- ANY_SNs.sjlife$new.diagnosis.cleaned
ANY_SNs.sjlife$new.diagnosis.cleaned2[grepl("Soft tissue|Dermatofibroma", ANY_SNs.sjlife$new.diagnosis.cleaned2, ignore.case = T)] <- "Other"
ANY_SNs.sjlife$new.diagnosis.cleaned2[is.na(ANY_SNs.sjlife$new.diagnosis.cleaned)] <- "Other"

ANY_SNs.sjlife$new.diagnosis.cleaned.malignant <- ANY_SNs.sjlife$new.diagnosis.cleaned2
ANY_SNs.sjlife$new.diagnosis.cleaned.malignant <- paste0(ANY_SNs.sjlife$new.diagnosis.cleaned.malignant, "_", ifelse(grepl("^Malignant", ANY_SNs.sjlife$icdo3behavior), "M", "B"))
ANY_SNs.sjlife$new.diagnosis.cleaned.malignant[grepl("NA_M|NA_B",ANY_SNs.sjlife$new.diagnosis.cleaned.malignant)] <- NA


other.anySN.sjlife <- ANY_SNs.sjlife[ANY_SNs.sjlife$new.diagnosis.cleaned2 =="Other",]

View(table(ANY_SNs.sjlife$new.diagnosis.cleaned.malignant))
# other.anySN.sjlife <- ANY_SNs.sjlife[is.na(ANY_SNs.sjlife$new.diagnosis.cleaned),]

table(ifelse(grepl("^Malignant", other.anySN.sjlife$icdo3behavior), "M", "B"))
# B  M 
# 25 98

other.SMN.sjlife[!other.SMN.sjlife$KEY %in% other.sjlife.sn$KEY,]


# cc <- other.SMN.ccss[other.SMN.ccss$ccssid %in% NMSCs$ccssid,]

#############
# SMNs$smnTYPE <- "Unique"
# SMNs$smnTYPE [SMNs$KEY %in% NMSCs$KEY] <- "NMSCs"
# SMNs$smnTYPE [SMNs$KEY %in% BREASTcancer$KEY] <- "BREASTcancer"
# SMNs$smnTYPE [SMNs$KEY %in% THYROIDcancer$KEY] <- "THYROIDcancer"
# SMNs$smnTYPE [SMNs$KEY %in% MENINGIOMA$KEY] <- "MENINGIOMA"
# SMNs$smnTYPE [SMNs$KEY %in% SARCOMA$KEY] <- "SARCOMA"
# other <- SMNs[SMNs$smnTYPE == "Unique",]
# 
# table(SMNs$smnTYPE)
# 
# table(NMSCs$KEY %in% SMNs$KEY)
# # FALSE  TRUE 
# # 225   238
# table(BREASTcancer$sjlid %in% SMNs$sjlid)
# # FALSE  TRUE 
# # 14    62 
# table(THYROIDcancer$sjlid %in% SMNs$sjlid)
# # FALSE  TRUE 
# # 2    85 
# table(MENINGIOMA$sjlid %in% SMNs$sjlid)
# # FALSE  TRUE 
# # 103    46 
# table(SARCOMA$sjlid %in% SMNs$sjlid)
# # FALSE  TRUE 
# # 1    32
# 
# table(all.malignants %in% SMNs$sjlid)

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

###########################
## Get SN and SMN tables ##
###########################

SMNs.ccss <- SMNs 
AnySN.ccss <- ANY_SNs

## Malignant status from Qi
malignantStatus <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/ExportedCCSS_data_update_malignant.txt", header = T, stringsAsFactors = F)
malignantStatus$KEY <- paste0(malignantStatus$ccssid, ":", malignantStatus$count)


SMNs.ccss$KEY <- paste0(SMNs.ccss$ccssid, ":", SMNs.ccss$count)
table(SMNs.ccss$KEY %in% malignantStatus$KEY)

SMNs.ccss$d_candx <- as.Date(SMNs.ccss$d_candx, format = "%d%b%Y")
SMNs.ccss$KEY <- paste0(SMNs.ccss$ccssid, ":", SMNs.ccss$d_candx)




library(readxl)
## Read file from Kyla
kyla.status <- read_excel("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Kyla/combinedsn_final_02_17_2023.xlsx")
kyla.status$d_candx <- as.numeric(kyla.status$d_candx)
kyla.status$d_candx <- as.Date(kyla.status$d_candx, origin = "1899-12-30") 
# kyla.status$KEY <- paste0(kyla.status$ccssid, ":", kyla.status$n_cond)
kyla.status$KEY <- paste0(kyla.status$ccssid, ":", kyla.status$d_candx)

table(SMNs.ccss$KEY %in% kyla.status$KEY)

###############
## CCSS SMNs ##
###############

# NMSC
NMSCs.ccss <- NMSCs
NMSCs.ccss$d_candx <- as.Date(NMSCs.ccss$d_candx, format = "%d%b%Y")
# NMSCs$KEY <- paste0(NMSCs$ccssid, ":", NMSCs$count)
NMSCs.ccss$KEY <- paste0(NMSCs.ccss$ccssid, ":", NMSCs.ccss$d_candx)
NMSCs.ccss <- cbind(diag=kyla.status$candxo3[match(NMSCs.ccss$KEY, kyla.status$KEY)],NMSCs.ccss)

# Breast
BREASTcancer.ccss <- BREASTcancer
BREASTcancer.ccss$d_candx <- as.Date(BREASTcancer.ccss$d_candx, format = "%d%b%Y")

# MENINGIOMA
MENINGIOMA.ccss <- MENINGIOMA
MENINGIOMA.ccss$d_candx <- as.Date(MENINGIOMA.ccss$d_candx, format = "%d%b%Y")

# Thyroid
THYROIDcancer.ccss <- THYROIDcancer 
THYROIDcancer.ccss$d_candx <- as.Date(THYROIDcancer.ccss$d_candx, format = "%d%b%Y")

SARCOMA.ccss <- SARCOMA
SARCOMA.ccss$d_candx <- as.Date(SARCOMA.ccss$d_candx, format = "%d%b%Y")
# SARCOMA$KEY <- paste0(SARCOMA$ccssid, ":", SARCOMA$count)
SARCOMA.ccss$KEY <- paste0(SARCOMA.ccss$ccssid, ":", SARCOMA.ccss$d_candx)
SARCOMA.ccss <- cbind(diag=kyla.status$candxo3[match(SARCOMA.ccss$KEY, kyla.status$KEY)],SARCOMA.ccss)




SMNs.ccss <- cbind(diag=kyla.status$candxo3 [match(SMNs.ccss$KEY, kyla.status$KEY)], SMNs.ccss)
SMNs.ccss <- cbind(diaggrp=SMNs.ccss$groupdx3, SMNs.ccss)
# cc <- cbind.data.frame(SMNs.ccss$candxo3, SMNs.ccss$groupdx3)
# table(SMNs.ccss$ccssid %in% to.remove)

SMNs.ccss$diag[is.na(SMNs.ccss$diag)] <- SMNs.ccss$groupdx3[is.na(SMNs.ccss$diag)]


check.sarcoma <- SMNs.ccss[grepl("Sarcoma", SMNs.ccss$groupdx3, ignore.case = T),]


SMNs.ccss$new.diagnosis <- NA

# SMNs.ccss$new.diagnosis[grepl("basal cell|squamous cell|squamous", SMNs.ccss$diag, ignore.case = T)] <- "NMSCs"
SMNs.ccss$new.diagnosis[grepl("basal cell|squamous cell|squamous", SMNs.ccss$diag, ignore.case = T)] <- "Other"
SMNs.ccss$new.diagnosis[grepl("breast", SMNs.ccss$diaggrp, ignore.case = T)] <- "Breast"
SMNs.ccss$new.diagnosis[grepl("thyroid", SMNs.ccss$diaggrp, ignore.case = T)] <- "Thyroid"
SMNs.ccss$new.diagnosis[grepl("meningioma", SMNs.ccss$diaggrp, ignore.case = T)] <- "Meningioma"
SMNs.ccss$new.diagnosis[grepl("Sarcoma", SMNs.ccss$diaggrp, ignore.case = T)] <- "Sarcoma"
SMNs.ccss$new.diagnosis[is.na(SMNs.ccss$new.diagnosis)] <- "Other"

other.smn.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis=="Other",]
malignant.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis!="Other",]

SMNs.ccss$new.diagnosis.cleaned <- NA
SMNs.ccss$new.diagnosis.cleaned[grepl("basal cell", SMNs.ccss$diag, ignore.case = T)] <- "Other"
SMNs.ccss$new.diagnosis.cleaned[grepl("squamous", SMNs.ccss$diag, ignore.case = T)] <- "Other"
SMNs.ccss$new.diagnosis.cleaned[grepl("melanoma", SMNs.ccss$diag, ignore.case = T)] <- "Malignant melanoma"

SMNs.ccss$new.diagnosis.cleaned[grepl("ductal|Mammary Carcinoma|Infiltrating duct carcinoma|Infiltrating duct", SMNs.ccss$diag, ignore.case = T)] <- "Infiltrating ductal carcinoma of breast"
SMNs.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Secretory", SMNs.ccss$diag, ignore.case = T)] <- "Secretory carcinoma of breast"
SMNs.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Metaplastic, Breast", SMNs.ccss$diag, ignore.case = T)] <- "Metaplastic carcinoma of breast"
SMNs.ccss$new.diagnosis.cleaned[grepl("Cribriform carcinoma|Cribriform", SMNs.ccss$diag, ignore.case = T)] <- "Cribriform carcinoma of breast"
SMNs.ccss$new.diagnosis.cleaned[grepl("Lobular carcinoma,", SMNs.ccss$diag, ignore.case = T)] <- "Lobular carcinoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Mucinous adenocarcinoma", SMNs.ccss$diag, ignore.case = T)] <- "Mucinous adenocarcinoma of breast"
SMNs.ccss$new.diagnosis.cleaned[grepl("Medullary carcinoma", SMNs.ccss$diag, ignore.case = T)] <- "Medullary carcinoma of breast"


cc <- SMNs.ccss[SMNs.ccss$new.diagnosis == "Breast",]

other.smn.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis=="Other",]
malignant.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis!="Other",]
malignant.ccss <- cbind.data.frame(new.diagnosis.cleaned=malignant.ccss$new.diagnosis.cleaned, malignant.ccss)

SMNs.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Papillary|Carcinoma, Hurthle Cell and Papillary|Papillary carcinoma", SMNs.ccss$diag, ignore.case = T)] <- "Carcinoma, Papillary, Thyroid"
SMNs.ccss$new.diagnosis.cleaned[grepl("Adenocarcinoma, Papillary", SMNs.ccss$diag, ignore.case = T)] <- "Adenocarcinoma, Papillary, Thyroid"
SMNs.ccss$new.diagnosis.cleaned[grepl("Papillary adenocarcinoma", SMNs.ccss$diag, ignore.case = T)] <- "Adenocarcinoma, Papillary, Thyroid"
SMNs.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Follicular, Thyroid|Follicular carcinoma|Follicular", SMNs.ccss$diag, ignore.case = T)] <- "Adenocarcinoma, Follicular, Thyroid"
SMNs.ccss$new.diagnosis.cleaned[grepl("Microcarcinoma, Papillary|Papillary microcarcinoma", SMNs.ccss$diag, ignore.case = T)] <- "Microcarcinoma, Papillary, Thyroid"
SMNs.ccss$new.diagnosis.cleaned[grepl("Oxyphilic adenocarcinoma", SMNs.ccss$diag, ignore.case = T)] <- "Oxyphilic Adenocarcinoma, Thyroid"


SMNs.ccss$new.diagnosis.cleaned[grepl("meningioma|Meningeal sarcomatosis", SMNs.ccss$diag, ignore.case = T)] <- "Meningioma"

SMNs.ccss$new.diagnosis.cleaned[grepl("Osteosarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Osteosarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Chondrosarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Chondrosarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Dermatofibrosarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Dermatofibrosarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Leiomyosarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Leiomyosarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Ewing", SMNs.ccss$diag, ignore.case = T)] <- "Ewing's Sarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Histiocytic Sarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Histiocytic Sarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Epithelioid", SMNs.ccss$diag, ignore.case = T)] <- "Epithelioid Sarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Liposarcoma|Lipoarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Liposarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Haemangiosarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Haemangiosarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Fibrosarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Fibrosarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Kaposi sarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Kaposi sarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Giant cell sarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Giant cell sarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Haemangiosarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Haemangiosarcoma"

SMNs.ccss$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Sarcoma, High Grade Spindle Cell", SMNs.ccss$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
SMNs.ccss$new.diagnosis.cleaned[grepl("Sarcoma, High Grade, Adrenal Gland", SMNs.ccss$diag, ignore.case = T)] <- "Sarcoma, High Grade Adrenal Gland"
SMNs.ccss$new.diagnosis.cleaned[grepl("Sarcoma, Soft Tissue", SMNs.ccss$diag, ignore.case = T)] <- "Sarcoma, Soft Tissue"
SMNs.ccss$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", SMNs.ccss$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"

SMNs.ccss$diaggrp[grepl("breast", SMNs.ccss$diaggrp, ignore.case = T)] <- "Breast"
SMNs.ccss$diaggrp[grepl("meningioma", SMNs.ccss$diaggrp, ignore.case = T)] <- "Meningioma"
SMNs.ccss$diaggrp[grepl("thyroid", SMNs.ccss$diaggrp, ignore.case = T)] <- "Thyroid"
SMNs.ccss$diaggrp[grepl("Sarcoma", SMNs.ccss$diaggrp, ignore.case = T)] <- "Sarcoma"



other.smn.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis=="Other",]
malignant.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis!="Other",]
malignant.ccss <- cbind.data.frame(new.diagnosis.cleaned=malignant.ccss$new.diagnosis.cleaned, malignant.ccss)


# SMNs.ccss$diag[grepl("sarcoma" , SMNs.ccss$diag, ignore.case = T)][!SMNs.ccss$diag[grepl("sarcoma" , SMNs.ccss$diag, ignore.case = T)] %in% SMNs.ccss$diag[grepl("sarcoma" , SMNs.ccss$new.diagnosis.cleaned, ignore.case = T)]]
cc <- SMNs.ccss[SMNs.ccss$new.diagnosis == "Breast",]

SMNs.ccss$new.diagnosis.cleaned2 <- SMNs.ccss$new.diagnosis.cleaned
SMNs.ccss$new.diagnosis.cleaned2[grepl("Sarcoma", SMNs.ccss$diaggrp, ignore.case = T)& grepl("^Malignant|Dermatofibro|Ganglioma|papillary|firbromatosis|fibromyxoid|Giant cell|Neurilemoma|lipoma|Haemangioblastoma|comedocarcinoma", SMNs.ccss$diag, ignore.case = T)] <- "Other"
SMNs.ccss$new.diagnosis.cleaned2[grepl("Malignant melanoma", SMNs.ccss$new.diagnosis.cleaned, ignore.case = T)] <- "Other"
SMNs.ccss$new.diagnosis.cleaned2[is.na(SMNs.ccss$new.diagnosis.cleaned2)] <- "Other"

other.smn.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis=="Other",]
malignant.ccss <- SMNs.ccss[SMNs.ccss$new.diagnosis!="Other",]
malignant.ccss <- cbind.data.frame(new.diagnosis.cleaned2=malignant.ccss$new.diagnosis.cleaned2, malignant.ccss)

# Malignant Peripheral Nerve Sheath Tumor, Neck, Left

other.SMN.ccss <- SMNs.ccss[(SMNs.ccss$new.diagnosis.cleaned2=="Other"),]
table(ifelse(grepl("Yes", other.SMN.ccss$seersmn), "M", "B"))
# M 
# 342

SMNs.ccss.table <- as.data.frame(table(SMNs.ccss$new.diagnosis.cleaned2, SMNs.ccss$diaggrp))
SMNs.ccss.table <- SMNs.ccss.table[SMNs.ccss.table$Freq !=0,]

View(table(SMNs.ccss$new.diagnosis.cleaned2))

#############
## Any SNs ##
#############


AnySN.ccss <- ANY_SNs

AnySN.ccss$KEY <- paste0(AnySN.ccss$ccssid, ":", AnySN.ccss$count)
table(AnySN.ccss$KEY %in% malignantStatus$KEY)

AnySN.ccss$d_candx <- as.Date(AnySN.ccss$d_candx, format = "%d%b%Y")
AnySN.ccss$KEY <- paste0(AnySN.ccss$ccssid, ":", AnySN.ccss$d_candx)


AnySN.ccss <- cbind(diag=kyla.status$candxo3 [match(AnySN.ccss$KEY, kyla.status$KEY)], AnySN.ccss)
AnySN.ccss <- cbind(diaggrp=AnySN.ccss$groupdx3, AnySN.ccss)
# cc <- cbind.data.frame(AnySN.ccss$candxo3, AnySN.ccss$groupdx3)
# table(AnySN.ccss$ccssid %in% to.remove)

AnySN.ccss$diag[is.na(AnySN.ccss$diag)] <- AnySN.ccss$groupdx3[is.na(AnySN.ccss$diag)]


check.sarcoma <- AnySN.ccss[grepl("Sarcoma", AnySN.ccss$groupdx3, ignore.case = T),]


AnySN.ccss$new.diagnosis <- NA

# AnySN.ccss$new.diagnosis[grepl("basal cell|squamous cell|squamous", AnySN.ccss$diag, ignore.case = T)] <- "NMSCs"
AnySN.ccss$new.diagnosis[grepl("basal cell|squamous cell|squamous", AnySN.ccss$diag, ignore.case = T)] <- "Other"
AnySN.ccss$new.diagnosis[grepl("breast", AnySN.ccss$diaggrp, ignore.case = T)] <- "Breast"
AnySN.ccss$new.diagnosis[grepl("thyroid", AnySN.ccss$diaggrp, ignore.case = T)] <- "Thyroid"
AnySN.ccss$new.diagnosis[grepl("meningioma", AnySN.ccss$diaggrp, ignore.case = T)] <- "Meningioma"
AnySN.ccss$new.diagnosis[grepl("Sarcoma", AnySN.ccss$diaggrp, ignore.case = T)] <- "Sarcoma"
AnySN.ccss$new.diagnosis[is.na(AnySN.ccss$new.diagnosis)] <- "Other"

other.sn.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis=="Other",]
malignant.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis!="Other",]

AnySN.ccss$new.diagnosis.cleaned <- NA
AnySN.ccss$new.diagnosis.cleaned[grepl("basal cell", AnySN.ccss$diag, ignore.case = T)] <- "Other"
AnySN.ccss$new.diagnosis.cleaned[grepl("squamous", AnySN.ccss$diag, ignore.case = T)] <- "Other"
AnySN.ccss$new.diagnosis.cleaned[grepl("melanoma", AnySN.ccss$diag, ignore.case = T)] <- "Malignant melanoma"

AnySN.ccss$new.diagnosis.cleaned[grepl("ductal|Mammary Carcinoma|Infiltrating duct carcinoma|Infiltrating duct", AnySN.ccss$diag, ignore.case = T)] <- "Infiltrating ductal carcinoma of breast"
AnySN.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Secretory", AnySN.ccss$diag, ignore.case = T)] <- "Secretory carcinoma of breast"
AnySN.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Metaplastic, Breast", AnySN.ccss$diag, ignore.case = T)] <- "Metaplastic carcinoma of breast"
AnySN.ccss$new.diagnosis.cleaned[grepl("Cribriform carcinoma|Cribriform", AnySN.ccss$diag, ignore.case = T)] <- "Cribriform carcinoma of breast"
AnySN.ccss$new.diagnosis.cleaned[grepl("Lobular carcinoma,", AnySN.ccss$diag, ignore.case = T)] <- "Lobular carcinoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Mucinous adenocarcinoma", AnySN.ccss$diag, ignore.case = T)] <- "Mucinous adenocarcinoma of breast"
AnySN.ccss$new.diagnosis.cleaned[grepl("Medullary carcinoma", AnySN.ccss$diag, ignore.case = T)] <- "Medullary carcinoma of breast"


cc <- AnySN.ccss[AnySN.ccss$new.diagnosis == "Breast",]

other.sn.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis=="Other",]
malignant.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis!="Other",]
malignant.ccss <- cbind.data.frame(new.diagnosis.cleaned=malignant.ccss$new.diagnosis.cleaned, malignant.ccss)

AnySN.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Papillary|Carcinoma, Hurthle Cell and Papillary|Papillary carcinoma", AnySN.ccss$diag, ignore.case = T)] <- "Carcinoma, Papillary, Thyroid"
AnySN.ccss$new.diagnosis.cleaned[grepl("Adenocarcinoma, Papillary", AnySN.ccss$diag, ignore.case = T)] <- "Adenocarcinoma, Papillary, Thyroid"
AnySN.ccss$new.diagnosis.cleaned[grepl("Papillary adenocarcinoma", AnySN.ccss$diag, ignore.case = T)] <- "Adenocarcinoma, Papillary, Thyroid"
AnySN.ccss$new.diagnosis.cleaned[grepl("Carcinoma, Follicular, Thyroid|Follicular carcinoma|Follicular", AnySN.ccss$diag, ignore.case = T)] <- "Adenocarcinoma, Follicular, Thyroid"
AnySN.ccss$new.diagnosis.cleaned[grepl("Microcarcinoma, Papillary|Papillary microcarcinoma", AnySN.ccss$diag, ignore.case = T)] <- "Microcarcinoma, Papillary, Thyroid"
AnySN.ccss$new.diagnosis.cleaned[grepl("Oxyphilic adenocarcinoma", AnySN.ccss$diag, ignore.case = T)] <- "Oxyphilic Adenocarcinoma, Thyroid"


AnySN.ccss$new.diagnosis.cleaned[grepl("meningioma|Meningeal sarcomatosis", AnySN.ccss$diag, ignore.case = T)] <- "Meningioma"

AnySN.ccss$new.diagnosis.cleaned[grepl("Osteosarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Osteosarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Chondrosarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Chondrosarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Dermatofibrosarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Dermatofibrosarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Leiomyosarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Leiomyosarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Ewing", AnySN.ccss$diag, ignore.case = T)] <- "Ewing's Sarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Histiocytic Sarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Histiocytic Sarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Epithelioid", AnySN.ccss$diag, ignore.case = T)] <- "Epithelioid Sarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Liposarcoma|Lipoarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Liposarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Haemangiosarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Haemangiosarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Fibrosarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Fibrosarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Kaposi sarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Kaposi sarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Giant cell sarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Giant cell sarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Haemangiosarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Haemangiosarcoma"

AnySN.ccss$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Sarcoma, High Grade Spindle Cell", AnySN.ccss$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"
AnySN.ccss$new.diagnosis.cleaned[grepl("Sarcoma, High Grade, Adrenal Gland", AnySN.ccss$diag, ignore.case = T)] <- "Sarcoma, High Grade Adrenal Gland"
AnySN.ccss$new.diagnosis.cleaned[grepl("Sarcoma, Soft Tissue", AnySN.ccss$diag, ignore.case = T)] <- "Sarcoma, Soft Tissue"
AnySN.ccss$new.diagnosis.cleaned[grepl("Spindle Cell Sarcoma", AnySN.ccss$diag, ignore.case = T)] <- "Spindle Cell Sarcoma"

AnySN.ccss$diaggrp[grepl("breast", AnySN.ccss$diaggrp, ignore.case = T)] <- "Breast"
AnySN.ccss$diaggrp[grepl("meningioma", AnySN.ccss$diaggrp, ignore.case = T)] <- "Meningioma"
AnySN.ccss$diaggrp[grepl("thyroid", AnySN.ccss$diaggrp, ignore.case = T)] <- "Thyroid"
AnySN.ccss$diaggrp[grepl("Sarcoma", AnySN.ccss$diaggrp, ignore.case = T)] <- "Sarcoma"



other.sn.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis=="Other",]
malignant.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis!="Other",]
malignant.ccss <- cbind.data.frame(new.diagnosis.cleaned=malignant.ccss$new.diagnosis.cleaned, malignant.ccss)


# AnySN.ccss$diag[grepl("sarcoma" , AnySN.ccss$diag, ignore.case = T)][!AnySN.ccss$diag[grepl("sarcoma" , AnySN.ccss$diag, ignore.case = T)] %in% AnySN.ccss$diag[grepl("sarcoma" , AnySN.ccss$new.diagnosis.cleaned, ignore.case = T)]]
cc <- AnySN.ccss[AnySN.ccss$new.diagnosis == "Breast",]

AnySN.ccss$new.diagnosis.cleaned2 <- AnySN.ccss$new.diagnosis.cleaned
AnySN.ccss$new.diagnosis.cleaned2[grepl("Sarcoma", AnySN.ccss$diaggrp, ignore.case = T)& grepl("^Malignant|Dermatofibro|Ganglioma|papillary|firbromatosis|fibromyxoid|Giant cell|Neurilemoma|lipoma|Haemangioblastoma|comedocarcinoma", AnySN.ccss$diag, ignore.case = T)] <- "Other"
AnySN.ccss$new.diagnosis.cleaned2[grepl("Malignant melanoma", AnySN.ccss$new.diagnosis.cleaned, ignore.case = T)] <- "Other"
AnySN.ccss$new.diagnosis.cleaned2[is.na(AnySN.ccss$new.diagnosis.cleaned2)] <- "Other"

other.AnySN.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis=="Other",]
malignant.ccss <- AnySN.ccss[AnySN.ccss$new.diagnosis!="Other",]
malignant.ccss <- cbind.data.frame(new.diagnosis.cleaned2=malignant.ccss$new.diagnosis.cleaned2, malignant.ccss)

# Malignant Peripheral Nerve Sheath Tumor, Neck, Left
other.AnySN.ccss <- AnySN.ccss[(AnySN.ccss$new.diagnosis.cleaned2=="Other"),]
table(ifelse(grepl("Yes", other.SMN.ccss$seersmn), "M", "B"))
# B   M 
# 727 303 

AnySN.ccss.table <- as.data.frame(table(AnySN.ccss$new.diagnosis.cleaned2, AnySN.ccss$diaggrp))
AnySN.ccss.table <- AnySN.ccss.table[AnySN.ccss.table$Freq !=0,]

View(table(AnySN.ccss$new.diagnosis.cleaned2))



# Malignant Peripheral Nerve Sheath Tumor, Neck, Left

other.AnySN.ccss <- AnySN.ccss[(AnySN.ccss$new.diagnosis.cleaned=="Other"),]
table(ifelse(grepl("Yes", other.AnySN.ccss$seersmn), "M", "B"))
# M 
# 96

AnySN.ccss$new.diagnosis.cleaned.malignant <- AnySN.ccss$new.diagnosis.cleaned2
AnySN.ccss$new.diagnosis.cleaned.malignant <- paste0(AnySN.ccss$new.diagnosis.cleaned.malignant, "_", ifelse(grepl("Yes", AnySN.ccss$seersmn), "M", "B"))
# AnySN.ccss$new.diagnosis.cleaned.malignant[grepl("NA_M|NA_B",AnySN.ccss$new.diagnosis.cleaned.malignant)] <- NA


View(table(AnySN.ccss$new.diagnosis.cleaned.malignant))
View(table(AnySN.ccss$new.diagnosis.cleaned2))

# cc <- as.data.frame(table(AnySN.ccss$new.diagnosis.cleaned2))
# 
# ## How many of AnySN.ccss are NMSCs
AnySN.ccss.NMSC <- AnySN.ccss[AnySN.ccss$KEY %in% NMSCs.ccss$KEY,]

AnySN.ccss.table <- as.data.frame(table(AnySN.ccss$new.diagnosis.cleaned, AnySN.ccss$diaggrp))
AnySN.ccss.table <- AnySN.ccss.table[AnySN.ccss.table$Freq !=0,]

other.AnySN.ccss <- AnySN.ccss[(AnySN.ccss$new.diagnosis.cleaned2=="Other"),]
other.AnySN.ccss <- other.AnySN.ccss[!other.AnySN.ccss$KEY %in% AnySN.ccss.NMSC$KEY]
# sum(grepl("Basal", other.AnySN.ccss$diag))

cc <- AnySN.ccss[grepl("Basal", AnySN.ccss$diag),]

# save(other.SMN.sjlife, other.anySN.sjlife, other.SMN.ccss, other.AnySN.ccss, 
#     file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/other_SN_and_SMN_objects.RData")

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/other_SN_and_SMN_objects.RData")

## Create table for other SNs

sub.other.any.sn.sjlife <- cbind.data.frame(KEY=other.anySN.sjlife$KEY,sjlid=other.anySN.sjlife$sjlid, diag=other.anySN.sjlife$diag, 
                                            diaggrp=other.anySN.sjlife$diaggrp, icdo3morph=other.anySN.sjlife$icdo3morph,
                                            icdo3behavior=other.anySN.sjlife$icdo3behavior)

sub.other.smn.sjlife <- cbind.data.frame(KEY=other.SMN.sjlife$KEY, sjlid=other.SMN.sjlife$sjlid, diag=other.SMN.sjlife$diag, 
                                            diaggrp=other.SMN.sjlife$diaggrp, icdo3morph=other.SMN.sjlife$icdo3morph, 
                                            icdo3behavior=other.SMN.sjlife$icdo3behavior)


merged.sjlife.other <- merge(sub.other.any.sn.sjlife, sub.other.smn.sjlife, by="KEY", all = T, suffixes = c(".SN", ".SMN"))

# merged.sjlife.other$diag.SN[grepl("Dermatofibrosarcoma", merged.sjlife.other$diag.SN, ignore.case = T)] <- "Malignant peripheral nerve sheath tumor"
merged.sjlife.other$diag.SMN[grepl("Dermatofibrosarcoma", merged.sjlife.other$diag.SMN, ignore.case = T)] <- "Dermatofibroma"

merged.sjlife.other$icdo3morph.SMN[grepl("sarcoma", merged.sjlife.other$icdo3morph.SMN, ignore.case = T)] <- merged.sjlife.other$diag.SMN[grepl("sarcoma", merged.sjlife.other$icdo3morph.SMN, ignore.case = T)]
merged.sjlife.other$icdo3morph.SN[grepl("Sarcoma", merged.sjlife.other$icdo3morph.SN, ignore.case = F)] <- merged.sjlife.other$diag.SN[grepl("Sarcoma", merged.sjlife.other$icdo3morph.SN, ignore.case = F)]
merged.sjlife.other$icdo3morph.SN[grepl("Sarcoma, Soft Tissue, Undifferentiated", merged.sjlife.other$icdo3morph.SN, ignore.case = T)] <- "Malignant peripheral nerve sheath tumor"
merged.sjlife.other$icdo3morph.SMN[grepl("Sarcoma, Soft Tissue, Undifferentiated", merged.sjlife.other$icdo3morph.SMN, ignore.case = T)] <- "Malignant peripheral nerve sheath tumor"
merged.sjlife.other$icdo3morph.SN[grepl("Squamous cell carcinoma in situ", merged.sjlife.other$icdo3morph.SN, ignore.case = T)] <- merged.sjlife.other$diag.SN[grepl("Squamous cell carcinoma in situ", merged.sjlife.other$icdo3morph.SN, ignore.case = T)] 
merged.sjlife.other$icdo3morph.SMN[grepl("Squamous cell carcinoma in situ", merged.sjlife.other$icdo3morph.SMN, ignore.case = T)] <- merged.sjlife.other$diag.SMN[grepl("Squamous cell carcinoma in situ", merged.sjlife.other$icdo3morph.SMN, ignore.case = T)] 


merged.sjlife.other$icdo3morph.SN <- gsub("\\, NOS|NOS|lymphocyte-rich|Neck, Left|, Retroperitoneum|, nodular lymphocyte predominance|diffuse", "",merged.sjlife.other$icdo3morph.SN, ignore.case = T)
merged.sjlife.other$icdo3morph.SMN <- gsub("\\, NOS|NOS|lymphocyte-rich|Neck, Left|, Retroperitoneum|, nodular lymphocyte predominance|diffuse", "",merged.sjlife.other$icdo3morph.SMN, ignore.case = T)

merged.sjlife.other$icdo3morph.SN[grepl("Malignant peripheral", merged.sjlife.other$icdo3morph.SN, ignore.case = T)] <- "Malignant peripheral nerve sheath tumor"
merged.sjlife.other$icdo3morph.SMN[grepl("Malignant peripheral", merged.sjlife.other$icdo3morph.SMN, ignore.case = T)] <- "Malignant peripheral nerve sheath tumor"

merged.sjlife.other$icdo3morph.SN <- gsub(", $", "",merged.sjlife.other$icdo3morph.SN, ignore.case = T)
merged.sjlife.other$icdo3morph.SMN <- gsub(", $", "",merged.sjlife.other$icdo3morph.SMN, ignore.case = T)

cc.smn <- as.data.frame(table(merged.sjlife.other$icdo3morph.SMN))
cc.sn <- as.data.frame(table(merged.sjlife.other$icdo3morph.SN))

merged.sjlife.other$icdo3morph.SMN.B.M <- paste0(merged.sjlife.other$icdo3morph.SMN, "_", merged.sjlife.other$icdo3behavior.SMN)
merged.sjlife.other$icdo3morph.SN.B.M <- paste0(merged.sjlife.other$icdo3morph.SN, "_", merged.sjlife.other$icdo3behavior.SN)


cc <- merge (cc.smn, cc.sn, by = "Var1", all = T)
View(cc)
View(table(merged.sjlife.other$icdo3morph.SN.B.M))

merged.sjlife.other$icdo3morph.SN[!grepl("Malignant", merged.sjlife.other$icdo3behavior.SN, ignore.case = T)]

View(table(SMN=merged.sjlife.other$icdo3morph.SMN, SN=merged.sjlife.other$icdo3morph.SN))


## CCSS
sub.other.any.sn.ccss <- cbind.data.frame(KEY=other.AnySN.ccss$KEY,ccssid=other.AnySN.ccss$ccssid, diag=other.AnySN.ccss$diag, 
                                            diaggrp=other.AnySN.ccss$diaggrp, icdo3behavior=other.AnySN.ccss$seersmn)

sub.other.smn.ccss <- cbind.data.frame(KEY=other.SMN.ccss$KEY,ccssid=other.SMN.ccss$ccssid, diag=other.SMN.ccss$diag, 
                                          diaggrp=other.SMN.ccss$diaggrp, icdo3behavior=other.SMN.ccss$seersmn)


sub.other.any.sn.ccss$diag_clean <- gsub("^[0-9]+/[0-9]+\\s+", "", sub.other.any.sn.ccss$diag)  # Remove leading number/number  
sub.other.any.sn.ccss$diag_clean <- gsub("\\s*,?\\s*NOS.*", "", sub.other.any.sn.ccss$diag_clean)  # Remove everything after 'NOS' or ', NOS'  

sub.other.any.sn.ccss$diag_clean <- gsub("^[0-9]+/[0-9]+\\s+", "", sub.other.any.sn.ccss$diag)  # Remove leading number/number  
sub.other.any.sn.ccss$diag_clean <- gsub("\\s*,?\\s*NOS.*", "", sub.other.any.sn.ccss$diag_clean)  # Remove everything after 'NOS' or ', NOS'  
sub.other.any.sn.ccss$diag_clean <- gsub("\\s*\\(.*?\\)", "", sub.other.any.sn.ccss$diag_clean)  # Remove everything inside parentheses  



sub.other.smn.ccss$diag_clean <- gsub("^[0-9]+/[0-9]+\\s+", "", sub.other.smn.ccss$diag)  # Remove leading number/number  
sub.other.smn.ccss$diag_clean <- gsub("\\s*,?\\s*NOS.*", "", sub.other.smn.ccss$diag_clean)  # Remove everything after 'NOS' or ', NOS'  
sub.other.smn.ccss$diag_clean <- gsub("\\s*\\(.*?\\)", "", sub.other.smn.ccss$diag_clean)  # Remove everything inside parentheses  

# sub.other.any.sn.ccss$diag_clean[grepl("lobular", sub.other.any.sn.ccss$diag_clean, ignore.case = T)] <- "Carcinoma"
# sub.other.smn.ccss$diag_clean[grepl("lobular", sub.other.smn.ccss$diag_clean, ignore.case = T)] <- "Carcinoma"

sub.other.any.sn.ccss$diag_clean[grepl("Dermatofibrosarcoma", sub.other.any.sn.ccss$diag_clean, ignore.case = T)] <- "Neurilemoma"
sub.other.smn.ccss$diag_clean[grepl("Dermatofibrosarcoma", sub.other.smn.ccss$diag_clean, ignore.case = T)] <- "Malignant peripheral nerve sheath tumor"

sub.other.any.sn.ccss$diag_clean[grepl("Basal|Squamous", sub.other.any.sn.ccss$diag_clean, ignore.case = T)] <- sub.other.any.sn.ccss$diaggrp[grepl("Basal|Squamous", sub.other.any.sn.ccss$diag_clean, ignore.case = T)]
sub.other.smn.ccss$diag_clean[grepl("Basal|Squamous", sub.other.smn.ccss$diag_clean, ignore.case = T)] <-  sub.other.any.smn.ccss$diaggrp[grepl("Basal|Squamous", sub.other.any.smn.ccss$diag_clean, ignore.case = T)]

sub.other.any.sn.ccss$diag_clean[grepl("8340.2000000000007", sub.other.any.sn.ccss$diag_clean, ignore.case = T)] <- "Unknown diagnosis"
sub.other.smn.ccss$diag_clean[grepl("8340.2000000000007", sub.other.smn.ccss$diag_clean, ignore.case = T)] <-  "Unknown diagnosis"

View(table(SN=sub.other.any.sn.ccss$diag_clean, SMN=sub.other.any.sn.ccss$diag_clean))

count.sn <- as.data.frame(table(sub.other.any.sn.ccss$diag_clean))
count.smn <- as.data.frame(table(sub.other.smn.ccss$diag_clean))

View(table(SN=sub.other.any.sn.ccss$diag_clean, ICDO=sub.other.any.sn.ccss$icdo3behavior))

cc <-merge (count.smn, count.sn, by = "Var1", all = T, suffixes = c("SMN", "SN"))

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
toyadav<- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/toyadav.sas7bdat")
toyadav$age_at_last_contact <- as.numeric(difftime(toyadav$lstcondt, toyadav$dob, units = "days")) / 365.25
toyadav$SNdt <- as.numeric(toyadav$diagdt - toyadav$dob) / 365.25

cc <- cbind.data.frame(sjlid=toyadav$sjlid, sndx=toyadav$sndx, agedx=toyadav$agedx, d_entry=toyadav$d_entry, sn=toyadav$sn, agefup = toyadav$agefup, ageend=toyadav$ageend, age_at_last_contact=toyadav$age_at_last_contact, age.at.SN=toyadav$SNdt, dob=toyadav$dob, lstcondt=toyadav$lstcondt)
unique_sn_data <- toyadav %>%
  select(sjlid, sndx, age_at_last_contact) %>%
  distinct()
median(unique_sn_data$age_at_last_contact, na.rm = TRUE)

cc$agefupAN <- as.numeric(difftime(cc$lstcondt, cc$dob, units = "days")) / 365.25  # Convert to years

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

## Age at follow up
ccss$count <- as.numeric(ccss$count)
ccss$agelstcontact[!is.na(ccss$AGE.ANY_SN)] <- ccss$AGE.ANY_SN[!is.na(ccss$AGE.ANY_SN)]


unique_samples <- ccss %>%
  distinct(ccssid, .keep_all = TRUE) %>%  # Keep only unique samples based on sjlid
  select(ccssid, agelstcontact)  # Select relevant columns

# Step 2: Calculate median follow-up age and IQR
median_agefup <- median(unique_samples$agelstcontact, na.rm = TRUE)  # Calculate median
iqr_agefup <- IQR(unique_samples$agelstcontact, na.rm = TRUE)  # Calculate IQR

# Display results
median_agefup
iqr_agefup



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








# sj.cc <- cbind.data.frame(sj$sjlid, sj$agedx, sj$agelstcontact, sj$sncount, sj$agedx)
data <- sj
data$AGE.ANY_SN <- time_length(interval(as.Date(data$dob), as.Date(data$gradedt)), "years")
## Age at last contact for cases is SN diagnosis data
# data$agelstcontact[!is.na(data$AGE.ANY_SN)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]
data$fuyears <- data$agelstcontact - data$AGE.ANY_SN
data$fuyears[is.na(data$fuyears)] <- data$agelstcontact[is.na(data$fuyears)]

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

# final$sjlid[duplicated(final$sjlid)]
# cc.final <- final[, c(1,(ncol(final)-11):ncol(final))]
# cc.data <- data[, c(1,(ncol(data)-11):ncol(data))]

minimum  <-  min(final$start, na.rm = TRUE) - 1
if (minimum < 0) minimum <- 0
maximum  <-  max(final$end, na.rm = TRUE) + 1
SNs_py <-  survSplit(final, cut=seq(minimum, maximum, 1), end="end",start="start",event="event") 
table(SNs_py$event)   ## Double check event numebr is correct
#### If you need the rows for first event analysis, take evt1=1
length(unique(SNs_py$sjlid[SNs_py$event==1]))
length(unique(SNs_py$sjlid[SNs_py$event==1 & SNs_py$evt1==1 ]))

SNs_py$PY <- SNs_py$end-SNs_py$start

sum(SNs_py$PY)

pp <- data [is.na(data$sncount)| data$sncount==1,]

## Attained age
attained.age <- round(median(SNs_py$agelstcontact, na.rm= T),1)
attained.age.IQR <- paste0(unname(round((quantile(SNs_py$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
attained.age.IQR <- gsub(" ", "-", attained.age.IQR)
attained.age <- paste0(attained.age, " (", attained.age.IQR, ")")
# "36.2 (26.2-46.9)"


## Follow up age
follow.up <- round(median(data$fuyears, na.rm= T),1)
follow.up.IQR <- paste0(unname(round((quantile(data$fuyears, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
follow.up.IQR <- gsub(" ", "-", follow.up.IQR)
follow.up <- paste0(follow.up, " (", follow.up.IQR, ")")
follow.up
"24.2 (11.7-35.4)"

##########
## CCSS ##
##########

# sj.cc <- cbind.data.frame(sj$sjlid, sj$agedx, sj$agelstcontact, sj$sncount, sj$agedx)
data <- ccss
data$AGE.ANY_SN <- as.numeric(data$a_candx)

data$gradeage <- data$gradedt
data$gradedt <- as.Date(data$d_candx, format = "%d%b%Y")

data$dob <- data$gradedt - as.numeric(data$gradeage)
## Age at last contact for cases is SN diagnosis data
# data$agelstcontact[!is.na(data$AGE.ANY_SN)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

data$fuyears <- data$agelstcontact - data$AGE.ANY_SN
data$fuyears[is.na(data$fuyears)] <- data$agelstcontact[is.na(data$fuyears)]

pp <- cbind.data.frame(data$ccssid, data$fuyears, data$AGE.ANY_SN, data$agelstcontact, data$a_candx)

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
sum(SNs_py$PY)


## Attained age
attained.age <- round(median(SNs_py$agelstcontact, na.rm= T),1)
attained.age.IQR <- paste0(unname(round((quantile(SNs_py$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
attained.age.IQR <- gsub(" ", "-", attained.age.IQR)
attained.age <- paste0(attained.age, " (", attained.age.IQR, ")")
attained.age


data$count <- as.numeric(data$count)
pp <- data [is.na(data$count)| data$count==1,]

## Follow up age
follow.up <- round(median(data$fuyears, na.rm= T),1)
follow.up.IQR <- paste0(unname(round((quantile(data$fuyears, prob=c(.25,.5,.75), type=1, na.rm = T))[c(1,3)], 1)), collapse = " ")
follow.up.IQR <- gsub(" ", "-", follow.up.IQR)
follow.up <- paste0(follow.up, " (", follow.up.IQR, ")")
follow.up
# 28 (8.9-37.2)

######





























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



