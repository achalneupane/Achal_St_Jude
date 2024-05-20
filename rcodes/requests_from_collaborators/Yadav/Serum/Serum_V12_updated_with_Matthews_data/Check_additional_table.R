# # Yutaka's email on 04/14/2024:
# Would it be possible to make Table 1 without restricting to people free from CMP at the first visit?
# I just want to see the numbers in view of the case-cohort design.

library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")

source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/get_matching_rows_from_CTCAE.R")
# df <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
# dim(df) # 20174

## add TB from Mathew
# TB <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/ForAchal_Survivors.txt", header = T, sep = "\t") ## March version
TB <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/Achal_survivors_04.04.2024.txt", header = T, sep = "\t") ## Updated by Matt in April

df <- TB

df$ageatsample <- floor(df$ageatsample * 10) / 10
df.original <- df
df <- df %>%
  distinct()
dim(df)  # 20137

## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
SERUM <- df[grepl("Plasma", df$aliquot_type, ignore.case = T),]

SERUM.original <- SERUM
dim(SERUM.original)
# [1] 10383     6
# 10514  7 ## April version

## I see age at serum sample is <18 yrs. Could you identify the samples among 18
#or higher only? Everyone needs to be adults at serum sample.
SERUM <- SERUM[which(SERUM$ageatsample >= 18),]

# remove vial zero
SERUM <- SERUM[SERUM$num_vials > 0 ,]


## 171 cases to keep
all.wanted.df.1200.to.update <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v11_output//plasma_data_batch1_1200_samples.txt", header = T, sep = "\t")
keep.171.cases <- all.wanted.df.1200.to.update$sjlid[all.wanted.df.1200.to.update$selection_group == "171_CMP_cases"]
keep.171.cases <- SERUM[SERUM$sjlid %in% keep.171.cases,]
SERUM <- SERUM[!SERUM$sjlid %in% keep.171.cases$sjlid,] # exclude 171 cases
dim(SERUM)
# 8536    7
## Now keep only those that have more than one vial and are alive
SERUM <- SERUM[SERUM$num_vials > 1,]
SERUM <- rbind.data.frame(SERUM, keep.171.cases)
dim(SERUM)
# [1] 7948    6


dim(SERUM)
# 8879
# 7948    7 ## April 4 version

## Keep max num vials and alive over dead
SERUM <- SERUM %>%
  group_by(sjlid, ageatsample) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(vitalstatus == "Alive" | all(vitalstatus == "Deceased"))

dim(SERUM)
# [1] 7934    6

## Also removing samples that have duplicate ageevent in CTCAE data
# SERUM <- SERUM[!SERUM$sjlid %in% c("SJL1225801", "SJL1265801", "SJL1430801", "SJL4730101", "SJL5134305", "SJL5146506"),]
dim(SERUM)
# 7934    6
## Could you please work on it by merging with the CTCAE grades for cardiomyopathy from the most recent data freeze? Please look at serum and plasma separately.
# read CTCAE
# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
CTCAE <- CTCAE[grepl("Cardiomyopathy", CTCAE$condition),]
# CTCAE.original <- CTCAE
CTCAE <- CTCAE.original

CTCAE <- CTCAE.original[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
# CTCAE.original.2 <- CTCAE.original[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
## Since Trans-omics CMP ageatsample is in one decimal, I am coverting CTCAE age also to one decimal place.
CTCAE$ageevent <- round(CTCAE$ageevent,1)

# # CTCAE.cc <- CTCAE[grepl("SJL0253301|SJL1063101", CTCAE$sjlid),]
# # Condition 1: Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher. If there are two ageevent that are miniumum value, we still apply this filter (but this did not make any difference)
# CTCAE <- CTCAE %>%
#   group_by(sjlid) %>%
#   filter(
#     !(grade[which(ageevent == min(ageevent))[1]] != 0)
#   ) %>%
#   ungroup() %>%
#   arrange(sjlid)
dim(CTCAE)
# 9218    9
# 9218    9

# V10 & V11
# 9218    9

CTCAE <- CTCAE[CTCAE$sjlid %in% unique(SERUM$sjlid),]
dim(CTCAE)
# v12 7148

# CTCAE <- get_rows_with_smaller_sample_age(CTCAE, SERUM, 0)
CTCAE <- get_rows_with_smaller_sample_age.all(CTCAE, SERUM, 7)
dim(CTCAE)
# V12
dim(CTCAE)
# 8090   13

## Add event number
CTCAE <- CTCAE %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid)  # Restore the original order

## first event should have sample
CTCAE <- CTCAE[!(CTCAE$grade == 0 & is.na(CTCAE$Sample_age)),]
CTCAE <- CTCAE[CTCAE$grade != -9,]


## Remove rows once grades 2 or higher are seen in ordered df by sjlid and event_number
CTCAE <- filter_rows_by_condition(CTCAE, "sjlid", "grade")
sum(CTCAE$grade !=0 & CTCAE$event_number==1)
table(CTCAE$grade)
## V12
# 0    2    3    5 
# 4540  200   53    1 

CTCAE.2 <- CTCAE


## Add event number
CTCAE.2 <- CTCAE.2 %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(new_event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid)  # Restore the original order

# # SJL5553107
# # Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher
# CTCAE.2 <- CTCAE.2 %>%
#   dplyr::group_by(sjlid) %>%
#   dplyr::filter(
#     !(grade[which(ageevent == min(ageevent))[1]] != 0)
#   ) %>%
#   dplyr::ungroup()


## Add event number
CTCAE.2 <- CTCAE.2 %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(new_event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid)  # Restore the original order


CTCAE.2$grade_2_or_higher <- ifelse(CTCAE.2$grade >= 2, "grade_2_or_higher", "grade_0")
table(CTCAE.2$event_number,CTCAE.2$grade)
## V12
# 0    2    3    5
# 1 2549    0    0    0
# 2 1286   72   22    0
# 3  508   49    9    0
# 4  151   19    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

table(CTCAE.2$event_number,CTCAE.2$grade_2_or_higher)
## V12
# grade_0 grade_2_or_higher
# 1    2549                 0
# 2    1286                94
# 3     508                58
# 4     151                27
# 5      29                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0

table(CTCAE.2$new_event_number,CTCAE.2$grade)
## V12
# 0    2    3    5
# 1 3173    0    0    0
# 2 1038   92   28    0
# 3  275   42   10    1
# 4   52   10    2    0
# 5    2    1    0    0

table(CTCAE.2$new_event_number,CTCAE.2$grade_2_or_higher)
## V12
# grade_0 grade_2_or_higher
# 1    3173                 0
# 2    1038               120
# 3     275                53
# 4      52                12
# 5       2                 1


####################################
## Email from Yadav on 04/30/2024 ##
####################################
# Yutaka and I discussed about these samples and would like to ask you provide the following variables for 1100 survivors you have identified:
# Age at primary cancer diagnosis
# Year of cancer diagnosis
# Age at plasma sample
# Cancer diagnosis type
# An indicator variable showing if they had cardiomyopathy (Grade 2 or higher) at baseline plasma sample
# Yes/no variable indicating if they were exposed to anthracyclines or heart RT
# Category (incident cardiomyopathy, HL or not HL survivors)

## Adding requested variables
# Age at primary cancer diagnosis & Year of cancer diagnosis
all.wanted.df.1200.to.update <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_data_batch1_1200_samples.txt", header = T, sep = "\t")
diagnosis <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat')
diagnosis$KEY <- paste0(diagnosis$sjlid, "_", diagnosis$diaggrp)
all.wanted.df.1200.to.update$KEY <- paste0(all.wanted.df.1200.to.update$sjlid, "_", all.wanted.df.1200.to.update$diaggrp)
all.wanted.df.1200.to.update$agedx <- diagnosis$agedx[match(all.wanted.df.1200.to.update$KEY, diagnosis$KEY)]
all.wanted.df.1200.to.update$diagdt <- diagnosis$diagdt[match(all.wanted.df.1200.to.update$KEY, diagnosis$KEY)]


chemo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat")
chemo <- chemo[c("sjlid", "anthracyclines_dose_any",	"anthracyclines_dose_prim",	"anthracyclines_dose_5",	"anthracyclines_dose_10")]

radiation <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat")
radiation <- radiation[c("sjlid", "maxchestrtdose")]

all.wanted.df.1200.to.update <- cbind.data.frame(all.wanted.df.1200.to.update, chemo[match(all.wanted.df.1200.to.update$sjlid, chemo$sjlid), "anthracyclines_dose_any"])
all.wanted.df.1200.to.update <- cbind.data.frame(all.wanted.df.1200.to.update, radiation[match(all.wanted.df.1200.to.update$sjlid, radiation$sjlid),"maxchestrtdose"])

all.wanted.df.1200.to.update$anthra_or_chestRT <- ifelse((all.wanted.df.1200.to.update$anthracyclines_dose_any > 0 | all.wanted.df.1200.to.update$maxchestrtdose > 200), "Yes", "No")

  
table(CTCAE.2$new_event_number,CTCAE.2$grade)
sum(CTCAE.2$new_event_number >= 2 & CTCAE.2$grade >=2)
# 186
remove.sjl <- CTCAE.2$sjlid[CTCAE.2$new_event_number >= 2 & CTCAE.2$grade >=2]
CTCAE.3 <- CTCAE.2[!CTCAE.2$sjlid %in% remove.sjl,]
table(CTCAE.3$new_event_number,CTCAE.3$grade)
# 0    2    3    4
# 1 2987  262  125    4
# 2  972    0    0    0
# 3  262    0    0    0
# 4   51    0    0    0
# 5    2    0    0    0
ge.grade2.baseline <- unique(CTCAE.3$sjlid[CTCAE.3$new_event_number == 1 & CTCAE.3$grade >=2])
length(ge.grade2.baseline)
# 391

all.wanted.df.1200.to.update$CMP_at_baseline [all.wanted.df.1200.to.update$sjlid %in% ge.grade2.baseline] <- "Yes"
all.wanted.df.1200.to.update$CMP_at_baseline[all.wanted.df.1200.to.update$CMP_status == "Yes"] <- "No"

write.table(all.wanted.df.1200.to.update, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma.df.1100.with_additional_variables_toYadav.txt", col.names = T, row.names = F, sep = "\t", quote = F)



diagnosis <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat')
diagnosis <- diagnosis[which(diagnosis$primdx ==1),] # primary diagnosis
diagnosis$KEY <- paste0(diagnosis$sjlid, "_", diagnosis$diaggrp)
# diagnosis$dups <- duplicated(diagnosis$KEY)
SERUM.not.CA171 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/SERUM.not.CA171.txt", sep = "\t", header = T)
# SERUM.not.CA171$diagrp2 <- diagnosis$diaggrp[match(SERUM.not.CA171$sjlid, diagnosis$sjlid)]

SERUM.not.CA171$agedx <- diagnosis$agedx[match(SERUM.not.CA171$sjlid, diagnosis$sjlid)]
SERUM.not.CA171$diagdt <- diagnosis$diagdt[match(SERUM.not.CA171$sjlid, diagnosis$sjlid)]

SERUM.not.CA171$CMP_at_baseline [SERUM.not.CA171$sjlid %in% ge.grade2.baseline] <- "Yes"
# SERUM.not.CA171$CMP_at_baseline[is.na(SERUM.not.CA171$CMP_at_baseline)] <- "No"

SERUM.not.CA171 <- cbind.data.frame(SERUM.not.CA171, chemo[match(SERUM.not.CA171$sjlid, chemo$sjlid), "anthracyclines_dose_any"])
SERUM.not.CA171 <- cbind.data.frame(SERUM.not.CA171, radiation[match(SERUM.not.CA171$sjlid, radiation$sjlid),"maxchestrtdose"])

SERUM.not.CA171$anthra_or_chestRT <- ifelse((SERUM.not.CA171$anthracyclines_dose_any > 0 | SERUM.not.CA171$maxchestrtdose > 200), "Yes", "No")

SERUM.not.CA171 <- SERUM.not.CA171[c("tb_number", "sjlid", "num_vials", "ageevent", "vitalstatus", "agedx", "diagdt", "Sample_age", "diaggrp", "CMP_at_baseline", "anthra_or_chestRT")]

all.hodgkin <- SERUM.not.CA171[grepl("^Hodgkin", SERUM.not.CA171$diaggrp, ignore.case = T),]
all.non.hodgkin <- SERUM.not.CA171[!SERUM.not.CA171$sjlid %in% c(all.hodgkin$sjlid),] # exclude those in CA.171 + 200 hodgkins samples from the original table
all.non.hodgkin <- all.non.hodgkin[all.non.hodgkin$diaggrp!="",]

write.table(all.hodgkin, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/hodgkin_all.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(all.non.hodgkin, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/nonhodgkin_all.txt", col.names = T, row.names = F, sep = "\t", quote = F)
