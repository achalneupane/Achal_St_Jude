## Yadav's email on 10/23/2023 
# Could you provide the summary as in “Sheet1” for the following subset of survivors?

# Exposed to either anthracyclines or chest radiation
# With WGS (no restriction on treatment exposures)
# Exposed to either anthracyclines or chest radiation with WGS


library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")
CTCAE.serum <- read.table("Serum_data_processed_v6.txt", header = T, sep = "\t")

demographic <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
chemo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat")
chemo <- chemo[c("sjlid", "anthracyclines_dose_any",	"anthracyclines_dose_prim",	"anthracyclines_dose_5",	"anthracyclines_dose_10")]

radiation <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat")
radiation <- radiation[c("sjlid", "maxchestrtdose")]
radiation.2 <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiationsum_yn.sas7bdat")
radiation.2 <- radiation.2[c("sjlid", "Chest", "Chest_prim", "Chest_5", "Chest_10")]

have.WGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr8.preQC_biallelic_renamed_ID_edited.vcf.gz.fam")

# demographic <- cbind.data.frame(demographic, chemo[match(demographic$sjlid, chemo$sjlid),])
demographic <- cbind.data.frame(demographic, chemo[match(demographic$sjlid, chemo$sjlid), "anthracyclines_dose_any"])
# demographic <- cbind.data.frame(demographic, radiation[match(demographic$sjlid, radiation$sjlid),])
demographic <- cbind.data.frame(demographic, radiation[match(demographic$sjlid, radiation$sjlid),"maxchestrtdose"])
# demographic <- cbind.data.frame(demographic, radiation.2[match(demographic$sjlid, radiation.2$sjlid),])
demographic <- cbind.data.frame(demographic, radiation.2[match(demographic$sjlid, radiation.2$sjlid),c("Chest", "Chest_prim", "Chest_5", "Chest_10")])
demographic$WGS <- ifelse(demographic$sjlid %in% have.WGS$V2, "Yes", 'No')

demographic <- demographic[match(CTCAE.serum$sjlid, demographic$sjlid),]
demographic <- demographic[!colnames(demographic) %in% colnames(CTCAE.serum)]

CTCAE.serum <- cbind.data.frame(CTCAE.serum,demographic)
table(CTCAE.serum$new_event_number, CTCAE.serum$grade)

# Exposed to either anthracyclines or chest radiation
CTCAE.serum.1 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 | CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade)
table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade >= 2)

## exposed to anthracyclines with chest RT
CTCAE.serum.4 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade)
table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade >= 2)

## exposed to anthracyclines without chest RT
CTCAE.serum.5 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose <= 200,]
table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade)
table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade >= 2)

## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0,]
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade)
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade >= 2)

## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0,]
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade)
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade >= 2)


## break down exposed to anthracyclines regardless of chest RT
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="White"] <- "White"
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="Black"] <- "Black"
CTCAE.serum.6$race_group[is.na(CTCAE.serum.6$race_group)] <- "Other"
table(CTCAE.serum.6$race_group)

## Only in White
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="White",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)

## Only in Black
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Black",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)

## Only in other races
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Other",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)

######################################################
## With WGS (no restriction on treatment exposures) ##
######################################################
CTCAE.serum.2 <- CTCAE.serum[CTCAE.serum$WGS == "Yes",]
table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade)
table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade >= 2)
# Exposed to either anthracyclines or chest radiation with WGS
CTCAE.serum.3 <- CTCAE.serum.2[CTCAE.serum.2$anthracyclines_dose_any > 0 | CTCAE.serum.2$maxchestrtdose > 200,]
table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade)
table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade >= 2)






# # CTCAE.original <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
# CTCAE.original.2 <- CTCAE.original[grepl("Cardiomyopathy", CTCAE.original$condition),]
# CTCAE.original.2 <- CTCAE.original.2[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
# # Also, we previously performed proteomics/metabolomics experiments on 200
# # survivors (100 with cardiomyopathy); they are attached here (see Sheet3 for
# # relevant information). Could you check how many of these overlap with the
# # current ~3500 or so (based on grade, age at assessment and sample age)? Also,
# # check if these 200 samples fulfill the same criteria applied to the 3500
# # samples (and if not, how many meet and how many do not).
# library(readxl)
# # Read data from Sheet3
# proteomics <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/pilot_study_200_samples.xlsx", sheet = "Sheet3")
# colnames(proteomics) <- c("sjlid", "grade", "ageevent", "TBID", "Sample_age", "status")
# sum(proteomics$sjlid %in% CTCAE.serum$sjlid)
# # 153
# 
# 
# # Define a function to check if two values are within X-days window
# within_window <- function(value1, value2) {
#   return(abs(value1 - value2) <= 20/365.25)
# }
# 
# # Perform the merge on sjlid, grade, and ageevent
# merged_data <- merge(proteomics, CTCAE.serum, by = c("sjlid", "grade"), suffixes = c(".x", ".y"))
# 
# # Subset the merged data based on your criteria
# overlap_data <- subset(merged_data, within_window(Sample_age.x, Sample_age.y) & within_window(ageevent.x, ageevent.y))
# 
# 
# # Overlaps
# dim(overlap_data)
# # 57
# 
# proteomics.ca <- proteomics[proteomics$status == 1,]
# proteomics.co <- proteomics[proteomics$status == 0,]
# 
# sum(proteomics.ca$sjlid %in% CTCAE.serum$sjlid) ## 59
# sum(proteomics.co$sjlid %in% CTCAE.serum$sjlid) ## 98
# 
# SERUM.kyla <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
# SERUM.kyla <- SERUM.kyla[SERUM.kyla$aliquot_type == "Serum",]
# sum(proteomics.ca$sjlid %in% SERUM.kyla$sjlid) ## 100
# sum(proteomics.co$sjlid %in% SERUM.kyla$sjlid) ## 99
# 
# 
# table(CTCAE.serum$new_event_number,CTCAE.serum$grade)
# table(CTCAE.serum$new_event_number,CTCAE.serum$grade >= 2)
# 
# # Calculate age different from first visit to first CMP
# result <- calculate_ageevent_difference(CTCAE.serum)
# print(result)


## Read data from excel to exclude previous study
previous.study <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/Serum_data_processed_v6_corrected_10_5_2023.xlsx", sheet = "Serum_data_processed_v6")

# Yadav On 11/28/2023:
# Can you please add the summary Table for survivors exposed to anthracyclines (with or without chest RT)? No need to restrict to with WGS only.
CTCAE.serum.saved <- CTCAE.serum
CTCAE.serum <- cbind.data.frame(previous.study, CTCAE.serum[match(previous.study$sjlid, CTCAE.serum$sjlid),c("race", "anthracyclines_dose_any", "maxchestrtdose", "Chest", "Chest_prim", "Chest_5", "Chest_10", "WGS")])
CTCAE.serum$new_event_number <- CTCAE.serum$Order
## Exclude previous ones
CTCAE.serum <- CTCAE.serum[!grepl("Yes", CTCAE.serum$Excluded),]

table(CTCAE.serum$new_event_number, CTCAE.serum$grade)

# Exposed to either anthracyclines or chest radiation
CTCAE.serum.1 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 | CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade)
table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade >= 2)

## exposed to anthracyclines with chest RT
CTCAE.serum.4 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade)
table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade >= 2)

## exposed to anthracyclines without chest RT
CTCAE.serum.5 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose <= 200,]
table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade)
table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade >= 2)

## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0,]
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade)
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade >= 2)


## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="White"] <- "White"
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="Black"] <- "Black"
CTCAE.serum.6$race_group[is.na(CTCAE.serum.6$race_group)] <- "Other"
table(CTCAE.serum.6$race_group)

## Only in White
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="White",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)

## Only in Black
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Black",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)

## Only in other races
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Other",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)

######################################################
## With WGS (no restriction on treatment exposures) ##
######################################################
CTCAE.serum.2 <- CTCAE.serum[CTCAE.serum$WGS == "Yes",]
table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade)
table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade >= 2)
# Exposed to either anthracyclines or chest radiation with WGS
CTCAE.serum.3 <- CTCAE.serum.2[CTCAE.serum.2$anthracyclines_dose_any > 0 | CTCAE.serum.2$maxchestrtdose > 200,]
table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade)
table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade >= 2)
