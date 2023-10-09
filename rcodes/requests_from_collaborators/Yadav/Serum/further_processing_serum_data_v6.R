## Yadav's email on 10/23/2023 
# Could you provide the summary as in “Sheet1” for the following subset of survivors?

# Exposed to either anthracyclines or chest radiation
# With WGS (no restriction on treatment exposures)
# Exposed to either anthracyclines or chest radiation with WGS


library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")
CTCAE.serum <- read.table("Serum_data_processed_v7.txt", header = T, sep = "\t")

demographic <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
chemo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat")
chemo <- chemo[c("sjlid", "anthracyclines_dose_any",	"anthracyclines_dose_prim",	"anthracyclines_dose_5",	"anthracyclines_dose_10")]

radiation <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat")
radiation <- radiation[c("sjlid", "maxchestrtdose")]
radiation.2 <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiationsum_yn.sas7bdat")
radiation.2 <- radiation.2[c("sjlid", "Chest", "Chest_prim", "Chest_5", "Chest_10")]

have.WGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr8.preQC_biallelic_renamed_ID_edited.vcf.gz.fam")

demographic <- cbind.data.frame(demographic, chemo[match(demographic$sjlid, chemo$sjlid),])
demographic <- cbind.data.frame(demographic, radiation[match(demographic$sjlid, radiation$sjlid),])
demographic <- cbind.data.frame(demographic, radiation.2[match(demographic$sjlid, radiation.2$sjlid),])
demographic$WGS <- ifelse(demographic$sjlid %in% have.WGS$V2, "Yes", 'No')

CTCAE.serum <- cbind.data.frame(CTCAE.serum,demographic[match(CTCAE.serum$sjlid, demographic$sjlid),])
# Exposed to either anthracyclines or chest radiation
# With WGS (no restriction on treatment exposures)
# Exposed to either anthracyclines or chest radiation with WGS



# CTCAE.original <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
CTCAE.original.2 <- CTCAE.original[grepl("Cardiomyopathy", CTCAE.original$condition),]
CTCAE.original.2 <- CTCAE.original.2[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
# Also, we previously performed proteomics/metabolomics experiments on 200
# survivors (100 with cardiomyopathy); they are attached here (see Sheet3 for
# relevant information). Could you check how many of these overlap with the
# current ~3500 or so (based on grade, age at assessment and sample age)? Also,
# check if these 200 samples fulfill the same criteria applied to the 3500
# samples (and if not, how many meet and how many do not).
library(readxl)
# Read data from Sheet3
proteomics <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/pilot_study_200_samples.xlsx", sheet = "Sheet3")
colnames(proteomics) <- c("sjlid", "grade", "ageevent", "TBID", "Sample_age", "status")
sum(proteomics$sjlid %in% CTCAE.serum$sjlid)
# 157


# Define a function to check if two values are within X-days window
within_window <- function(value1, value2) {
  return(abs(value1 - value2) <= 20/365.25)
}

# Perform the merge on sjlid, grade, and ageevent
merged_data <- merge(proteomics, CTCAE.serum, by = c("sjlid", "grade"), suffixes = c(".x", ".y"))

# Subset the merged data based on your criteria
overlap_data <- subset(merged_data, within_window(Sample_age.x, Sample_age.y) & within_window(ageevent.x, ageevent.y))


# Overlaps
dim(overlap_data)
# 57

proteomics.ca <- proteomics[proteomics$status == 1,]
proteomics.co <- proteomics[proteomics$status == 0,]

sum(proteomics.ca$sjlid %in% CTCAE.serum$sjlid) ## 59
sum(proteomics.co$sjlid %in% CTCAE.serum$sjlid) ## 98

SERUM.kyla <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
SERUM.kyla <- SERUM.kyla[SERUM.kyla$aliquot_type == "Serum",]
sum(proteomics.ca$sjlid %in% SERUM.kyla$sjlid) ## 100
sum(proteomics.co$sjlid %in% SERUM.kyla$sjlid) ## 99


table(CTCAE.serum$new_event_number,CTCAE.serum$grade)
table(CTCAE.serum$new_event_number,CTCAE.serum$grade >= 2)

# Calculate age different from first visit to first CMP
result <- calculate_ageevent_difference(CTCAE.serum)
print(result)


