library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")

source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/get_matching_rows_from_CTCAE_cleaned.R")
df <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
dim(df) # 20174
# df$ageatsample <- floor(df$ageatsample * 10) / 10
df.original <- df
df <- df %>%
  distinct()
dim(df)  # 20137


## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
SERUM <- df[grepl("Serum", df$aliquot_type, ignore.case = T),]
SERUM.original <- SERUM
## Keep max num vials and alive over dead
SERUM <- SERUM %>%
  group_by(sjlid, ageatsample) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(vitalstatus == "Alive" | all(vitalstatus == "Deceased"))

## work on it by merging with the CTCAE grades for cardiomyopathy from the most recent data freeze.
# read CTCAE
# CTCAE.original <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
CTCAE <- CTCAE.original[grepl("Cardiomyopathy", CTCAE.original$condition),]
# CTCAE.original <- CTCAE
CTCAE <- CTCAE[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
## Since Trans-omics CMP ageatsample is in one decimal, I am coverting CTCAE age also to one decimal place.
# CTCAE$ageevent <- round(CTCAE$ageevent,1)

# Condition 1: Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher. If there are two ageevent that are miniumum value, we still apply this filter (but this did not make any difference)
CTCAE <- CTCAE %>%
  group_by(sjlid) %>%
  filter(
    !(grade[which(ageevent == min(ageevent))[1]] != 0)
  ) %>%
  ungroup() %>%
arrange(sjlid)
dim(CTCAE)
# 9218    9

CTCAE <- CTCAE[CTCAE$sjlid %in% unique(SERUM$sjlid),]
dim(CTCAE)
# 8332

# CTCAE <- get_rows_with_smaller_sample_age(CTCAE, SERUM, 0)
## If we use 21 days (0.05 years); we still don't have "SJL5309110" "SJL4836507". serum for SJL5309110 was collected after ageevent and some (SJL4836507, SJL5087912) were not in Kyla's serum list
CTCAE <- get_rows_with_smaller_sample_age.all(CTCAE, SERUM, 21) 
dim(CTCAE)
# 8332

## Order by agevent and add event number
CTCAE <- CTCAE %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid)  # Restore the original order

dim(CTCAE)
# 8332
## Remove rows once grades 2 or higher are seen in ordered df by sjlid and event_number
CTCAE.2 <- filter_rows_by_condition(CTCAE, "sjlid", "grade")
sum(CTCAE.2$grade !=0 & CTCAE.2$event_number==1)
table(CTCAE.2$grade)

dim(CTCAE.2)
# 8124

## baseline should have sample
CTCAE.2 <- CTCAE.2[!(CTCAE.2$grade == 0 & is.na(CTCAE.2$Sample_age)),]
CTCAE.2 <- CTCAE.2[CTCAE.2$grade != -9,]
dim(CTCAE.2)

# SJL5553107
# Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher
CTCAE.2 <- CTCAE.2 %>%
  dplyr::group_by(sjlid) %>%
  dplyr::filter(
    !(grade[which(ageevent == min(ageevent))[1]] != 0)
  ) %>%
  dplyr::ungroup()


dim(CTCAE.2)
# 6194


## Add event number
CTCAE.2 <- CTCAE.2 %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(new_event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid)  # Restore the original order


table(CTCAE.2$new_event_number,CTCAE.2$grade)
table(CTCAE.2$new_event_number,CTCAE.2$grade >= 2)

# Calculate age different from first visit to first CMP
result <- calculate_ageevent_difference(CTCAE.2)
print(result)


write.table(CTCAE.2, "Serum_data_processed_v7.txt", col.names = T, row.names = F, sep = "\t")



