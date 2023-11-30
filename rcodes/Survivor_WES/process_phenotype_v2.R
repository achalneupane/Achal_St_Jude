library(dplyr)
library(tidyr)
library(haven)

## Phenotype for all samples sequenced
all.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/sample_mapping.txt", header = F)

# CTCAE.data <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
sum(is.na(CTCAE.data$grade))
# 8075
# remove missing and unknown grades
CTCAE.data <- CTCAE.data[!is.na(CTCAE.data$grade),]
CTCAE.data <- CTCAE.data[-which(CTCAE.data$grade==-9),]

cc <- CTCAE.data[1:30, 1:11]


# Filter rows where all grades are zero: if all grades are zero, get the row with the max gradedt (i.e., latest grade date)
result_zero <- CTCAE.data %>%
  group_by(sjlid, condition) %>%
  filter(all(grade == 0)) %>%
  slice(which.max(gradedt)) %>%
  ungroup()

# Filter rows where any grade is greater than 0: 
## If the any grade is 1 or higher, then get the row with max grade (Also note, if there are multiple grades >=1, get the row with first gradedt with max grades) 
result_non_zero <- CTCAE.data %>%
  group_by(sjlid, condition) %>%
  filter(any(grade > 0)) %>%
  filter(grade == max(grade)) %>%
  slice(which.min(gradedt)) %>%
  ungroup()

# Combine the results
CTCAE.data.1 <- bind_rows(result_zero, result_non_zero)

###############################
## Convert it to wide format ##
###############################
CTCAE.data.2 <- CTCAE.data.1 %>%
  select(MRN, sjlid, studypop, sjlife_cohort, gender, condition, grade, ageevent) 

# Define the columns to keep at the beginning
first_columns <- c("MRN", "sjlid", "studypop", "sjlife_cohort", "gender")

# Pivot the dataframe to wide format
CTCAE.data.3 <- CTCAE.data.2 %>%
  pivot_wider(
    id_cols = c(MRN, sjlid, studypop, sjlife_cohort, gender),
    names_from = condition,
    values_from = c(grade, ageevent),
    names_glue = "{condition}_{.value}"
  ) %>%
  select(c(all_of(first_columns), sort(names(.))))

## Display the wide format result
# CTCAE.data.3

CTCAE.data.4 <- CTCAE.data.3 %>%
  mutate(across(contains("_grade"), 
                ~case_when(
                  . == 0 ~ "0",
                  . >= 2 ~ "1",
                  TRUE ~ NA_character_
                ),
                .names = "{str_remove(.col, '_grade')}_status_0_vs_2")) %>%
  select(c(all_of(first_columns), sort(names(.))))


CTCAE.data.4 <- CTCAE.data.4 %>%
  mutate(across(contains("_grade"), 
                ~case_when(
                  . == 0 ~ "0",
                  . >= 3 ~ "1",
                  TRUE ~ NA_character_
                ),
                .names = "{str_remove(.col, '_grade')}_status_0_vs_3")) %>%
  select(c(all_of(first_columns), sort(names(.))))
