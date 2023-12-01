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
                  . > 0 ~ "1",
                  TRUE ~ NA_character_
                ),
                .names = "{str_remove(.col, '_grade')}_status_gt_0")) %>%
  select(c(all_of(first_columns), sort(names(.))))


CTCAE.data.4 <- CTCAE.data.4 %>%
  mutate(across(contains("_grade"), 
                ~case_when(
                  . == 0 ~ "0",
                  . > 2 ~ "1",
                  TRUE ~ NA_character_
                ),
                .names = "{str_remove(.col, '_grade')}_status_gt_2")) %>%
  select(c(all_of(first_columns), sort(names(.))))


CTCAE.data.4 <- CTCAE.data.4 %>%
  mutate(across(contains("_grade"), 
                ~case_when(
                  . == 0 ~ "0",
                  . > 3 ~ "1",
                  TRUE ~ NA_character_
                ),
                .names = "{str_remove(.col, '_grade')}_status_gt_3")) %>%
  select(c(all_of(first_columns), sort(names(.))))

###########################################
## Add clinical and demographic features ##
###########################################
lstcondt <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data/lstcondt.sas7bdat')
CTCAE.data.4$agelstcont <- lstcondt$agelstcontact[match(CTCAE.data.4$sjlid, lstcondt$sjlid)]

demographics <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat')
# CTCAE.data.4$gender1 <- demographics$gender[match(CTCAE.data.4$sjlid, demographics$sjlid)]
CTCAE.data.4$race <- demographics$race[match(CTCAE.data.4$sjlid, demographics$sjlid)]

diagnosis <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat')
CTCAE.data.4$agedx <- diagnosis$agedx[match(CTCAE.data.4$sjlid, diagnosis$sjlid)]

chemo <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat')
radiation <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat')


## MaxsegRT
CTCAE.data.4$maxsegrtdose <- radiation$maxsegrtdose[match(CTCAE.data.4$sjlid, radiation$sjlid)]

## Maxchest
CTCAE.data.4$maxchestrtdose <- radiation$maxchestrtdose[match(CTCAE.data.4$sjlid, radiation$sjlid)]

## Maxneck
CTCAE.data.4$maxneckrtdose <- radiation$maxneckrtdose[match(CTCAE.data.4$sjlid, radiation$sjlid)]

## Maxabdomen
CTCAE.data.4$maxabdrtdose <- radiation$maxabdrtdose[match(CTCAE.data.4$sjlid, radiation$sjlid)]

## Maxpelvis
CTCAE.data.4$maxpelvisrtdose <- radiation$maxpelvisrtdose[match(CTCAE.data.4$sjlid, radiation$sjlid)]

## aa_class_dose_5 **
CTCAE.data.4$aa_class_dose_5 <- chemo$alkylating_dose_5[match(CTCAE.data.4$sjlid, chemo$sjlid)]
CTCAE.data.4$alkylating_dose_any <- chemo$alkylating_dose_any[match(CTCAE.data.4$sjlid, chemo$sjlid)]

## anthra_jco_dose_5 **
CTCAE.data.4$anthra_jco_dose_5 <- chemo$anthracyclines_dose_5[match(CTCAE.data.4$sjlid, chemo$sjlid)]
CTCAE.data.4$anthracyclines_dose_any <- chemo$anthracyclines_dose_any[match(CTCAE.data.4$sjlid, chemo$sjlid)]

## epitxn_dose_5 **
CTCAE.data.4$epitxn_dose_5 <- chemo$epipodophyllotoxins_dose_5[match(CTCAE.data.4$sjlid, chemo$sjlid)]
CTCAE.data.4$epipodophyllotoxins_dose_any <- chemo$epipodophyllotoxins_dose_any[match(CTCAE.data.4$sjlid, chemo$sjlid)]




# check
cc <- CTCAE.data.4[c(1:6,1256:1269, grep("Cardiomyopathy", colnames(CTCAE.data.4)))]

#############################################################################################
## Now check how many of the status columns have cases vs controls to association analysis ##
#############################################################################################


# Define the columns that represent grades
grade_columns <- CTCAE.data.4 %>%
  select(contains("status_")) %>%
  names()

# Initialize an empty list to store tables
table_list <- list()

# Loop through each column in grade_columns
for (col in grade_columns) {
  # Create a table for the current column
  tab <- table(CTCAE.data.4[[col]])
  
  # Ensure both 0 and 1 are present in the table
  if (!("0" %in% names(tab))) {
    tab[["0"]] <- 0
  }
  
  if (!("1" %in% names(tab))) {
    tab[["1"]] <- 0
  }
  
  # Rename the columns to 0 and 1
  names(tab) <- c("0", "1")
  
  # Convert the table to a data frame and add variable column
  tab_df <- data.frame(value = as.numeric(names(tab)), count = as.numeric(tab))
  tab_df$variable <- col
  
  # Append the data frame to the list
  table_list <- append(table_list, list(tab_df))
}

# Combine the list of tables into a single data frame
table_df <- do.call(rbind, table_list)

# Rename the columns
colnames(table_df) <- c("Status", "count", "variable")

# Display the result
print(table_df)


table_df <- table_df %>%
  pivot_wider(
    id_cols = variable,
    names_from = Status,
    values_from = count,
    names_prefix = "status_"
  )

# Keep status with at least 25 cases and controls
filtered_table_df <- table_df %>%
  filter(status_0 >= 25, status_1 >= 25)

