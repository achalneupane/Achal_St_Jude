library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")

source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/get_matching_rows_from_CTCAE.R")
df <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
dim(df) # 20174
df$ageatsample <- floor(df$ageatsample * 10) / 10
df.original <- df
df <- df %>%
  distinct()
dim(df)  # 20137

# ## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
# PLASMA <- df[grepl("Plasma", df$aliquot_type, ignore.case = T),]
# # PLASMA.within.7.days.age.eevnt <- get_matching_rows(PLASMA, CTCAE, 7/365.25) # within 7 days
# # PLASMA.within.7.days.age.eevnt <- PLASMA.within.7.days.age.eevnt[!is.na(PLASMA.within.7.days.age.eevnt$grade),]
# PLASMA <- get_matching_rows(PLASMA, CTCAE, 7) # on same day; use 7/365.25 for a 7 days window
# PLASMA <- PLASMA[!is.na(PLASMA$grade) & PLASMA$grade != -9,]

## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
SERUM <- df[grepl("Serum", df$aliquot_type, ignore.case = T),]
SERUM.original <- SERUM
## Keep max num vials and alive over dead
SERUM <- SERUM %>%
  group_by(sjlid, ageatsample) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(vitalstatus == "Alive" | all(vitalstatus == "Deceased"))

## Could you please work on it by merging with the CTCAE grades for cardiomyopathy from the most recent data freeze? Please look at serum and plasma separately.
# read CTCAE
# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
CTCAE <- CTCAE[grepl("Cardiomyopathy", CTCAE$condition),]
CTCAE.original <- CTCAE
CTCAE <- CTCAE.original[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
# CTCAE.original.2 <- CTCAE.original[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
## Since Trans-omics CMP ageatsample is in one decimal, I am coverting CTCAE age also to one decimal place.
CTCAE$ageevent <- round(CTCAE$ageevent,1)

# CTCAE.cc <- CTCAE[grepl("SJL0253301|SJL1063101", CTCAE$sjlid),]
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
CTCAE <- get_rows_with_smaller_sample_age.all(CTCAE, SERUM, 7)
dim(CTCAE)
# 8332

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

CTCAE.2 <- CTCAE


# SJL5553107
# Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher
CTCAE.2 <- CTCAE.2 %>%
  dplyr::group_by(sjlid) %>%
  dplyr::filter(
    !(grade[which(ageevent == min(ageevent))[1]] != 0)
  ) %>%
  dplyr::ungroup()


## Add event number
CTCAE.2 <- CTCAE.2 %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(new_event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid)  # Restore the original order


CTCAE.2$grade_2_or_higher <- ifelse(CTCAE.2$grade >= 2, "grade_2_or_higher", "grade_0")
table(CTCAE.2$event_number,CTCAE.2$grade)
table(CTCAE.2$event_number,CTCAE.2$grade_2_or_higher)

table(CTCAE.2$new_event_number,CTCAE.2$grade)
table(CTCAE.2$new_event_number,CTCAE.2$grade_2_or_higher)

CTCAE.3 <- CTCAE.2[!is.na(CTCAE.2$Sample_age),]
table(CTCAE.3$event_number,CTCAE.3$grade_2_or_higher)


table(is.na(CTCAE.2$Sample_age) & CTCAE.2$grade >=2)
table(CTCAE.2$grade >=2)

write.table(CTCAE.2, "Serum_data_processed_v6.txt", col.names = T, row.names = F, sep = "\t")



