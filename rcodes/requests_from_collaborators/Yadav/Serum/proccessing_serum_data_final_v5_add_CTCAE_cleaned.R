library(dplyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")

## Could you please work on it by merging with the CTCAE grades for cardiomyopathy from the most recent data freeze? Please look at serum and plasma separately.
# read CTCAE
# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
CTCAE <- CTCAE[grepl("Cardiomyopathy", CTCAE$condition),]
# CTCAE.original <- CTCAE
CTCAE <- CTCAE.original[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
## Since Trans-omics CMP ageatsample is in one decimal, I am coverting CTCAE age also to one decimal place.
CTCAE$ageevent <- round(CTCAE$ageevent,1)

# CTCAE.cc <- CTCAE[grepl("SJL0253301|SJL1063101", CTCAE$sjlid),]
# Condition 1: Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher. If there are two ageevent that are miniumum value, we still apply this filter (but this did not make any difference)
CTCAE <- CTCAE %>%
  group_by(sjlid) %>%
  filter(
    !(grade[which(ageevent == min(ageevent))[1]] != 0)
  ) %>%
  ungroup()


source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/get_matching_rows_from_CTCAE.R")
df <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
df$ageatsample <- floor(df$ageatsample * 10) / 10
df.original <- df

## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
PLASMA <- df[grepl("Plasma", df$aliquot_type, ignore.case = T),]
# PLASMA.within.7.days.age.eevnt <- get_matching_rows(PLASMA, CTCAE, 7/365.25) # within 7 days
# PLASMA.within.7.days.age.eevnt <- PLASMA.within.7.days.age.eevnt[!is.na(PLASMA.within.7.days.age.eevnt$grade),]
PLASMA <- get_matching_rows(PLASMA, CTCAE, 7) # on same day; use 7/365.25 for a 7 days window
PLASMA <- PLASMA[!is.na(PLASMA$grade) & PLASMA$grade != -9,]

## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
SERUM <- df[grepl("Serum", df$aliquot_type, ignore.case = T),]
SERUM.original <- SERUM
SERUM <- get_matching_rows(SERUM, CTCAE, 7) # on same day
SERUM <- SERUM[!is.na(SERUM$grade) & SERUM$grade != -9,]

## Conditions::
# 1. Remove rows where ageevent is not greater than (or within 1 week of) sample age
# 2. within each sample if all rows have grade 2 or higher, or if minimum ageevent has grade 2 or higher, skip this sample. Do not extract any rows from that sample.
# 3. Withing each sample, remove rows where grades are smaller than the grades previously seen in the ordered rows.
# 4. Then, within each tb_number, keep rows with max num_vials; then keep ALIVE over DEAD (something like this: filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD")))

## SJL5192513, SJL1599510, SJL1659216, SJL1675616

##############################
# ## Count by SJLID samples
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/Count_functions_by_sjlid.R")

## Count by serum samples
# source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/Count_functions_by_tb_number.R")

########### 
## Serum ##
###########
df.serum <- SERUM
df.serum$Sample_age <- df.serum$ageatsample


## Processing...
# Condition 1: Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher. If there are two ageevent that are miniumum value, we still apply this filter (but this did not make any difference)
step1 <- df.serum %>%
  group_by(sjlid) %>%
  filter(
    !(grade[which(ageevent == min(ageevent))[1]] != 0)
  ) %>%
  ungroup()

dim(step1)
# anti_join(df.serum, step1)

# Condition 2: Proceed to remove rows where grades are smaller than the grades previously seen in the ordered rows. In simpler terms, we keeps the rows where the grade value is the highest seen so far within each group. This effectively retains only the rows with the highest grade value within each sjlid group.
step2 <- step1 %>% 
  arrange(sjlid, ageevent) %>% 
  group_by(sjlid) %>% 
  filter(cummax(grade) == grade) %>% 
  ungroup()

dim(step2)
step2$tb_number <- paste0(step2$sjlid, ":", step2$Sample_age)
# Condition 3: Keep rows with max num_vials and vitalstatus
step3 <- step2 %>%
  group_by(tb_number) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(vitalstatus == "Alive" | all(vitalstatus == "Deceased"))

dim(step3)
step3$tb_number[duplicated(step3$tb_number)] # check duplicates
# removed_rows <- anti_join(step2, step3)

# # # Remove rows with duplicate grade, ageevent, sample_age combinations (it only removes 23 rows; need to discuss with Yadav!)
step3 <- step3 %>%
  distinct(sjlid, grade, ageevent, Sample_age, .keep_all = TRUE)

dim(step3) # 5992   10
# removed_rows <- anti_join(step2, step3)

FINAL.1 <- step3
FINAL.SERUM <- step3




CTCAE.in.serum.with.1st.event.0 <-  CTCAE[CTCAE$sjlid %in% unique(FINAL.SERUM$sjlid ),]
table(CTCAE.in.serum.with.1st.event.0$grade)
# -9    0    2    3    4    5 
# 1 7518  180   97    1    3 

## Add event number
CTCAE.in.serum.with.1st.event.0 <- CTCAE.in.serum.with.1st.event.0 %>%
  group_by(sjlid) %>%
  arrange(ageevent) %>%
  mutate(event_number = row_number()) %>%
  ungroup() %>%
  arrange(sjlid)  # Restore the original order

table(CTCAE.in.serum.with.1st.event.0$event_number, CTCAE.in.serum.with.1st.event.0$grade)
# Ev# -9    0    2    3    4    5
# 1    0 4105    0    0    0    0
# 2    1 2075   82   27    0    0
# 3    0  880   58   22    0    0
# 4    0  321   26   16    1    2
# 5    0   83    8   15    0    0
# 6    0   32    2    6    0    1
# 7    0   15    2    5    0    0
# 8    0    6    1    3    0    0
# 9    0    1    1    3    0    0

## Check how many overlap from previous round of 165 cases in Kylas data
previous.164 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/serum_data_processed_final_v3_cleaned.txt", header = T)
## How many in total overlap
length(unique(SERUM.original$sjlid)) # The one Kyla gave
# 5479
length(unique(previous.164$sjlid))
# 3747

sum(unique(SERUM.original$sjlid) %in% unique(previous.164$sjlid))
# 3747

result.164 <- previous.164 %>%
  group_by(sjlid) %>%
  filter(any(grade == 0) & any(grade >= 2))

length(unique(FINAL.SERUM$sjlid))
## 4105
sum(unique(FINAL.SERUM$sjlid) %in% unique(result.164$sjlid))
# 146

## Add event number:
FINAL.SERUM <- FINAL.SERUM %>%
  group_by(sjlid) %>%
  arrange(ageevent) %>%
  mutate(event_number = row_number()) %>%
  ungroup() %>%
  arrange(sjlid)  # Restore the original order

table(FINAL.SERUM$event_number, FINAL.SERUM$grade)
#     0    2    3
# 1 4105    0    0
# 2 1301   59   19
# 3  384   26    8
# 4   73    9    4
# 5    3    0    1



## So now, we want to work with this dataset: CTCAE.in.serum.with.1st.event.0
table(CTCAE.in.serum.with.1st.event.0$grade)
CTCAE.in.serum.with.1st.event.0 <- CTCAE.in.serum.with.1st.event.0[CTCAE.in.serum.with.1st.event.0$grade != -9,] # remove grade -9
# grade grade 0 and the first cardiomypathy



# first_CMP_event <- CTCAE.in.serum.with.1st.event.0 %>%
#   group_by(sjlid) %>%
#   arrange(sjlid, event_number) %>%
#   filter(!cumany(grade >= 2) | row_number() <= min(which(grade >= 2))) %>%
#   ungroup()
## Note: the warnings basically means there are o grade 2 or greater for those sjlids



first_CMP_event <- filter_rows_by_condition(CTCAE.in.serum.with.1st.event.0, "sjlid", "grade")

## Add Sample_age
V5_FINAL_res <- get_rows_with_smaller_sample_age(first_CMP_event, SERUM.original, 7)
# If you don’t have a sample, just keep the first row with grade 2 or higher, if available.
V5_FINAL_res.missing <- V5_FINAL_res[(V5_FINAL_res$grade == 0 & is.na(V5_FINAL_res$Sample_age)),]
V5_FINAL_res <- V5_FINAL_res[!(V5_FINAL_res$grade == 0 & is.na(V5_FINAL_res$Sample_age)),]
# V5_FINAL_res <- rbind.data.frame(V5_FINAL_res, V5_FINAL_res.missing[V5_FINAL_res.missing$event_number == 1,])


V5_FINAL_res$grade_2_or_higher <- ifelse(V5_FINAL_res$grade >= 2, "grade_2_or_higher", "grade_0")



table(V5_FINAL_res$event_number,V5_FINAL_res$grade)
table(V5_FINAL_res$event_number,V5_FINAL_res$grade >= 2)

table(V5_FINAL_res$grade == 0 & V5_FINAL_res$event_number == 1)
table(first_CMP_event$grade == 0 & first_CMP_event$event_number == 1)

##########################################


























############ 
## Plasma ##
############
df <- PLASMA
df$Sample_age <- df$ageatsample

## Get counts in non-missing events
check_grades_eq_or_higher_than(df, 2)
check_grades_eq_or_higher_than(df, 3)
check_grades_eq_or_higher_than(df, 4)
check_grades_eq_or_higher_than(df, 5)

# counts for 0 to X transitions with any grade at min ageevent
check_grades_transition(df, 0, 2)
check_grades_transition(df, 0, 3)
check_grades_transition(df, 0, 4)
check_grades_transition(df, 0, 5)

check_grades_eq_or_higher_than.min.agevent.grade.0(df, 2)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 3)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 4)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 5)

# counts for 0 to X transitions with min ageevent 0
check_grades_transition.agevent.grade.0(df, 0, 2)
check_grades_transition.agevent.grade.0(df, 0, 3)
check_grades_transition.agevent.grade.0(df, 0, 4)
check_grades_transition.agevent.grade.0(df, 0, 5)

## Processing...
# Condition 1: Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher. If there are two ageevent that are miniumum value, we still apply this filter (but this did not make any difference)
step1 <- df %>%
  group_by(sjlid) %>%
  filter(
    !(grade[which(ageevent == min(ageevent))[1]] != 0)
  ) %>%
  ungroup()

dim(step1)
# anti_join(df, step1)

# Condition 2: Proceed to remove rows where grades are smaller than the grades previously seen in the ordered rows. In simpler terms, we keeps the rows where the grade value is the highest seen so far within each group. This effectively retains only the rows with the highest grade value within each sjlid group.
step2 <- step1 %>% 
  arrange(sjlid, ageevent) %>% 
  group_by(sjlid) %>% 
  filter(cummax(grade) == grade) %>% 
  ungroup()

dim(step2)
step2$tb_number <- paste0(step2$sjlid, ":", step2$Sample_age)
# Condition 3: Keep rows with max num_vials and vitalstatus
step3 <- step2 %>%
  group_by(tb_number) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(vitalstatus == "Alive" | all(vitalstatus == "Deceased"))

dim(step3)
step3$tb_number[duplicated(step3$tb_number)] # check duplicates
# removed_rows <- anti_join(step2, step3)

# # # Remove rows with duplicate grade, ageevent, sample_age combinations (it only removes 23 rows; need to discuss with Yadav!)
step3 <- step3 %>%
  distinct(sjlid, grade, ageevent, Sample_age, .keep_all = TRUE)

dim(step3) # 5992   10
# removed_rows <- anti_join(step2, step3)

FINAL.1 <- step3
FINAL.PLASMA <- step3

# ## Get counts in non-missing events
check_grades_eq_or_higher_than(FINAL.1, 2)
check_grades_eq_or_higher_than(FINAL.1, 3)
check_grades_eq_or_higher_than(FINAL.1, 4)
check_grades_eq_or_higher_than(FINAL.1, 5)

check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 2)
check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 3)
check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 4)
check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 5)

# counts for 0 to X transitions with any grade at min ageevent
check_grades_transition(FINAL.1, 0, 2)
check_grades_transition(FINAL.1, 0, 3)
check_grades_transition(FINAL.1, 0, 4)
check_grades_transition(FINAL.1, 0, 5)

# counts for 0 to X transitions with min ageevent 0
check_grades_transition.agevent.grade.0(FINAL.1, 0, 2)
check_grades_transition.agevent.grade.0(FINAL.1, 0, 3)
check_grades_transition.agevent.grade.0(FINAL.1, 0, 4)
check_grades_transition.agevent.grade.0(FINAL.1, 0, 5)

##################################
table(FINAL.PLASMA$grade)
table(FINAL.SERUM$grade)

dim(FINAL.PLASMA)
# [1] 6039    9
# dim(FINAL.SERUM)
# [1] 5992    9

sum(unique(FINAL.SERUM$sjlid) %in% unique(FINAL.PLASMA$sjlid))
# 4066
sum(!unique(FINAL.SERUM$sjlid) %in% unique(FINAL.PLASMA$sjlid))
# 39

serum.grade.2.or.higher <- FINAL.SERUM[FINAL.SERUM$grade > 0,]
plasma.grade.2.or.higher <- FINAL.PLASMA[FINAL.PLASMA$grade > 0,]

in.both <- plasma.grade.2.or.higher$sjlid[plasma.grade.2.or.higher$sjlid %in% serum.grade.2.or.higher$sjlid]
only.in.plasma <- plasma.grade.2.or.higher$sjlid[!plasma.grade.2.or.higher$sjlid %in% serum.grade.2.or.higher$sjlid]
# "SJL1066301" "SJL1277801" "SJL1293701" "SJL1795407" "SJL2542201" "SJL4765507" "SJL4828001" "SJL5087912" "SJL5383017"
only.in.serum <- serum.grade.2.or.higher$sjlid[!serum.grade.2.or.higher$sjlid %in% plasma.grade.2.or.higher$sjlid]
# "SJL5007201" "SJL5032107"


## Add event number:
FINAL.PLASMA <- FINAL.PLASMA %>%
  group_by(sjlid) %>%
  arrange(ageevent) %>%
  mutate(event_number = row_number()) %>%
  ungroup() %>%
  arrange(sjlid)  # Restore the original order

## Add event number:
FINAL.SERUM <- FINAL.SERUM %>%
  group_by(sjlid) %>%
  arrange(ageevent) %>%
  mutate(event_number = row_number()) %>%
  ungroup() %>%
  arrange(sjlid)  # Restore the original order


table(FINAL.PLASMA$event_number, FINAL.PLASMA$grade)
table(FINAL.SERUM$event_number, FINAL.SERUM$grade)

FINAL.PLASMA <- FINAL.PLASMA[c("sjlid", "aliquot_type", "num_vials", "vitalstatus", "grade", "Sample_age", "ageevent", "event_number")]
FINAL.SERUM <- FINAL.SERUM[c("sjlid", "aliquot_type", "num_vials", "vitalstatus", "grade", "Sample_age", "ageevent", "event_number")]
write.table(FINAL.PLASMA, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/Plasma_data_processed_v4.txt", col.names = T, row.names = F, quote = F)
write.table(FINAL.SERUM, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/Serum_data_processed_v4.txt", col.names = T, row.names = F, quote = F)


## Yadav: do you know how many of the 3853 had data at multiple timepoints and what is the breakdown by 2 times, 3 times, etc
df.3853$Sample_age
# cc <- df.3853[grepl("SJL0253301|SJL1063201|SJL1063401|SJL1063501|SJL1063101|TB-12-5258|TB-14-3725", df.3853$sjlid), c("sjlid", "tb_number", "num_vials", "aliquot_type", "volume", "Sample_age", "Survival_Status", "grade", "ageevent")]
# dput(cc)



#####################################
## Sample age grouped by tb_number ##
#####################################
## Sample age
timepoint_counts <- df.3498 %>%
  group_by(tb_number) %>%
  summarize(unique_timepoints = n_distinct(Sample_age)) %>%
  ungroup()

# Count the occurrences of each unique_timepoints count
breakdown <- timepoint_counts %>%
  count(unique_timepoints)

View(breakdown)


## agevent
timepoint_counts <- df.3498 %>%
  group_by(tb_number) %>%
  summarize(unique_timepoints = n_distinct(ageevent)) %>%
  ungroup()

# Count the occurrences of each unique_timepoints count
breakdown <- timepoint_counts %>%
  count(unique_timepoints)

View(breakdown)


# Group data by sjlid and grade, then count occurrences
grade_counts <- df.3498 %>%
  group_by(tb_number, grade) %>%
  tally() %>%
  ungroup()

# Count the occurrences of each grade count
breakdown <- data.frame(table(grade_counts$grade))
View(breakdown)


#################################
## Sample age grouped by SJLID ##
#################################
timepoint_counts <- df.3498 %>%
  group_by(sjlid) %>%
  summarize(unique_timepoints = n_distinct(Sample_age)) %>%
  ungroup()

# Count the occurrences of each unique_timepoints count
breakdown <- timepoint_counts %>%
  count(unique_timepoints)

View(breakdown)


## agevent
timepoint_counts <- df.3498 %>%
  group_by(sjlid) %>%
  summarize(unique_timepoints = n_distinct(ageevent)) %>%
  ungroup()

# Count the occurrences of each unique_timepoints count
breakdown <- timepoint_counts %>%
  count(unique_timepoints)

View(breakdown)


# Group data by sjlid and grade, then count occurrences
grade_counts <- df.3498 %>%
  group_by(sjlid, grade) %>%
  tally() %>%
  ungroup()

# Count the occurrences of each grade count
breakdown <- data.frame(table(grade_counts$grade))
View(breakdown)
