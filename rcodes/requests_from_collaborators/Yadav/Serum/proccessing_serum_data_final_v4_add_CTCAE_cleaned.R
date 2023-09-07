library(dplyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")

## Could you please work on it by merging with the CTCAE grades for cardiomyopathy from the most recent data freeze? Please look at serum and plasma separately.
# read CTCAE
# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
CTCAE <- CTCAE[grepl("Cardiomyopathy", CTCAE$condition),]
CTCAE.original <- CTCAE
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


## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
PLASMA <- df[grepl("Plasma", df$aliquot_type, ignore.case = T),]
# PLASMA.within.7.days.age.eevnt <- get_matching_rows(PLASMA, CTCAE, 7/365.25) # within 7 days
# PLASMA.within.7.days.age.eevnt <- PLASMA.within.7.days.age.eevnt[!is.na(PLASMA.within.7.days.age.eevnt$grade),]
PLASMA <- get_matching_rows(PLASMA, CTCAE, 0) # on same day
PLASMA <- PLASMA[!is.na(PLASMA$grade) & PLASMA$grade != -9,]

# PLASMA <- PLASMA %>%
#   group_by(sjlid) %>%
#   arrange(desc(ageatsample), desc(vitalstatus), desc(num_vials)) %>%
#   slice(1) %>%
#   ungroup()

## Note: Since we are rounding age down to one decimal place, a 7-days window did not make any difference in terms of the number of rows or samples
SERUM <- df[grepl("Serum", df$aliquot_type, ignore.case = T),]
# # SERUM <- SERUM[grepl("SJL1527107|SJL2100910|SJL2279209", SERUM$sjlid),]
# SERUM.within.7.days.age.eevnt <- get_matching_rows(SERUM, CTCAE, 7/365.25) # within 7 days
# SERUM.within.7.days.age.eevnt <- SERUM.within.7.days.age.eevnt[!is.na(SERUM.within.7.days.age.eevnt$grade),]
# SERUM.within.0.days.age.eevnt <- get_matching_rows(SERUM, CTCAE, 0) # within same days
# SERUM.within.0.days.age.eevnt <- SERUM.within.0.days.age.eevnt[!is.na(SERUM.within.0.days.age.eevnt$grade),]

SERUM <- get_matching_rows(SERUM, CTCAE, 0) # on same day
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


# ## Get counts in non-missing events
# # counts for any grade at min ageevent
# check_grades_eq_or_higher_than(df, 2)
# check_grades_eq_or_higher_than(df, 3)
# check_grades_eq_or_higher_than(df, 4)
# check_grades_eq_or_higher_than(df, 5)
# 
# # counts for 0 to X transitions with any grade at min ageevent
# check_grades_transition(df, 0, 2)
# check_grades_transition(df, 0, 3)
# check_grades_transition(df, 0, 4)
# check_grades_transition(df, 0, 5)
# 
# # counts for 0 grade at min ageevent
# check_grades_eq_or_higher_than.min.agevent.grade.0(df, 2)
# check_grades_eq_or_higher_than.min.agevent.grade.0(df, 3)
# check_grades_eq_or_higher_than.min.agevent.grade.0(df, 4)
# check_grades_eq_or_higher_than.min.agevent.grade.0(df, 5)
# 
# # counts for 0 to X transitions with min ageevent 0
# check_grades_transition.agevent.grade.0(df, 0, 2)
# check_grades_transition.agevent.grade.0(df, 0, 3)
# check_grades_transition.agevent.grade.0(df, 0, 4)
# check_grades_transition.agevent.grade.0(df, 0, 5)

df <- SERUM
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
