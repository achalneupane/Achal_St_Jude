## Conditions::
# 1. Remove rows where ageevent is not greater than (or within 1 week of) sample age
# 2. within each sample if all rows have grade 2 or higher, or if minimum ageevent has grade 2 or higher, skip this sample. Do not extract any rows from that sample.
# 3. Withing each sample, remove rows where grades are smaller than the grades previously seen in the ordered rows.
# 4. Then, within each tb_number, keep rows with max num_vials; then keep ALIVE over DEAD (something like this: filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD")))

setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")
library(dplyr)
library(haven)

##############################
# ## Count by SJLID samples
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/Count_functions_by_sjlid.R")

## Count by serum samples
# source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/Count_functions_by_tb_number.R")

all.df <- read.delim("All_serum_samples_for_R01_22Aug2023_for_Achal_edited.txt", sep = "\t", header = T, check.names = F)
dim(all.df)

head(all.df)

sum(is.na(all.df$tb_number))

all.df$original_Sample_age <- all.df$Sample_age
all.df$original_ageevent <- all.df$ageevent

all.df$Sample_age <- as.numeric(all.df$Sample_age)

# # Round age down to one decimal place, so easier to compare
# all.df$Sample_age <- floor(all.df$Sample_age * 10) / 10
# all.df$ageevent <- floor(all.df$ageevent * 10) / 10


## read echo data
echo <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/echo_machine.sas7bdat')
cc <- echo[c(1:10,73)]
tracking <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data/tracking.sas7bdat')

##################################################
## 1. First process df without ageevent missing ##
##################################################
# Remove rows with missing ageevent and process them separately
missing.age.samples <- all.df$sjlid[is.na(all.df$ageevent)]
## Process them separately
df <- all.df[!all.df$sjlid %in% missing.age.samples,]
df.3853 <- df
## Get counts in non-missing events
# counts for any grade at min ageevent
check_grades_eq_or_higher_than(df, 2)
check_grades_eq_or_higher_than(df, 3)
check_grades_eq_or_higher_than(df, 4)
check_grades_eq_or_higher_than(df, 5)

# counts for 0 to X transitions with any grade at min ageevent
check_grades_transition(df, 0, 2)
check_grades_transition(df, 0, 3)
check_grades_transition(df, 0, 4)
check_grades_transition(df, 0, 5)

# counts for 0 grade at min ageevent
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 2)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 3)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 4)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 5)

# counts for 0 to X transitions with min ageevent 0
check_grades_transition.agevent.grade.0(df, 0, 2)
check_grades_transition.agevent.grade.0(df, 0, 3)
check_grades_transition.agevent.grade.0(df, 0, 4)
check_grades_transition.agevent.grade.0(df, 0, 5)


## Next, filter rows where ageevent is greater than (or within 1 week of) sample age
# df <- df[!is.na(df$Sample_age) & df$ageevent >= df$Sample_age, ]
df <- df[!is.na(df$Sample_age) & df$ageevent+7/365.25 >= df$Sample_age, ]
df.3782 <- df

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

# Condition 3: Keep rows with max num_vials and Survival_Status
step3 <- step2 %>%
  group_by(tb_number) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD"))

dim(step3)
# removed_rows <- anti_join(step2, step3)

# # # Remove rows with duplicate grade, ageevent, sample_age combinations (it only removes 23 rows; need to discuss with Yadav!)
# step3 <- step3 %>%
#   distinct(sjlid, grade, ageevent, Sample_age, .keep_all = TRUE)

dim(step3)
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


##################################################
## 2. Now, work on those with missing age event ##
##################################################
# if ageevent is missing and grade is zero, keep (eg. SJL5237316, SJL5257302)
df <- all.df[all.df$sjlid %in% missing.age.samples,]
sum(is.na(df$ageevent)) # 323
sum(is.na(df$ageevent) & df$grade == 0) # 323  ## all are grade 0 and with missing ageevent, so we can run the same code

## No need to process anything further
# ## Remove rows with duplicate grade, ageevent, sample_age combinations (it only removes 23 rows; need to discuss with Yadav!)
# df <- df %>%
#   distinct(sjlid, grade, ageevent, Sample_age, .keep_all = TRUE)

dim(df)
FINAL.2 <- df


## Merge 1 and 2.
FINAL <- rbind.data.frame(FINAL.1, FINAL.2) 
dim(FINAL)

write.table(FINAL, file = "serum_data_processed_final_v3_cleaned.txt", sep = "\t",  row.names = FALSE, col.names = TRUE, quote = F)
df.3498 <- FINAL.1




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
