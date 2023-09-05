## Conditions::
# 1. within each sample if all rows have grade 2 or higher, or if minimum ageevent has grade 2 or higher, skip this sample. Do not extract any rows from that sample.
# 2. Then, within each tb_number, keep rows with max num_vials; then keep ALIVE over DEAD (something like this: filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD")))
# 3. Then, within each sample, if all grades are zero for all tb_number, keep rows with unique ageevent per sample and also with unique ageevent (making sure all controls are included). If, there are grades 2 or higher (including grade 0) in any tb_number, then continue with step number 4.
# 4. Within each sample, if the minimum age event has has zero grade, then look for 0 or higher grades (grade 0, 2,3,4,5, etc). Once that done, among the extracted rows (within each sample), make sure all grade 0 has different ageevent; also make sure any grade 2, has ageevent greater than any grade 0 ageevent; any grade 3 has ageevent greater than any ageevent of grade 2 or lower grade; any grade 4 has ageevent greater than any ageevent of grade 3 or lower grades; any grade 5 has ageevent greater than any ageevent of grade 4 or lower grades; and so forth; Thus, make sure all grades have different ageevent. 
# 5. Then output the dataframe with the extracted rows

setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")
library(dplyr)





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

# Round age down to one decimal place, so easier to compare
all.df$Sample_age <- floor(all.df$Sample_age * 10) / 10
all.df$ageevent <- floor(all.df$ageevent * 10) / 10


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
check_grades_eq_or_higher_than(df, 6)

# counts for 0 to X transitions with any grade at min ageevent
check_grades_transition(df, 0, 2)
check_grades_transition(df, 0, 3)
check_grades_transition(df, 0, 4)
check_grades_transition(df, 0, 5)
check_grades_transition(df, 0, 6)

# counts for 0 grade at min ageevent
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 2)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 3)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 4)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 5)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 6)

# counts for 0 to X transitions with min ageevent 0
check_grades_transition.agevent.grade.0(df, 0, 2)
check_grades_transition.agevent.grade.0(df, 0, 3)
check_grades_transition.agevent.grade.0(df, 0, 4)
check_grades_transition.agevent.grade.0(df, 0, 5)
check_grades_transition.agevent.grade.0(df, 0, 6)


## Next, filter rows where ageevent is greater than sample age
df <- df[!is.na(df$Sample_age) & df$ageevent >= df$Sample_age, ]
df.3724 <- df

## Get counts in non-missing events
check_grades_eq_or_higher_than(df, 2)
check_grades_eq_or_higher_than(df, 3)
check_grades_eq_or_higher_than(df, 4)
check_grades_eq_or_higher_than(df, 5)
check_grades_eq_or_higher_than(df, 6)

# counts for 0 to X transitions with any grade at min ageevent
check_grades_transition(df, 0, 2)
check_grades_transition(df, 0, 3)
check_grades_transition(df, 0, 4)
check_grades_transition(df, 0, 5)
check_grades_transition(df, 0, 6)

check_grades_eq_or_higher_than.min.agevent.grade.0(df, 2)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 3)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 4)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 5)
check_grades_eq_or_higher_than.min.agevent.grade.0(df, 6)

# counts for 0 to X transitions with min ageevent 0
check_grades_transition.agevent.grade.0(df, 0, 2)
check_grades_transition.agevent.grade.0(df, 0, 3)
check_grades_transition.agevent.grade.0(df, 0, 4)
check_grades_transition.agevent.grade.0(df, 0, 5)
check_grades_transition.agevent.grade.0(df, 0, 6)

## Processing...

# Condition 1: Skip samples with grade 2 or higher or minimum ageevent with grade 2 or higher
step1 <- df %>%
  group_by(sjlid) %>%
  filter(
    !(grade[which.min(ageevent)] != 0)
  ) %>%
  ungroup()

dim(step1)

# Condition 2: Keep rows with max num_vials and Survival_Status
step2 <- step1 %>%
  group_by(tb_number) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD"))

dim(step2)

# Step 3: Keep rows with unique ageevent if all grades are zero
# If there are grades 2 or higher (including grade 0), continue to Step 4
step3 <- step2 %>%
  group_by(sjlid, ageevent) %>%
  mutate(
    all_zero = all(grade == 0),
    any_non_zero = any(grade != 0)
  ) %>%
  filter(!all_zero | (all_zero & !any_non_zero & row_number() == 1)) %>%
  ungroup() %>%
  select(-all_zero, -any_non_zero)

dim(step3)

## remove duplicate events if all grades are zero
step3 <- step3 %>%
  group_by(sjlid) %>%
  mutate(
    has_non_zero = any(grade != 0)
  ) %>%
  filter(!has_non_zero | grade == 0 | !duplicated(ageevent)) %>%
  ungroup() %>%
  select(-has_non_zero) %>%
  group_by(sjlid, ageevent) %>%
  distinct(.keep_all = TRUE) %>%
  ungroup()

dim(step3)

# Step 4: Continue with further processing
transformed_df <- step3 %>%
  group_by(sjlid) %>%
  mutate(
    has_grade_2_or_higher = any(grade >= 2),
    min_ageevent_zero = min(ageevent[grade == 0])
  ) %>%
  filter(
    !has_grade_2_or_higher |
      (grade > 0 & ageevent > min_ageevent_zero) |
      (grade == 0 & ageevent == min_ageevent_zero)
  ) %>%
  ungroup() %>%
  select(-has_grade_2_or_higher, -min_ageevent_zero) %>%
  arrange(sjlid, tb_number, grade, ageevent) %>%
  group_by(sjlid, tb_number, grade) %>%
  mutate(
    rank_age = dense_rank(ageevent)
  ) %>%
  ungroup()

dim(transformed_df)


# Remove rows with duplicate grade and ageevent combinations
transformed_df <- transformed_df %>%
  distinct(sjlid, grade, ageevent, .keep_all = TRUE)

dim(transformed_df)
FINAL.1 <- transformed_df


# ## Get counts in non-missing events
check_grades_eq_or_higher_than(FINAL.1, 2)
check_grades_eq_or_higher_than(FINAL.1, 3)
check_grades_eq_or_higher_than(FINAL.1, 4)
check_grades_eq_or_higher_than(FINAL.1, 5)
check_grades_eq_or_higher_than(FINAL.1, 6)

# counts for 0 to X transitions with any grade at min ageevent
check_grades_transition(FINAL.1, 0, 2)
check_grades_transition(FINAL.1, 0, 3)
check_grades_transition(FINAL.1, 0, 4)
check_grades_transition(FINAL.1, 0, 5)
check_grades_transition(FINAL.1, 0, 6)

check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 2)
check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 3)
check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 4)
check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 5)
check_grades_eq_or_higher_than.min.agevent.grade.0(FINAL.1, 6)


# counts for 0 to X transitions with min ageevent 0
check_grades_transition.agevent.grade.0(FINAL.1, 0, 2)
check_grades_transition.agevent.grade.0(FINAL.1, 0, 3)
check_grades_transition.agevent.grade.0(FINAL.1, 0, 4)
check_grades_transition.agevent.grade.0(FINAL.1, 0, 5)
check_grades_transition.agevent.grade.0(FINAL.1, 0, 6)


##################################################
## 2. Now, work on those with missing age event ##
##################################################
# if ageevent is missing and grade is zero, keep (eg. SJL5237316, SJL5257302)
df <- all.df[all.df$sjlid %in% missing.age.samples,]
sum(is.na(df$ageevent)) # 323
sum(is.na(df$ageevent) & df$grade == 0) # 323  ## all are grade 0 and with missing ageevent, so we can run the same code

## Processing...

# Condition 1: Skip samples with grade 2 or higher or minimum sample age with grade 2 or higher
step1 <- df %>%
  group_by(sjlid) %>%
  filter(
    !(grade[which.min(Sample_age)] != 0)
  ) %>%
  ungroup()

dim(step1)

# Condition 2: Keep rows with max num_vials and Survival_Status
step2 <- step1 %>%
  group_by(tb_number) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD"))

dim(step2)

# Step 3: Keep rows with unique ageevent if all grades are zero
# If there are grades 2 or higher (including grade 0), continue to Step 4
step3 <- step2 %>%
  group_by(sjlid, ageevent) %>%
  mutate(
    all_zero = all(grade == 0),
    any_non_zero = any(grade != 0)
  ) %>%
  filter(!all_zero | (all_zero & !any_non_zero & row_number() == 1)) %>%
  ungroup() %>%
  select(-all_zero, -any_non_zero)

dim(step3)

## remove duplicate events if all grades are zero
step3 <- step3 %>%
  group_by(sjlid) %>%
  mutate(
    has_non_zero = any(grade != 0)
  ) %>%
  filter(!has_non_zero | grade == 0 | !duplicated(ageevent)) %>%
  ungroup() %>%
  select(-has_non_zero) %>%
  group_by(sjlid, ageevent) %>%
  distinct(.keep_all = TRUE) %>%
  ungroup()

dim(step3)

# Step 4: Continue with further processing
transformed_df <- step3 %>%
  group_by(sjlid) %>%
  mutate(
    has_grade_2_or_higher = any(grade >= 2),
    min_ageevent_zero = min(ageevent[grade == 0])
  ) %>%
  filter(
    !has_grade_2_or_higher |
      (grade > 0 & ageevent > min_ageevent_zero) |
      (grade == 0 & ageevent == min_ageevent_zero)
  ) %>%
  ungroup() %>%
  select(-has_grade_2_or_higher, -min_ageevent_zero) %>%
  arrange(sjlid, tb_number, grade, ageevent) %>%
  group_by(sjlid, tb_number, grade) %>%
  mutate(
    rank_age = dense_rank(ageevent)
  ) %>%
  ungroup()

dim(transformed_df)


# Remove rows with duplicate grade and ageevent combinations
transformed_df <- transformed_df %>%
  distinct(sjlid, grade, ageevent, .keep_all = TRUE)

dim(transformed_df)
FINAL.2 <- transformed_df


## Merge 1 and 2.
FINAL <- rbind.data.frame(FINAL.1, FINAL.2) 


write.table(FINAL, file = "serum_data_processed_final.txt", sep = "\t",  row.names = FALSE, col.names = TRUE, quote = F)
df.3447 <- FINAL



samples_with_grade_2 <- transformed_df %>%
  filter(grade == 2|grade == 3) %>%
  distinct(sjlid) %>%
  nrow()

samples_with_grade_2_or_higher <- transformed_df %>%
  filter(grade >= 2) %>%
  distinct(sjlid) %>%
  nrow()

cat("Number of samples with grade 2:", samples_with_grade_2, "\n")



##################################
## Yadav: do you know how many of the 3853 had data at multiple timepoints and what is the breakdown by 2 times, 3 times, etc
df.3853$Sample_age
cc <- df.3853[grepl("SJL0253301|SJL1063201|SJL1063401|SJL1063501|SJL1063101|TB-12-5258|TB-14-3725", df.3853$sjlid), c("sjlid", "tb_number", "num_vials", "aliquot_type", "volume", "Sample_age", "Survival_Status", "grade", "ageevent")]
dput(cc)


## Sample age
timepoint_counts <- df.3447 %>%
  group_by(sjlid) %>%
  summarize(unique_timepoints = n_distinct(Sample_age)) %>%
  ungroup()

# Count the occurrences of each unique_timepoints count
breakdown <- timepoint_counts %>%
  count(unique_timepoints)

View(breakdown)


## agevent
timepoint_counts <- df.3447 %>%
  group_by(sjlid) %>%
  summarize(unique_timepoints = n_distinct(ageevent)) %>%
  ungroup()

# Count the occurrences of each unique_timepoints count
breakdown <- timepoint_counts %>%
  count(unique_timepoints)

View(breakdown)


# Group data by sjlid and grade, then count occurrences
grade_counts <- df.3853 %>%
  group_by(sjlid, grade) %>%
  tally() %>%
  ungroup()

# Count the occurrences of each grade count
breakdown <- data.frame(table(grade_counts$grade))
View(breakdown)
