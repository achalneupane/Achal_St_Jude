library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")

source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/get_matching_rows_from_CTCAE.R")
df <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
dim(df) # 20174

## add TB from Mathew
TB <- read.table("InventoryforAchal_from_Mathew.txt", header = T, sep = "\t")
TB$KEY <- paste0(TB$sjlid, TB$num_vials, TB$aliquot_type, TB$ageatsample)

df$Key <- paste0(df$sjlid, df$num_vials, df$aliquot_type, df$ageatsample)
sum(df$Key %in% TB$KEY)
df$tb_number <- TB$tb_number[match(df$Key, TB$KEY)]


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
# SERUM <- df[grepl("Serum", df$aliquot_type, ignore.case = T),]
SERUM <- df[grepl("Plasma", df$aliquot_type, ignore.case = T),]

SERUM.original <- SERUM

## I see age at serum sample is <18 yrs. Could you identify the samples among 18
#or higher only? Everyone needs to be adults at serum sample.
SERUM <- SERUM[which(SERUM$ageatsample >= 18),]

# remove vial zero
SERUM <- SERUM[SERUM$num_vials > 0 ,]

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
# 7221 # Serum, 18 or older, vial > 0
# 7367 Plasma, 18 or older, vial > 0

# CTCAE <- get_rows_with_smaller_sample_age(CTCAE, SERUM, 0)
CTCAE <- get_rows_with_smaller_sample_age.all(CTCAE, SERUM, 7)
dim(CTCAE)
# 8332
# 7221 SERUM with > 0 vial
# 7367 PLASMA with > 0 vial

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
# 0    2    3    5 
# 5120  209   56    1 

## removing num vial 0
# 0    2    3    5 
# 4388  207   54    1 
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
# 0    2    3    5
# 1 3563    0    0    0
# 2 1495   69   25    0
# 3  571   50   11    0
# 4  164   19    8    1
# 5   31    5    2    0
# 6   11    0    1    0
# 7    5    1    1    0
# 8    1    0    0    0

# 0    2    3    5
# 1 3056    0    0    0
# 2 1342   69   24    0
# 3  520   48    9    0
# 4  156   18    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

# removing num vial zero
# 0    2    3    5
# 1 2352    0    0    0
# 2 1320   44   14    0
# 3  514   46    8    0
# 4  156   18    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

# PLASMA with > 0 vial
# 0    2    3    5
# 1 3029    0    0    0
# 2 1337   73   22    0
# 3  516   49    9    0
# 4  155   19    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0
table(CTCAE.2$event_number,CTCAE.2$grade_2_or_higher)
# 1    3563                 0
# 2    1495                94
# 3     571                61
# 4     164                28
# 5      31                 7
# 6      11                 1
# 7       5                 2
# 8       1                 0

# grade_0 grade_2_or_higher  # after removing younger serum < 18 yo
# 1    3056                 0
# 2    1342                93
# 3     520                57
# 4     156                26
# 5      29                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0

# Serum removing num vial zero
# grade_0 grade_2_or_higher   # after removing num vial= zero
# 1    2352                 0
# 2    1320                58
# 3     514                54
# 4     156                26
# 5      29                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0

# Plasma removing num vial zero
# grade_0 grade_2_or_higher
# 1    3029                 0
# 2    1337                95
# 3     516                58
# 4     155                27
# 5      29                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0
table(CTCAE.2$new_event_number,CTCAE.2$grade)
# 0    2    3    5
# 1 4085    0    0    0
# 2 1296   89   34    0
# 3  383   43   13    0
# 4   74   11    1    1
# 5    3    1    0    0

# after removing younger serum < 18 yo
# 0    2    3    5
# 1 3402    0    0    0
# 2 1258   86   28    0
# 3  383   42   13    0
# 4   74   11    1    1
# 5    3    1    0    0

# SERUM after removing num vial zero
# 0    2    3    5  
# 1 3174    0    0    0
# 2  960   82   24    0
# 3  226   23    7    0
# 4   27    8    0    1
# 5    1    0    0    0

# Plasma removing num vial zero
# 0    2    3    5
# 1 3384    0    0    0
# 2 1258   93   28    0
# 3  367   42   10    1
# 4   71   10    2    0
# 5    3    1    0    0
table(CTCAE.2$new_event_number,CTCAE.2$grade_2_or_higher)
# 1    4085                 0
# 2    1296               123
# 3     383                56
# 4      74                13
# 5       3                 1

# grade_0 grade_2_or_higher # after removing younger serum < 18 yo
# 1    3402                 0
# 2    1258               114
# 3     383                55
# 4      74                13
# 5       3                 1

# SERUM after removing num vial = 0
# grade_0 grade_2_or_higher
# 1    3174                 0
# 2     960               106
# 3     226                30
# 4      27                 9
# 5       1                 0

# PLASMA after removing num vial = 0
# grade_0 grade_2_or_higher
# 1    3384                 0
# 2    1258               121
# 3     367                53
# 4      71                12
# 5       3                 1

CTCAE.3 <- CTCAE.2[!is.na(CTCAE.2$Sample_age),]
table(CTCAE.3$event_number,CTCAE.3$grade_2_or_higher)
# 1    3563                 0
# 2    1495                60
# 3     571                30
# 4     164                13
# 5      31                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0

# after removing younger serum < 18 yo
# grade_0 grade_2_or_higher
# 1    3056                 0
# 2    1342                59
# 3     520                28
# 4     156                12
# 5      29                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0

# After removing num vial = zero
# grade_0 grade_2_or_higher
# 1    2352                 0
# 2    1320                32
# 3     514                26
# 4     156                12
# 5      29                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0

table(is.na(CTCAE.2$Sample_age) & CTCAE.2$grade >=2)
table(CTCAE.2$grade >=2)
# FALSE  TRUE 
# 5841   193 

# FALSE  TRUE 
# 5120   183 

# FALSE  TRUE 
# 4460    73

# after removing vial zero
# FALSE  TRUE 
# 4388   145 
# write.table(CTCAE.2, "Serum_data_processed_v8.txt", col.names = T, row.names = F, sep = "\t")
# write.table(CTCAE.2, "Plasma_data_processed_v8_after_removing_numvial_0.txt", col.names = T, row.names = F, sep = "\t") ## based on Plasma data
write.table(CTCAE.2, "Serum_data_processed_v8_after_removing_numvial_0.txt", col.names = T, row.names = F, sep = "\t") ## based on Serum data

