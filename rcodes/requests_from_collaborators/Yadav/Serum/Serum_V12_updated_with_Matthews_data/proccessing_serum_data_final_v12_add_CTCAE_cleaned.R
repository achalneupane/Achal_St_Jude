## This version is esentially the same as V11. Only update is made to proccessing_serum_data_final_v11_add_CTCAE_cleaned.R file where I removed TB with last vials, with vial 1 when subsequent visits have more number of vials.

library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/")

source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/get_matching_rows_from_CTCAE.R")
# df <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
# dim(df) # 20174

## add TB from Mathew
# TB <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/ForAchal_Survivors.txt", header = T, sep = "\t") ## March version
TB <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/Achal_survivors_04.04.2024.txt", header = T, sep = "\t") ## Updated by Matt in April

df <- TB

df$ageatsample <- floor(df$ageatsample * 10) / 10
df.original <- df
df <- df %>%
  distinct()
dim(df)  # 20137
# 20779 # Mathew ## March
# 21041 ## April version
# ## Remove those with vials less than 2
# morethan.2.vials <- df %>%
#   dplyr::group_by(sjlid, ageatsample) %>%
#   dplyr::mutate(total_num_vials = sum(num_vials)) %>%
#   ungroup()
# morethan.2.vials <- morethan.2.vials[morethan.2.vials$total_num_vials >= 2,]
# dim(morethan.2.vials)
# # [1] 20104     8
# morethan.2.vials <- morethan.2.vials[morethan.2.vials$vitalstatus == "Alive",]
# dim(morethan.2.vials)
# # [1] 19161     8
# morethan.2.vials <- morethan.2.vials[grepl("Plasma", morethan.2.vials$aliquot_type),]
# morethan.2.vials <- morethan.2.vials[morethan.2.vials$num_vials > 0,]

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
dim(SERUM.original)
# [1] 10383     6
# 10514  7 ## April version

## I see age at serum sample is <18 yrs. Could you identify the samples among 18
#or higher only? Everyone needs to be adults at serum sample.
SERUM <- SERUM[which(SERUM$ageatsample >= 18),]

# remove vial zero
SERUM <- SERUM[SERUM$num_vials > 0 ,]


## 171 cases to keep
all.wanted.df.1200.to.update <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v11_output//plasma_data_batch1_1200_samples.txt", header = T, sep = "\t")
keep.171.cases <- all.wanted.df.1200.to.update$sjlid[all.wanted.df.1200.to.update$selection_group == "171_CMP_cases"]
keep.171.cases <- SERUM[SERUM$sjlid %in% keep.171.cases,]
SERUM <- SERUM[!SERUM$sjlid %in% keep.171.cases$sjlid,] # exclude 171 cases
dim(SERUM)
# 8536    7
## Now keep only those that have more than one vial and are alive
SERUM <- SERUM[SERUM$num_vials > 1,]
SERUM <- rbind.data.frame(SERUM, keep.171.cases)
dim(SERUM)
# [1] 7948    6


dim(SERUM)
# 8879
# 7948    7 ## April 4 version

## Keep max num vials and alive over dead
SERUM <- SERUM %>%
  group_by(sjlid, ageatsample) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(vitalstatus == "Alive" | all(vitalstatus == "Deceased"))

dim(SERUM)
# [1] 7934    6

## Also removing samples that have duplicate ageevent in CTCAE data
# SERUM <- SERUM[!SERUM$sjlid %in% c("SJL1225801", "SJL1265801", "SJL1430801", "SJL4730101", "SJL5134305", "SJL5146506"),]
dim(SERUM)
# 7934    6
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
# 9218    9

# V10 & V11
# 9218    9

CTCAE <- CTCAE[CTCAE$sjlid %in% unique(SERUM$sjlid),]
dim(CTCAE)
# 8332
# 7221 # Serum, 18 or older, vial > 0
# 7367 Plasma, 18 or older, vial > 0

## V9 (after updating with Mathews data)
# 7428    9 # plasma
# 7286    9 # serum

# V10
# 6826    9

# V11
# 7406    9

# v12 7148

# CTCAE <- get_rows_with_smaller_sample_age(CTCAE, SERUM, 0)
CTCAE <- get_rows_with_smaller_sample_age.all(CTCAE, SERUM, 7)
dim(CTCAE)
# 8332
# 7221 SERUM with > 0 vial
# 7367 PLASMA with > 0 vial

## V9 (after updating with Mathews data)
# 7428   12 # plasma
# # 7286   12 # serum

# V10
# 6826   13

#V11
# 7406   13

# V12
dim(CTCAE)
# 7148   13

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


## V9 (after updating with Mathews data)
# Plasma
# 0    2    3    5 
# 5085  211   53    1 

# serum
# 0    2    3    5 
# 4390  209   55    1

## V10
# 0    2    3    5 
# 4413  194   52    1 

## V11
# 0    2    3    5 
# 5070  211   53    1 

## V12
# 0    2    3    5 
# 4540  200   53    1 

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

## V9 (after updating with Mathews data)
# 0    2    3    5 ## plasma
# 1 3030    0    0    0
# 2 1338   73   22    0
# 3  516   49    9    0
# 4  155   19    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

# 0    2    3    5 ## serum
# 1 2353    0    0    0
# 2 1321   44   14    0
# 3  514   46    8    0
# 4  156   18    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

## V10
# 0    2    3    5
# 1 2499    0    0    0
# 2 1238   72   22    0
# 3  485   49    9    0
# 4  145   19    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

## V11
# 0    2    3    5
# 1 3025    0    0    0
# 2 1332   73   22    0
# 3  515   49    9    0
# 4  153   19    7    1
# 5   28    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

## V12
# 0    2    3    5
# 1 2549    0    0    0
# 2 1286   72   22    0
# 3  508   49    9    0
# 4  151   19    7    1
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

## V9 (after updating with Mathews data)
# grade_0 grade_2_or_higher ## Plasma
# 1    3030                 0
# 2    1338                95
# 3     516                58
# 4     155                27
# 5      29                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0

# grade_0 grade_2_or_higher ## Serum
# 1    2353                 0
# 2    1321                58
# 3     514                54
# 4     156                26
# 5      29                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0

## V10
# grade_0 grade_2_or_higher
# 1    2499                 0
# 2    1238                94
# 3     485                58
# 4     145                27
# 5      29                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0

## V11
# grade_0 grade_2_or_higher
# 1    3025                 0
# 2    1332                95
# 3     515                58
# 4     153                27
# 5      28                 5
# 6      11                 1
# 7       5                 1
# 8       1                 0

## V12
# grade_0 grade_2_or_higher
# 1    2549                 0
# 2    1286                94
# 3     508                58
# 4     151                27
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

## V9 (after updating with Mathews data)
# 0    2    3    5 ## Plasma
# 1 3385    0    0    0
# 2 1259   93   28    0
# 3  367   42   10    1
# 4   71   10    2    0
# 5    3    1    0    0

# 0    2    3    5 ## Serum
# 1 3175    0    0    0
# 2  961   82   24    0
# 3  226   23    7    0
# 4   27    8    0    1
# 5    1    0    0    0

## V10
# 0    2    3    5
# 1 3037    0    0    0
# 2 1044   92   28    0
# 3  279   42   10    1
# 4   51   10    2    0
# 5    2    1    0    0

## V11
# 0    2    3    5
# 1 3379    0    0    0
# 2 1253   93   28    0
# 3  365   42   10    1
# 4   70   10    2    0
# 5    3    1    0    0

## V12
# 0    2    3    5
# 1 3173    0    0    0
# 2 1038   92   28    0
# 3  275   42   10    1
# 4   52   10    2    0
# 5    2    1    0    0

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

## V9 (after updating with Mathews data)
# grade_0 grade_2_or_higher ## Plasma
# 1    3385                 0
# 2    1259               121
# 3     367                53
# 4      71                12
# 5       3                 1

# grade_0 grade_2_or_higher ## Serum
# 1    3175                 0
# 2     961               106
# 3     226                30
# 4      27                 9
# 5       1                 0

## V10
# grade_0 grade_2_or_higher
# 1    3037                 0
# 2    1044               120
# 3     279                53
# 4      51                12
# 5       2                 1

## V11
# grade_0 grade_2_or_higher
# 1    3379                 0
# 2    1253               121
# 3     365                53
# 4      70                12
# 5       3                 1

## V12
# grade_0 grade_2_or_higher
# 1    3173                 0
# 2    1038               120
# 3     275                53
# 4      52                12
# 5       2                 1

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

## V9 (after updating with Mathews data)
# grade_0 grade_2_or_higher
# 1    3030                 0 ## Plasma
# 2    1338                62
# 3     516                29
# 4     155                13
# 5      29                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0

# grade_0 grade_2_or_higher ## serum
# 1    2353                 0
# 2    1321                32
# 3     514                26
# 4     156                12
# 5      29                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0


## V10 
# grade_0 grade_2_or_higher
# 1    2499                 0
# 2    1238                61
# 3     485                29
# 4     145                13
# 5      29                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0

## V11
# grade_0 grade_2_or_higher
# 1    3030                 0
# 2    1338                62
# 3     516                29
# 4     155                13
# 5      29                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0

## V12
# grade_0 grade_2_or_higher
# 1    2549                 0
# 2    1286                61
# 3     508                29
# 4     151                13
# 5      29                 2
# 6      11                 0
# 7       5                 0
# 8       1                 0

table(CTCAE.2$event_number,CTCAE.2$grade) ## Actual visit
# 0    2    3    5
# 1 2549    0    0    0
# 2 1286   72   22    0
# 3  508   49    9    0
# 4  151   19    7    1
# 5   29    4    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

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

## V9 (after updating with Mathews data)
# FALSE  TRUE ## Plasma 
# 5085   187
# 4390   145 ## serum

## V10
# FALSE  TRUE 
# 4413   186 

## V11
# FALSE  TRUE 
# 5070   187 

## V12
# FALSE  TRUE 
# 4540   186
dim(CTCAE.2)
# 4506   16
# 4726

write.table(CTCAE.2, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/Plasma_data_processed_v12_after_removing_numvial_0.txt", col.names = T, row.names = F, sep = "\t") ## based on Plasma data
# write.table(CTCAE.2, "Serum_data_processed_v9_after_removing_numvial_0.txt", col.names = T, row.names = F, sep = "\t") ## based on Serum data

# CTCAE.2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/Plasma_data_processed_v12_after_removing_numvial_0.txt", sep = "\t", header = T)