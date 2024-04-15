## Yadav's email on 10/23/2023 
# Could you provide the summary as in “Sheet1” for the following subset of survivors?

# Exposed to either anthracyclines or chest radiation
# With WGS (no restriction on treatment exposures)
# Exposed to either anthracyclines or chest radiation with WGS
# "SJL1239901" "SJL1261901" "SJL1527107" "SJL4769616" "SJL5015318" "SJL5103906"

library(dplyr)
library(tidyr)
library(haven)
setwd("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/")
# CTCAE.serum <- read.table("Serum_data_processed_v6.txt", header = T, sep = "\t")
# CTCAE.serum <- read.table("Serum_data_processed_v8.txt", header = T, sep = "\t") # removing 18 or younger

# CTCAE.serum <- read.table("Serum_data_processed_v9_after_removing_numvial_0.txt", header = T, sep = "\t") # Serum, 18 or older, vial > 0
CTCAE.serum <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/Plasma_data_processed_v12_after_removing_numvial_0.txt", header = T, sep = "\t") # Plasma, 18 or older, vial > 0

demographic <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
chemo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat")
chemo <- chemo[c("sjlid", "anthracyclines_dose_any",	"anthracyclines_dose_prim",	"anthracyclines_dose_5",	"anthracyclines_dose_10")]

radiation <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat")
radiation <- radiation[c("sjlid", "maxchestrtdose")]
radiation.2 <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiationsum_yn.sas7bdat")
radiation.2 <- radiation.2[c("sjlid", "Chest", "Chest_prim", "Chest_5", "Chest_10")]

have.WGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr8.preQC_biallelic_renamed_ID_edited.vcf.gz.fam")

# demographic <- cbind.data.frame(demographic, chemo[match(demographic$sjlid, chemo$sjlid),])
demographic <- cbind.data.frame(demographic, chemo[match(demographic$sjlid, chemo$sjlid), "anthracyclines_dose_any"])
# demographic <- cbind.data.frame(demographic, radiation[match(demographic$sjlid, radiation$sjlid),])
demographic <- cbind.data.frame(demographic, radiation[match(demographic$sjlid, radiation$sjlid),"maxchestrtdose"])
# demographic <- cbind.data.frame(demographic, radiation.2[match(demographic$sjlid, radiation.2$sjlid),])
demographic <- cbind.data.frame(demographic, radiation.2[match(demographic$sjlid, radiation.2$sjlid),c("Chest", "Chest_prim", "Chest_5", "Chest_10")])
demographic$WGS <- ifelse(demographic$sjlid %in% have.WGS$V2, "Yes", 'No')

demographic <- demographic[match(CTCAE.serum$sjlid, demographic$sjlid),]
demographic <- demographic[!colnames(demographic) %in% colnames(CTCAE.serum)]

CTCAE.serum <- cbind.data.frame(CTCAE.serum,demographic)
table(CTCAE.serum$new_event_number, CTCAE.serum$grade)
# 0    2    3    5
# 1 4085    0    0    0
# 2 1296   89   34    0
# 3  383   43   13    0
# 4   74   11    1    1
# 5    3    1    0    0

# 0    2    3    5
# 1 3402    0    0    0
# 2 1258   86   28    0
# 3  383   42   13    0
# 4   74   11    1    1
# 5    3    1    0    0

# Serum after removing num vial zero
# 0    2    3    5
# 1 3174    0    0    0
# 2  960   82   24    0
# 3  226   23    7    0
# 4   27    8    0    1
# 5    1    0    0    0

# Plasma after removing num vial zero
# 0    2    3    5
# 1 3384    0    0    0
# 2 1258   93   28    0
# 3  367   42   10    1
# 4   71   10    2    0
# 5    3    1    0    0


## V9 (after updating with Mathews data) 
# 0    2    3    5 ## plasma
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

# Exposed to either anthracyclines or chest radiation
CTCAE.serum.1 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 | CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade)
# 0    2    3    5
# 1 2799    0    0    0
# 2 1113   79   30    0
# 3  369   43   12    0
# 4   73   10    1    1
# 5    3    1    0    0

# 0    2    3    5 # after removing younger serum < 18 yo
# 1 2437    0    0    0
# 2 1076   76   24    0
# 3  369   42   12    0
# 4   73   10    1    1
# 5    3    1    0    0

# Serum Num vial greater than 0
     # 0    2    3    5
# 1 2228    0    0    0
# 2  787   72   19    0
# 3  218   23    7    0
# 4   27    7    0    1
# 5    1    0    0    0

# Plasma Num vial greater than 0
# 0    2    3    5
# 1 2432    0    0    0
# 2 1078   83   24    0
# 3  352   42    9    1
# 4   70    9    2    0
# 5    3    1    0    0

## V9 (after updating with Mathews data) 
# 0    2    3    5 ## plasma
# 1 2433    0    0    0
# 2 1079   83   24    0
# 3  352   42    9    1
# 4   70    9    2    0
# 5    3    1    0    0

# 0    2    3    5 ## Serum
# 1 2229    0    0    0
# 2  788   72   19    0
# 3  218   23    7    0
# 4   27    7    0    1
# 5    1    0    0    0

## V10 
# 0    2    3    5
# 1 2176    0    0    0
# 2  880   83   24    0
# 3  268   42    9    1
# 4   50    9    2    0
# 5    2    1    0    0

## V11 **
# 0    2    3    5
# 1 2430    0    0    0
# 2 1076   83   24    0
# 3  350   42    9    1
# 4   69    9    2    0
# 5    3    1    0    0

## v12
# 0    2    3    5
# 1 2276    0    0    0
# 2  866   83   24    0
# 3  262   42    9    1
# 4   51    9    2    0
# 5    2    1    0    0

table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade >= 2)
# FALSE TRUE
# 1  2799    0
# 2  1113  109
# 3   369   55
# 4    73   12
# 5     3    1

# FALSE TRUE
# 1  2437    0
# 2  1076  100
# 3   369   54
# 4    73   12
# 5     3    1

# Serum Num vial greater than 0
# FALSE TRUE
# 1  2228    0
# 2   787   91
# 3   218   30
# 4    27    8
# 5     1    0

# Plasma Num vial greater than 0 
# FALSE TRUE
# 1  2432    0
# 2  1078  107
# 3   352   52
# 4    70   11
# 5     3    1

## V9 (after updating with Mathews data) 
# FALSE TRUE ## plasma
# 1  2433    0
# 2  1079  107
# 3   352   52
# 4    70   11
# 5     3    1

# FALSE TRUE ## serum
# 1  2229    0
# 2   788   91
# 3   218   30
# 4    27    8
# 5     1    0

## V10
# FALSE TRUE
# 1  2176    0
# 2   880  107
# 3   268   52
# 4    50   11
# 5     2    1

## V11
# FALSE TRUE
# 1  2430    0
# 2  1076  107
# 3   350   52
# 4    69   11
# 5     3    1

## V12
# FALSE TRUE
# 1  2276    0
# 2   866  107
# 3   262   52
# 4    51   11
# 5     2    1

## For Table 2 (in Table_counts sheet of Serum_data_processed_v6_corrected_11_18_2023.xlxs)
table_2 <- CTCAE.serum.1[!is.na(CTCAE.serum.1$grade),]
dim(table_2)
# 4099   35
## 3628 35 ## V12
### write.table(table_2, "serum_table_2_v9.txt", row.names = F, col.names = T, quote = F, sep = "\t")  # Serum, 18 or older, vial > 0
write.table(table_2, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_table_2_v12.txt", row.names = F, col.names = T, quote = F, sep = "\t")  # Plasma, 18 or older, vial > 0
# table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_table_2_v12.txt", header = T, sep = "\t")

table(table_2$new_event_number,table_2$grade)
## V12
# 0    2    3    5
# 1 2276    0    0    0
# 2  866   83   24    0
# 3  262   42    9    1
# 4   51    9    2    0
# 5    2    1    0    0

table(table_2$event_number,table_2$grade) ## Actual CTCAE visit
# 0    2    3    5
# 1 1691    0    0    0
# 2 1088   66   19    0
# 3  486   47    8    0
# 4  146   18    6    1
# 5   29    3    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0

##############################################################
## Exclude any vial 1 if other visits have more than 1 vial ##
##############################################################
###############################################################################################
## We will also remove 1 vial samples if later visits have more vials (especially for cases) ##
###############################################################################################
cases <- table_2[table_2$grade_2_or_higher == "grade_2_or_higher",]
cases <- table_2[table_2$sjlid %in% cases$sjlid,]

controls <- table_2[!table_2$sjlid %in% cases$sjlid,]
table(controls$num_vials) ## should have zero with num_vials = 1
# 2    3    4    5    6    7 
# 860   26   76  312 1934    1 

num_vial1 <- cases[which(cases$num_vials == 1),]
cases.samples <- cases[cases$sjlid %in% num_vial1$sjlid,]


# Filter rows with grade 0
# Initialize a list to store tb_numbers
tb.to.remove <- list()

# Loop through each sjlid
for (sjlid in unique(cases.samples$sjlid)) {
  # Subset data for the current sjlid
  sjlid_data <- cases.samples[cases.samples$sjlid == sjlid, ]
  
  # Check if there are multiple rows with grade 0 before grade 2 or higher
  if((nrow(sjlid_data) >=3) & sjlid_data$num_vials[1] == 1){
    tb.to.remove.tmp <- sjlid_data$tb_number[1]
  } else {
    next
  }
  tb.to.remove <- c(tb.to.remove,tb.to.remove.tmp)
}

# Display the extracted tb_numbers
tb.to.remove <- unlist(tb.to.remove)



table_2.updated <- table_2[!table_2$tb_number %in% tb.to.remove,]
dim(table_2.updated)
# 3671   35 ## V11
# 3608   35 ## v12

table(table_2.updated$new_event_number,table_2.updated$grade)
## V12
# 0    2    3    5
# 1 2256    0    0    0
# 2  866   83   24    0
# 3  262   42    9    1
# 4   51    9    2    0
# 5    2    1    0    0

table(table_2.updated$event_number,table_2.updated$grade) # based on actual CTCAE data
# 0    2    3    5
# 1 1671    0    0    0
# 2 1088   66   19    0
# 3  486   47    8    0
# 4  146   18    6    1
# 5   29    3    1    0
# 6   11    0    1    0
# 7    5    1    0    0
# 8    1    0    0    0




table(table_2.updated$num_vials, table_2.updated$grade)
# 0    2    3    5
# 1   26    0    1    0
# 2  957    3    0    0
# 3   26    1    0    0
# 4   80    2    1    0
# 5  331   22    4    0
# 6 2016   56   10    0
# 7    1    0    0    0

write.table(table_2.updated, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_table_2_v12.updated.txt", row.names = F, col.names = T, quote = F, sep = "\t")  # Plasma, 18 or older, vial > 0

#######################################################################################################
## END
#######################################################################################################




















## break down exposed to either anthracyclines or chest radiation by race
CTCAE.serum.1$race_group [CTCAE.serum.1$race =="White"] <- "White"
CTCAE.serum.1$race_group [CTCAE.serum.1$race =="Black"] <- "Black"
CTCAE.serum.1$race_group[is.na(CTCAE.serum.1$race_group)] <- "Other"
table(CTCAE.serum.1$race_group)
# Black Other White 
# 577   189  3839 

# Black Other White 
# 423   135  2888

## Only in White
CTCAE.serum.1.w <- CTCAE.serum.1[CTCAE.serum.1$race_group =="White",]
table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade)
# 0    2    3    5
# 1 2340    0    0    0
# 2  964   67   22    0
# 3  321   39   10    0
# 4   64    6    1    1
# 5    3    1    0    0

# # After removing num vial = 0
# 0    2    3    5
# 1 1885    0    0    0
# 2  679   64   17    0
# 3  187   19    5    0
# 4   26    4    0    1
# 5    1    0    0    0

table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade >= 2)
# FALSE TRUE
# 1  2340    0
# 2   964   89
# 3   321   49
# 4    64    8
# 5     3    1

# after removing num vial = 0
# FALSE TRUE
# 1  1885    0
# 2   679   81
# 3   187   24
# 4    26    5
# 5     1    0

## Only in Black
CTCAE.serum.1.w <- CTCAE.serum.1[CTCAE.serum.1$race_group =="Black",]
table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade)
# 0   2   3
# 1 371   0   0
# 2 124  12   7
# 3  45   4   2
# 4   8   4   0

# after removing num vial zero
# 0   2   3
# 1 283   0   0
# 2  90   8   2
# 3  30   4   2
# 4   1   3   0

table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade >= 2)
# FALSE TRUE
# 1   371    0
# 2   124   19
# 3    45    6
# 4     8    4

## Only in other races
CTCAE.serum.1.w <- CTCAE.serum.1[CTCAE.serum.1$race_group =="Other",]
table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade)
# 0  3
# 1 88  0
# 2 25  1
# 3  3  0
# 4  1  0

table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade >= 2)
# FALSE TRUE
# 1    88    0
# 2    25    1
# 3     3    0
# 4     1    0

## exposed to anthracyclines with chest RT
CTCAE.serum.4 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade)
# 0   2   3
# 1 538   0   0
# 2 240  21  13
# 3  96  17   2
# 4  24   4   1

table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade >= 2)
# FALSE TRUE
# 1   538    0
# 2   240   34
# 3    96   19
# 4    24    5

## exposed to anthracyclines without chest RT
CTCAE.serum.5 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose <= 200,]
table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade)
# 0    2    3
# 1 1849    0    0
# 2  711   42   14
# 3  206   22    7
# 4   34    3    0
# 5    3    1    0

table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade >= 2)
# FALSE TRUE
# 1  1849    0
# 2   711   56
# 3   206   29
# 4    34    3
# 5     3    1

## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0,]
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade)
# 0    2    3
# 1 2425    0    0
# 2  957   64   27
# 3  304   39    9
# 4   58    8    1
# 5    3    1    0

table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade >= 2)
# FALSE TRUE
# 1  2425    0
# 2   957   91
# 3   304   48
# 4    58    9
# 5     3    1

## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0,]
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade)
# 0    2    3
# 1 2425    0    0
# 2  957   64   27
# 3  304   39    9
# 4   58    8    1
# 5    3    1    0

table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade >= 2)
# FALSE TRUE
# 1  2425    0
# 2   957   91
# 3   304   48
# 4    58    9
# 5     3    1

## break down exposed to anthracyclines regardless of chest RT by race
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="White"] <- "White"
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="Black"] <- "Black"
CTCAE.serum.6$race_group[is.na(CTCAE.serum.6$race_group)] <- "Other"
table(CTCAE.serum.6$race_group)
# Black Other White 
# 510   114  3290 


## Only in White
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="White",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
# 0    2    3
# 1 2022    0    0
# 2  828   53   21
# 3  263   35    8
# 4   50    5    1
# 5    3    1    0


table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)
# FALSE TRUE
# 1  2022    0
# 2   828   74
# 3   263   43
# 4    50    6
# 5     3    1

## Only in Black
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Black",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
# 0   2   3
# 1 332   0   0
# 2 109  11   5
# 3  38   4   1
# 4   7   3   0

table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)
# FALSE TRUE
# 1   332    0
# 2   109   16
# 3    38    5
# 4     7    3

## Only in other races
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Other",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
# 0  3
# 1 71  0
# 2 20  1
# 3  3  0
# 4  1  0

table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)
# FALSE TRUE
# 1    71    0
# 2    20    1
# 3     3    0
# 4     1    0

######################################################
## With WGS (no restriction on treatment exposures) ##
######################################################
CTCAE.serum.2 <- CTCAE.serum[CTCAE.serum$WGS == "Yes",]
table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade)
# 0    2    3    5
# 1 3401    0    0    0
# 2 1223   84   29    0
# 3  368   37   13    0
# 4   70   11    1    1
# 5    3    1    0    0

table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade >= 2)
# FALSE TRUE
# 1  3401    0
# 2  1223  113
# 3   368   50
# 4    70   13
# 5     3    1



# Exposed to either anthracyclines or chest radiation with WGS
CTCAE.serum.3 <- CTCAE.serum.2[CTCAE.serum.2$anthracyclines_dose_any > 0 | CTCAE.serum.2$maxchestrtdose > 200,]
table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade)
# 0    2    3    5
# 1 2328    0    0    0
# 2 1040   75   26    0
# 3  354   37   12    0
# 4   69   10    1    1
# 5    3    1    0    0




table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade >= 2)
# FALSE TRUE
# 1  2328    0
# 2  1040  101
# 3   354   49
# 4    69   12
# 5     3    1





# # CTCAE.original <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
# CTCAE.original.2 <- CTCAE.original[grepl("Cardiomyopathy", CTCAE.original$condition),]
# CTCAE.original.2 <- CTCAE.original.2[c("sjlid", "studypop", "sjlife_cohort", "gender", "organsys", "condition", "gradedt", "grade", "ageevent")]
# # Also, we previously performed proteomics/metabolomics experiments on 200
# # survivors (100 with cardiomyopathy); they are attached here (see Sheet3 for
# # relevant information). Could you check how many of these overlap with the
# # current ~3500 or so (based on grade, age at assessment and sample age)? Also,
# # check if these 200 samples fulfill the same criteria applied to the 3500
# # samples (and if not, how many meet and how many do not).
# library(readxl)
# # Read data from Sheet3
# proteomics <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/pilot_study_200_samples.xlsx", sheet = "Sheet3")
# colnames(proteomics) <- c("sjlid", "grade", "ageevent", "TBID", "Sample_age", "status")
# sum(proteomics$sjlid %in% CTCAE.serum$sjlid)
# # 153
# 
# 
# # Define a function to check if two values are within X-days window
# within_window <- function(value1, value2) {
#   return(abs(value1 - value2) <= 20/365.25)
# }
# 
# # Perform the merge on sjlid, grade, and ageevent
# merged_data <- merge(proteomics, CTCAE.serum, by = c("sjlid", "grade"), suffixes = c(".x", ".y"))
# 
# # Subset the merged data based on your criteria
# overlap_data <- subset(merged_data, within_window(Sample_age.x, Sample_age.y) & within_window(ageevent.x, ageevent.y))
# 
# 
# # Overlaps
# dim(overlap_data)
# # 57
# 
# proteomics.ca <- proteomics[proteomics$status == 1,]
# proteomics.co <- proteomics[proteomics$status == 0,]
# 
# sum(proteomics.ca$sjlid %in% CTCAE.serum$sjlid) ## 59
# sum(proteomics.co$sjlid %in% CTCAE.serum$sjlid) ## 98
# 
# SERUM.kyla <- read.table("Trans-omics CMP profiling Inventory 20230901.txt", header = T)
# SERUM.kyla <- SERUM.kyla[SERUM.kyla$aliquot_type == "Serum",]
# sum(proteomics.ca$sjlid %in% SERUM.kyla$sjlid) ## 100
# sum(proteomics.co$sjlid %in% SERUM.kyla$sjlid) ## 99
# 
# 
# table(CTCAE.serum$new_event_number,CTCAE.serum$grade)
# table(CTCAE.serum$new_event_number,CTCAE.serum$grade >= 2)
# 
# # Calculate age different from first visit to first CMP
# result <- calculate_ageevent_difference(CTCAE.serum)
# print(result)

library(readxl)
## Read data from excel to exclude previous study
# previous.study <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/Serum_data_processed_v6_corrected_10_5_2023.xlsx", sheet = "Serum_data_processed_v6")
previous.study <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/Serum_data_processed_v8_corrected_02_19_2024.xlsx", sheet = "Serum_data_processed_v8")

# Yadav On 11/28/2023:
# Can you please add the summary Table for survivors exposed to anthracyclines (with or without chest RT)? No need to restrict to with WGS only.
CTCAE.serum.saved <- CTCAE.serum
CTCAE.serum <- cbind.data.frame(previous.study, CTCAE.serum[match(previous.study$sjlid, CTCAE.serum$sjlid),c("race", "anthracyclines_dose_any", "maxchestrtdose", "Chest", "Chest_prim", "Chest_5", "Chest_10", "WGS")])
CTCAE.serum$new_event_number <- CTCAE.serum$Order
## Exclude previous ones
CTCAE.serum <- CTCAE.serum[!grepl("Yes", CTCAE.serum$Excluded),]

table(CTCAE.serum$new_event_number, CTCAE.serum$grade)
# 0    2    3    5
# 1 4018    0    0    0
# 2 1239   73   31    0
# 3  383   43   13    0
# 4   74   11    1    1
# 5    3    1    0    0

# Exposed to either anthracyclines or chest radiation
CTCAE.serum.1 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 | CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade)
# 0    2    3    5
# 1 2733    0    0    0
# 2 1056   63   27    0
# 3  369   43   12    0
# 4   73   10    1    1
# 5    3    1    0    0

table(CTCAE.serum.1$new_event_number,CTCAE.serum.1$grade >= 2)
# FALSE TRUE
# 1  2733    0
# 2  1056   90
# 3   369   55
# 4    73   12
# 5     3    1

## break down exposed to either anthracyclines or chest radiation by race
CTCAE.serum.1$race_group [CTCAE.serum.1$race =="White"] <- "White"
CTCAE.serum.1$race_group [CTCAE.serum.1$race =="Black"] <- "Black"
CTCAE.serum.1$race_group[is.na(CTCAE.serum.1$race_group)] <- "Other"
table(CTCAE.serum.1$race_group)
# Black Other White 
# 562   189  3712 

## Only in White
CTCAE.serum.1.w <- CTCAE.serum.1[CTCAE.serum.1$race_group =="White",]
table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade)
# 0    2    3    5
# 1 2286    0    0    0
# 2  909   52   19    0
# 3  321   39   10    0
# 4   64    6    1    1
# 5    3    1    0    0

table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade >= 2)
# FALSE TRUE
# 1  2286    0
# 2   909   71
# 3   321   49
# 4    64    8
# 5     3    1

## Only in Black
CTCAE.serum.1.w <- CTCAE.serum.1[CTCAE.serum.1$race_group =="Black",]
table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade)
# 0   2   3
# 1 359   0   0
# 2 122  11   7
# 3  45   4   2
# 4   8   4   0
table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade >= 2)
# FALSE TRUE
# 1   359    0
# 2   122   18
# 3    45    6
# 4     8    4

## Only in other races
CTCAE.serum.1.w <- CTCAE.serum.1[CTCAE.serum.1$race_group =="Other",]
table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade)
# 0  3
# 1 88  0
# 2 25  1
# 3  3  0
# 4  1  0

table(CTCAE.serum.1.w$new_event_number,CTCAE.serum.1.w$grade >= 2)
# FALSE TRUE
# 1    88    0
# 2    25    1
# 3     3    0
# 4     1    0

## exposed to anthracyclines with chest RT
CTCAE.serum.4 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose > 200,]
table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade)
# 0   2   3
# 1 525   0   0
# 2 224  18  12
# 3  96  17   2
# 4  24   4   1

table(CTCAE.serum.4$new_event_number,CTCAE.serum.4$grade >= 2)
# FALSE TRUE
# 1   525    0
# 2   224   30
# 3    96   19
# 4    24    5

## exposed to anthracyclines without chest RT
CTCAE.serum.5 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0 & CTCAE.serum$maxchestrtdose <= 200,]
table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade)
# 0    2    3
# 1 1796    0    0
# 2  671   29   12
# 3  206   22    7
# 4   34    3    0
# 5    3    1    0

table(CTCAE.serum.5$new_event_number,CTCAE.serum.5$grade >= 2)
# FALSE TRUE
# 1  1796    0
# 2   671   41
# 3   206   29
# 4    34    3
# 5     3    1

## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6 <- CTCAE.serum[CTCAE.serum$anthracyclines_dose_any > 0,]
table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade)
# 0    2    3
# 1 2359    0    0
# 2  900   48   24
# 3  304   39    9
# 4   58    8    1
# 5    3    1    0


table(CTCAE.serum.6$new_event_number,CTCAE.serum.6$grade >= 2)
# FALSE TRUE
# 1  2359    0
# 2   900   72
# 3   304   48
# 4    58    9
# 5     3    1

## exposed to anthracyclines regardless of chest RT
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="White"] <- "White"
CTCAE.serum.6$race_group [CTCAE.serum.6$race =="Black"] <- "Black"
CTCAE.serum.6$race_group[is.na(CTCAE.serum.6$race_group)] <- "Other"
table(CTCAE.serum.6$race_group)
# Black Other White 
# 495   114  3163 

## Only in White
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="White",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
# 0    2    3
# 1 1968    0    0
# 2  773   38   18
# 3  263   35    8
# 4   50    5    1
# 5    3    1    0

table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)
# FALSE TRUE
# 1  1968    0
# 2   773   56
# 3   263   43
# 4    50    6
# 5     3    1

## Only in Black
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Black",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
# 0   2   3
# 1 320   0   0
# 2 107  10   5
# 3  38   4   1
# 4   7   3   0

table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)
# FALSE TRUE
# 1   320    0
# 2   107   15
# 3    38    5
# 4     7    3

## Only in other races
CTCAE.serum.6.w <- CTCAE.serum.6[CTCAE.serum.6$race_group =="Other",]
table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade)
# 0  3
# 1 71  0
# 2 20  1
# 3  3  0
# 4  1  0

table(CTCAE.serum.6.w$new_event_number,CTCAE.serum.6.w$grade >= 2)
# FALSE TRUE
# 1    71    0
# 2    20    1
# 3     3    0
# 4     1    0

######################################################
## With WGS (no restriction on treatment exposures) ##
######################################################
CTCAE.serum.2 <- CTCAE.serum[CTCAE.serum$WGS == "Yes",]
table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade)
# 0    2    3    5
# 1 3335    0    0    0
# 2 1166   68   26    0
# 3  368   37   13    0
# 4   70   11    1    1
# 5    3    1    0    0

table(CTCAE.serum.2$new_event_number,CTCAE.serum.2$grade >= 2)
# FALSE TRUE
# 1  3335    0
# 2  1166   94
# 3   368   50
# 4    70   13
# 5     3    1

# Exposed to either anthracyclines or chest radiation with WGS
CTCAE.serum.3 <- CTCAE.serum.2[CTCAE.serum.2$anthracyclines_dose_any > 0 | CTCAE.serum.2$maxchestrtdose > 200,]
table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade)
# 0    2    3    5
# 1 2263    0    0    0
# 2  983   59   23    0
# 3  354   37   12    0
# 4   69   10    1    1
# 5    3    1    0    0

table(CTCAE.serum.3$new_event_number,CTCAE.serum.3$grade >= 2)
# FALSE TRUE
# 1  2263    0
# 2   983   82
# 3   354   49
# 4    69   12
# 5     3    1