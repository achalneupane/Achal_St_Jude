rm(list= ls())
library(haven)
library(dplyr)
table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_table_2_v12.updated.txt", header = T, sep = "\t") # after removing num vial zero and with 18 or older
dim(table_2)
# 3608   35

#################################################

diag <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat")

#########################
## Zhaoming's controls ##
#########################
# zhaoming.control <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/controls.selected.110.txt", header = F, sep = "\t")
# zhaoming.control$V1 %in% SERUM$sjlid
zhaoming.control.from.mathew <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/forachal_controls.txt", header = T, sep = "\t")
zhaoming.control.from.mathew <- zhaoming.control.from.mathew[grepl("Plasma", zhaoming.control.from.mathew$aliquot_type, ignore.case = T),]
colnames(zhaoming.control.from.mathew)[colnames(zhaoming.control.from.mathew) == "ageatsample"] <- "Sample_age"

## Keep samples with more than 10 vial and alive
num_vial.1 <- zhaoming.control.from.mathew$sjlid[zhaoming.control.from.mathew$num_vials < 2]
# "SJL5107399"
dim(zhaoming.control.from.mathew)
# 113   6
## select by younger age
# Find the 10 sjlid with the oldest sample_age of vials
sorted_data <- zhaoming.control.from.mathew %>%
  dplyr::arrange(Sample_age, )
remove.10.samples <- tail(unique(sorted_data$sjlid),8) # remove 10 older samples, we want everyone close to 18 years of age
remove.10.samples <- c(remove.10.samples, "SJL5125799", "SJL5107399") # this sample had other carcinoma

zhaoming.control.from.mathew <- zhaoming.control.from.mathew[!zhaoming.control.from.mathew$sjlid %in% remove.10.samples,]
length(unique(zhaoming.control.from.mathew$sjlid))
# 100
dim(zhaoming.control.from.mathew)
# 103
# Now get the first Sample_age
first_ageatsample.zhaoming.100 <- zhaoming.control.from.mathew %>%
  arrange(sjlid, Sample_age) %>%  # Sort by sjlid and ageatsample
  group_by(sjlid) %>%              # Group by sjlid
  slice(1)   

first_ageatsample.zhaoming.100$diaggrp <- diag$diaggrp[match(first_ageatsample.zhaoming.100$sjlid, diag$sjlid)]
first_ageatsample.zhaoming.100$ageevent <- diag$diagdt[match(first_ageatsample.zhaoming.100$sjlid, diag$sjlid)]
first_ageatsample.zhaoming.100$selection_group <- "100_community_controls"
first_ageatsample.zhaoming.100$CMP_status <- "No"
first_ageatsample.zhaoming.100 <- first_ageatsample.zhaoming.100[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status", "selection_group")]
dim(first_ageatsample.zhaoming.100)
# 100  9
###########################
## Add primary diagnosis ##
###########################
table_2$diaggrp <- diag$diaggrp[match(table_2$sjlid, diag$sjlid)]

###################
## Get 171 cases ##
###################
CA.171 <- table_2[table_2$grade_2_or_higher =="grade_2_or_higher",]
CA.171 <- unique(CA.171$sjlid)
table_2$CMP_status <- ifelse(table_2$sjlid %in% CA.171, "Yes", "No")
table(table_2$CMP_status)
# No  Yes 
# 3209  399 

## Get the frist ageevent for all samples
table_2.first.event <- table_2 %>%
  arrange(sjlid, ageevent) %>%  # Sort by sjlid and ageevent
  group_by(sjlid) %>%              # Group by sjlid
  slice(1) 

dim(table_2.first.event)
# 2276   37

CA.171 <- table_2.first.event[table_2.first.event$CMP_status == "Yes",]
CA.171$selection_group <- "171_CMP_cases"
CA.171 <- CA.171[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status", "selection_group")]
#################################################################################################
## randomly select 200 Hodgkin lymphoma survivors from all eligible Hodgkin lymphoma survivors ##
#################################################################################################
# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
## CTCAE <- CTCAE[grepl("Cardiomyopathy", CTCAE$condition),]
CTCAE$ageevent <- round(CTCAE$ageevent,1)
dim(CTCAE)
# [1] 913776     62

## Exclude CA.171 and community controls
CTCAE <- CTCAE[!CTCAE$sjlid %in% unique(c(CA.171$sjlid, first_ageatsample.zhaoming.100$sjlid)),] # exclude those in CA.171 samples from the original table
dim(CTCAE)
# [1] 865250     62


TB <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/Achal_survivors_04.04.2024.txt", header = T, sep = "\t") ## Updated by Matt in April
df <- TB

df$ageatsample <- floor(df$ageatsample * 10) / 10
df.original <- df
df <- df %>%
  distinct()
dim(df)  # 21041 ## April version
# 21041
# ## Remove those with vials less than 2
# df <- df %>%
#   dplyr::group_by(sjlid, ageatsample) %>%
#   dplyr::mutate(total_num_vials = sum(num_vials)) %>%
#   ungroup()
# df <- df[df$total_num_vials >= 2,]
# dim(df)
# # 20104
# SERUM <- df[grepl("Serum", df$aliquot_type, ignore.case = T),]
SERUM <- df[grepl("Plasma", df$aliquot_type, ignore.case = T),]
SERUM.original <- SERUM
dim(SERUM.original)
# 10514  7 ## April version
## I see age at serum sample is <18 yrs. Could you identify the samples among 18
#or higher only? Everyone needs to be adults at serum sample.
SERUM <- SERUM[which(SERUM$ageatsample >= 18),]
dim(SERUM)
# 9141    7
# remove vial zero 0 or 1 
SERUM <- SERUM[SERUM$num_vials > 1,] # nrow 7897
# SERUM <- SERUM[SERUM$num_vials >= 1 & SERUM$vitalstatus == "Alive",] # nrow 8502
dim(SERUM)
# 7897
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum/get_matching_rows_from_CTCAE.R")
CTCAE <- CTCAE[CTCAE$sjlid %in% SERUM$sjlid,]
dim(CTCAE)
# 693827     62

# CTCAE.SERUM <- get_rows_with_smaller_sample_age.all(CTCAE, SERUM, 7)
## saveRDS(CTCAE.SERUM, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/HL_non_HL_CTCAE_plasma_vial_ge_2_18yo_and_alive.rds")
## saveRDS(CTCAE.SERUM, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/HL_non_HL_CTCAE_plasma_vial_ge_1_18yo_and_alive.rds")
# saveRDS(CTCAE.SERUM, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/HL_non_HL_CTCAE_plasma_vial_ge_2_18yo.rds")

CTCAE.SERUM <- readRDS("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/HL_non_HL_CTCAE_plasma_vial_ge_2_18yo.rds")
## removing samples that have duplicate ageevent
# CTCAE.SERUM <- CTCAE.SERUM[!CTCAE.SERUM$sjlid %in% c("SJL1225801", "SJL1265801", "SJL1430801", "SJL4730101", "SJL5134305", "SJL5146506"),]
CTCAE.SERUM <- CTCAE.SERUM[CTCAE.SERUM$grade != -9,]
CTCAE.SERUM <- CTCAE.SERUM[!is.na(CTCAE.SERUM$grade),]
dim(CTCAE.SERUM)
# 665398     66

CTCAE.SERUM$diaggrp <- diag$diaggrp[match(CTCAE.SERUM$sjlid, diag$sjlid)]

gg <- CTCAE.SERUM[c("sjlid", "condition", "grade", "ageevent", "Sample_age", "diaggrp")]
CTCAE.SERUM.cardiomyopathy <- CTCAE.SERUM[grepl("cardiomyopathy", CTCAE.SERUM$condition, ignore.case = T),]
# cardio.gg <- CTCAE.SERUM.cardiomyopathy[c("sjlid", "condition", "grade", "ageevent", "Sample_age", "diaggrp")]

## add event number
CTCAE.SERUM.cardiomyopathy <- CTCAE.SERUM.cardiomyopathy %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid) 

# remove all rows as soon as row with max_grade (which is prior grade) is greater than 2
CTCAE.SERUM.cardiomyopathy <- CTCAE.SERUM.cardiomyopathy %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(sjlid, event_number) %>%
  dplyr::mutate(max_grade_prior = cummax(grade))

dim(CTCAE.SERUM.cardiomyopathy)
# 7506   69
gg <- CTCAE.SERUM.cardiomyopathy[c("sjlid", "condition", "grade", "ageevent", "Sample_age", "diaggrp", "num_vials", "max_grade_prior", "event_number")]

CTCAE.SERUM.cardiomyopathy <- CTCAE.SERUM.cardiomyopathy[!(is.na(CTCAE.SERUM.cardiomyopathy$Sample_age) & CTCAE.SERUM.cardiomyopathy$grade == 0),]
dim(CTCAE.SERUM.cardiomyopathy)
# 5235 69

## Remove rows once grades 2 or higher are seen in ordered df by sjlid and event_number
CTCAE.SERUM.cardiomyopathy <- filter_rows_by_condition(CTCAE.SERUM.cardiomyopathy, "sjlid", "max_grade_prior")
dim(CTCAE.SERUM.cardiomyopathy)
# 4698

# This is all
table(CTCAE.SERUM.cardiomyopathy$event_number, CTCAE.SERUM.cardiomyopathy$max_grade_prior)


CTCAE.SERUM.cardiomyopathy <- CTCAE.SERUM.cardiomyopathy %>%
  dplyr::group_by(sjlid) %>%
  dplyr::arrange(ageevent) %>%
  dplyr::mutate(new_event_number = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sjlid) 

gg <- CTCAE.SERUM.cardiomyopathy[c("sjlid", "condition", "grade", "ageevent", "Sample_age", "diaggrp", "num_vials", "max_grade_prior", "event_number", "new_event_number")]

# This is all
table(CTCAE.SERUM.cardiomyopathy$new_event_number, CTCAE.SERUM.cardiomyopathy$max_grade_prior)

## at baseline
CTCAE.SERUM.cardiomyopathy.baseline <- CTCAE.SERUM.cardiomyopathy %>%
  group_by(sjlid) %>%
  mutate(row_number = row_number()) %>%
  filter(max_grade_prior <= 0 | row_number <= which.max(max_grade_prior)) %>%
  select(-row_number)

remove15.at.baseline <- c("SJL1239901", "SJL1261901", "SJL1437801", "SJL1527107", "SJL4177213", "SJL4769616", "SJL5015318", "SJL5016118", "SJL5041805", "SJL5052915", "SJL5103906", "SJL5104317", "SJL5192513", "SJL5205006", "SJL5310713")
CTCAE.SERUM.cardiomyopathy.baseline <- CTCAE.SERUM.cardiomyopathy.baseline[!CTCAE.SERUM.cardiomyopathy.baseline$sjlid %in% remove15.at.baseline,]
# Table: Counts of samples by visit number and grades
table(CTCAE.SERUM.cardiomyopathy.baseline$new_event_number, CTCAE.SERUM.cardiomyopathy.baseline$max_grade_prior)
# 0    2    3    4
# 1 2987  262  125    4
# 2  972    0    0    0
# 3  262    0    0    0
# 4   51    0    0    0
# 5    2    0    0    0

# Table: Counts of samples by CTCAE visit number and max grades prior to or on the date of plasma sampling.
table(CTCAE.SERUM.cardiomyopathy.baseline$event_number, CTCAE.SERUM.cardiomyopathy.baseline$max_grade_prior)
# 0    2    3    4
# 1 2379  207  112    4
# 2 1220   43   10    0
# 3  484    8    1    0
# 4  147    3    1    0
# 5   27    1    1    0
# 6   11    0    0    0
# 7    5    0    0    0
# 8    1    0    0    0




# Yutaka on 04/09/2024: The case cohort design selects a "subcohort" randomly
# from a cohort and then enables comparison of the subcohort with each of
# different-disease "case" groups.  In our case, we wanted a subcohort
# consisting of HL and non HL survivors so that we can do multiple
# case-subcohort comparisons for different diseases.  If we are interested only
# in cardiomyopathy, your selection is good which is a nested case control
# design.  But we want to utilize these precious samples for different diseases
# in future, although there are some technical issues (batch to batch variations
# cannot be randomized to cases and comparison subjects).  But if we have many
# subjects with cardiomyopathy prior to plasma sampling, then we lose power as
# they cannot be used for proteomics/metabolomics analysis for cardiomyopathy.
# So, we need to assess how many survivors had cardiomyopathy before or at
# baseline vs. the table you provided.  I know this is confusing, but can you
# make the table for survivors who were excluded from your table?

## Excluding 171 cases
table_2.original <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_table_2_v12.txt", header = T, sep = "\t") # after removing num vial zero and with 18 or older
# exclude.CA.177 <- table_2.original$sjlid[table_2.original$grade_2_or_higher == "grade_2_or_higher"]
CTCAE.SERUM.cardiomyopathy.baseline.excluded <- CTCAE.SERUM.cardiomyopathy.baseline[!CTCAE.SERUM.cardiomyopathy.baseline$sjlid %in% table_2.original$sjlid,]

# Table: Samples not included in table 2 (all 171 cases and controls), but only in Table 7 max grades prior to or on the date of plasma sampling.
table(CTCAE.SERUM.cardiomyopathy.baseline.excluded$new_event_number, CTCAE.SERUM.cardiomyopathy.baseline.excluded$max_grade_prior)
# 0   2   3   4
# 1 882 262 125   4
# 2 170   0   0   0
# 3  12   0   0   0
# 4   1   0   0   0

# Table: Samples not included in table 1 (all 186 cases and controls), but only in Table 7 max grades prior to or on the date of plasma sampling with actual CTCAE visit number.
table(CTCAE.SERUM.cardiomyopathy.baseline.excluded$event_number, CTCAE.SERUM.cardiomyopathy.baseline.excluded$max_grade_prior)
# 0   2   3   4
# 1 845 207 112   4
# 2 196  43  10   0
# 3  20   8   1   0
# 4   4   3   1   0
# 5   0   1   1   0



## Exclude all grade or samples at baseline
CTCAE.2.original <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/Plasma_data_processed_v12_after_removing_numvial_0.txt", sep = "\t", header = T)
# exclude.CA.186 <- CTCAE.2.original$sjlid[CTCAE.2.original$grade_2_or_higher=="grade_2_or_higher"]
CTCAE.SERUM.cardiomyopathy.baseline.excluded186 <- CTCAE.SERUM.cardiomyopathy.baseline[!CTCAE.SERUM.cardiomyopathy.baseline$sjlid %in% CTCAE.2.original$sjlid,]
CTCAE.SERUM.cardiomyopathy.baseline.excluded186$sjlid[!CTCAE.SERUM.cardiomyopathy.baseline.excluded186$sjlid %in% CTCAE.2.original$sjlid]

sum(CTCAE.SERUM.cardiomyopathy.baseline.excluded186$new_event_number==1 & CTCAE.SERUM.cardiomyopathy.baseline.excluded186$max_grade_prior==0)
## Should be 0!!
CTCAE.SERUM.cardiomyopathy.baseline.excluded186$sjlid[CTCAE.SERUM.cardiomyopathy.baseline.excluded186$new_event_number==1 & CTCAE.SERUM.cardiomyopathy.baseline.excluded186$max_grade_prior==0]
# "SJL1225801" "SJL1265801" "SJL1430801" "SJL4730101" "SJL5134305" "SJL5146506"

table(CTCAE.SERUM.cardiomyopathy.baseline.excluded186$new_event_number, CTCAE.SERUM.cardiomyopathy.baseline.excluded186$max_grade_prior)
# 2   3   4
# 1 262 125   4

table(CTCAE.SERUM.cardiomyopathy.baseline.excluded186$event_number, CTCAE.SERUM.cardiomyopathy.baseline.excluded186$max_grade_prior)
# 2   3   4
# 1 207 112   4
# 2  43  10   0
# 3   8   1   0
# 4   3   1   0
# 5   1   1   0


dim(SERUM)
# 7897    7
length(unique(SERUM$sjlid))
# 4493
SERUM.not.CA171 <- SERUM[!SERUM$sjlid %in% CA.171$sjlid,]
dim(SERUM.not.CA171)
# 7493
length(unique(SERUM.not.CA171$sjlid))
# 4325
## Add diagnosis 
SERUM.not.CA171$diaggrp <- diag$diaggrp[match(SERUM.not.CA171$sjlid, diag$sjlid)]

colnames(SERUM.not.CA171)[colnames(SERUM.not.CA171) == "ageatsample"] <- "Sample_age"
SERUM.not.CA171$ageevent <- NA
SERUM.not.CA171$CMP_status <- "No"

write.table(SERUM.not.CA171, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/SERUM.not.CA171.txt", col.names = T, row.names = F, quote = F, sep = "\t", na="")
############################
## randomly select 200 HL ##
############################
all.hodgkin <- SERUM.not.CA171[grepl("^Hodgkin", SERUM.not.CA171$diaggrp, ignore.case = T),]
dim(all.hodgkin)
# 829  8
length(unique(all.hodgkin$sjlid))
# 447
table(all.hodgkin$diaggrp)

set.seed(54321)
# select 200
all.hodgkin.extract <- sample(all.hodgkin$tb_number, 200)
all.hodgkin.extract <- all.hodgkin[all.hodgkin$tb_number %in% all.hodgkin.extract,]
dim(all.hodgkin.extract)
# 200  37
all.hodgkin.extract$selection_group <- "200_Hodgkin_lymphoma"
all.hodgkin.extract <- all.hodgkin.extract[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status", "selection_group")]
####################################################
## randomly select the remaining 1100-171-200=729 ##
####################################################
all.non.hodgkin <- SERUM.not.CA171[!SERUM.not.CA171$sjlid %in% c(CA.171$sjlid,all.hodgkin$sjlid),] # exclude those in CA.171 + 200 hodgkins samples from the original table
# all.non.hodgkin <- all.non.hodgkin[all.non.hodgkin$diaggrp!="",]
# Randomly select 733 non- Hodgkin lymphoma survivors from all eligible non- Hodgkin lymphoma survivors.
dim(all.non.hodgkin)
# 6664   10
length(unique(all.non.hodgkin$sjlid))
# 3878
table(all.non.hodgkin$diaggrp)

set.seed(54321)
# select 733
all.non.hodgkin.extract <- sample(unique(all.non.hodgkin$tb_number), 729)
all.non.hodgkin.extract <- all.non.hodgkin[all.non.hodgkin$tb_number %in% all.non.hodgkin.extract,]
all.non.hodgkin.extract$selection_group <- "729_non_Hodgkin"
all.non.hodgkin.extract <- all.non.hodgkin.extract[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status", "selection_group")]
dim(all.non.hodgkin.extract)
# 729  9
#################################
## concatenate all wanted ones ##
#################################
all.wanted.df.1200 <- rbind.data.frame(CA.171, all.hodgkin.extract, all.non.hodgkin.extract, first_ageatsample.zhaoming.100)
dim(all.wanted.df.1200)
# 1200    9

write.table(all.wanted.df.1200, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_data_batch1_1200_samples.txt", col.names = T, row.names = F, quote = F, sep = "\t", na="")


# Note from Yadav: I think you should only provide the necessary information
# when you send these files. You would only need tb_number, sjlid, num_vials and
# the group to ECC people. When you send to the proteomics core, you should only
# send them tb_number and the group so that they can include all 4 groups of
# survivors in each experiment.



# Yadav on 3/26/2024: Achal, can you please provide sex, age at sample and race of these samples?
demographic <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
all.wanted.df.1200.to.update <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_data_batch1_1200_samples.txt", header = T, sep = "\t")
all.wanted.df.1200.to.update$Sex <- demographic$gender[match(all.wanted.df.1200.to.update$sjlid, demographic$sjlid)]
all.wanted.df.1200.to.update$racegrp <- demographic$racegrp[match(all.wanted.df.1200.to.update$sjlid, demographic$sjlid)]

all.wanted.df.1200.to.update.proteomics <- all.wanted.df.1200.to.update[c("tb_number", "Sample_age", "selection_group", "Sex", "racegrp")]
#############################
## Give this to proteomics ##
#############################
write.table(all.wanted.df.1200.to.update.proteomics, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_data_batch1_1200_samples_to_proteomics_core.txt", col.names = T, row.names = F, sep = "\t", quote = F)
table(all.wanted.df.1200.to.update$num_vials, all.wanted.df.1200.to.update$selection_group)
# 100_community_controls 171_CMP_cases 200_Hodgkin_lymphoma 729_non_Hodgkin
# 1                      0            26                    0               0
# 2                      0            96                   84             287
# 3                      2             0                    1               3
# 4                      9             2                   16               5
# 5                      5             6                   29              36
# 6                     84            41                   69             398
# 7                      0             0                    1               0

###############################
## all SERUM and PLASMA data ##
###############################
TB <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/Achal_survivors_04.04.2024.txt", header = T, sep = "\t") ## Updated by Matt in April
TB <- TB[c("sjlid", "num_vials", "aliquot_type", "ageatsample", "vitalstatus",  "tb_number")]
TB.control <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/forachal_controls.txt", header = T, sep = "\t")
TB <- rbind.data.frame(TB, TB.control)
df <- TB
df$ageatsample <- floor(df$ageatsample * 10) / 10
df.original <- df
df <- df %>%
  distinct()
dim(df)  # 20137
# 20779 # Mathew ## March
# 21267 ## April version

df <- df[which(df$num_vials > 0),]
total_vials.by.sjlid <- aggregate(num_vials ~ sjlid, data = df, FUN = sum)

df$KEY <- paste0(df$sjlid, ":", df$ageatsample)
PLASMA <- df[grepl("Plasma", df$aliquot_type, ignore.case = T),]
SERUM <- df[grepl("Serum", df$aliquot_type, ignore.case = T),]


total_plasma.vials.by.sjlid <- aggregate(num_vials ~ KEY, data = PLASMA, FUN = sum)
total_serum.vials.by.sjlid <- aggregate(num_vials ~ KEY, data = SERUM, FUN = sum)

##############################
## Give this to Kyla's team ##
##############################
all.wanted.df.1200.to.update.kyla <- all.wanted.df.1200.to.update[c("tb_number", "sjlid", "Sample_age", "num_vials", "vitalstatus","selection_group")]
write.table(all.wanted.df.1200.to.update.kyla, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_data_batch1_1200_samples_to_ECC.txt", col.names = T, row.names = F, sep = "\t", quote = F)



#########################
## add total num_vials ##
#########################
all.wanted.df.1200.to.update.kyla$total_vials_for_all_visits <- total_vials.by.sjlid$num_vials[match(all.wanted.df.1200.to.update.kyla$sjlid, total_vials.by.sjlid$sjlid)]
all.wanted.df.1200.to.update.kyla$KEY <- paste0(all.wanted.df.1200.to.update.kyla$sjlid, ":", all.wanted.df.1200.to.update.kyla$Sample_age)

all.wanted.df.1200.to.update.kyla$total_plasma_for_the_visit <- total_plasma.vials.by.sjlid$num_vials[match(all.wanted.df.1200.to.update.kyla$KEY, total_plasma.vials.by.sjlid$KEY)]
all.wanted.df.1200.to.update.kyla$total_serum_for_the_visit <- total_serum.vials.by.sjlid$num_vials[match(all.wanted.df.1200.to.update.kyla$KEY, total_serum.vials.by.sjlid$KEY)]

all.wanted.df.1200.to.update.kyla$total_plasma_for_the_visit[is.na(all.wanted.df.1200.to.update.kyla$total_plasma_for_the_visit)] <- 0
all.wanted.df.1200.to.update.kyla$total_serum_for_the_visit[is.na(all.wanted.df.1200.to.update.kyla$total_serum_for_the_visit)] <- 0
all.wanted.df.1200.to.update.kyla$visit_depleted_YN <- ifelse((all.wanted.df.1200.to.update.kyla$total_plasma_for_the_visit + all.wanted.df.1200.to.update.kyla$total_serum_for_the_visit) >= 2, "No", "Yes")
write.table(all.wanted.df.1200.to.update.kyla, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_data_batch1_1200_samples_to_ECC_with_vials_info.txt", col.names = T, row.names = F, sep = "\t", quote = F)


## EMAIL from Matt on 4/4/2024: There are a few discrepancies and there are 33
#samples (attached) that will now be depleted if used. The reason for this
#discrepancy is because the most updated biorepository (frozen on April 2nd)
#has a less vial count compared to the biorepository used to make this list
#(frozen on March 11th).

new.missing <- read_sas("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/diff.sas7bdat")
new.missing$selection_group <- all.wanted.df.1200.to.update$selection_group[match(new.missing$tb_number, all.wanted.df.1200.to.update$tb_number)]

# ## Test
# library("blockrand")
# randomized_samples <- block_ra(sample_data, n = 14, id_col = "Sample", block_col = "selection_group", strata_cols = c("Sample_age", "Sex", "racegrp"))
# 
# # View the randomized samples
# randomized_samples


