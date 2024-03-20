library(haven)
library(dplyr)
table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/plasma_table_2_v9.txt", header = T, sep = "\t") # after removing num vial zero and with 18 or older
diag <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat")

#########################
## Zhaoming's controls ##
#########################
# zhaoming.control <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/controls.selected.110.txt", header = F, sep = "\t")
# zhaoming.control$V1 %in% SERUM$sjlid
zhaoming.control.from.mathew <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/forachal_controls.txt", header = T, sep = "\t")
zhaoming.control.from.mathew <- zhaoming.control.from.mathew[grepl("Plasma", zhaoming.control.from.mathew$aliquot_type, ignore.case = T),]
colnames(zhaoming.control.from.mathew)[colnames(zhaoming.control.from.mathew) == "ageatsample"] <- "Sample_age"

## select by younger age
# Find the 10 sjlid with the oldest sample_age of vials
sorted_data <- zhaoming.control.from.mathew %>%
  arrange(Sample_age)
remove.10.samples <- tail(unique(sorted_data$sjlid),9) # remove 10 older samples, we want everyone close to 18 years of age
remove.10.samples <- c(remove.10.samples, "SJL5125799") # this sample had other carcinoma

zhaoming.control.from.mathew <- zhaoming.control.from.mathew[!zhaoming.control.from.mathew$sjlid %in% remove.10.samples,]
length(unique(zhaoming.control.from.mathew$sjlid))
# 100

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

## Get the frist ageevent for all samples
table_2.first.event <- table_2 %>%
  arrange(sjlid, ageevent) %>%  # Sort by sjlid and ageevent
  group_by(sjlid) %>%              # Group by sjlid
  slice(1) 

dim(table_2.first.event)
# 2433   37

CA.171 <- table_2.first.event[table_2.first.event$CMP_status == "Yes",]
CA.171$selection_group <- "171_CMP_cases"
CA.171 <- CA.171[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status", "selection_group")]
#################################################################################################
## randomly select 200 Hodgkin lymphoma survivors from all eligible Hodgkin lymphoma survivors ##
#################################################################################################
table_3 <- table_2.first.event[!table_2.first.event$sjlid %in% CA.171$sjlid,] # exclude those in CA.171 samples from the original table
all.hodgkin <- table_3[grepl("^Hodgkin", table_3$diaggrp, ignore.case = T),]
dim(all.hodgkin)
# 346  37

set.seed(54321)
# select 200
all.hodgkin.extract <- sample(all.hodgkin$tb_number, 200)
all.hodgkin.extract <- table_3[table_3$tb_number %in% all.hodgkin.extract,]
dim(all.hodgkin.extract)
# 200  37
all.hodgkin.extract$selection_group <- "200_Hodgkin_lymphoma"
all.hodgkin.extract <- all.hodgkin.extract[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status", "selection_group")]
####################################################
## randomly select the remaining 1100-171-200=729 ##
####################################################
table_3 <- table_2.first.event[!table_2.first.event$tb_number %in% c(CA.171$tb_number,all.hodgkin.extract$tb_number),] # exclude those in CA.171 + 200 hodgkins samples from the original table
# Randomly select 733 non- Hodgkin lymphoma survivors from all eligible non- Hodgkin lymphoma survivors.
all.non.hodgkin <- table_3[!table_3$tb_number %in% all.hodgkin$tb_number,]

set.seed(54321)
# select 733
all.non.hodgkin.extract <- sample(unique(all.non.hodgkin$tb_number), 729)
all.non.hodgkin.extract <- table_3[table_3$tb_number %in% all.non.hodgkin.extract,]
all.non.hodgkin.extract$selection_group <- "729_non_Hodgkin"
all.non.hodgkin.extract <- all.non.hodgkin.extract[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status", "selection_group")]
#################################
## concatenate all wanted ones ##
#################################
all.wanted.df.1200 <- rbind.data.frame(CA.171, all.hodgkin.extract, all.non.hodgkin.extract, first_ageatsample.zhaoming.100)
dim(all.wanted.df.1200)
# 1200    9
#########################
## Now, extract TB IDs ##
#########################
extract.1200.samples <- table_2.first.event[table_2.first.event$sjlid %in% all.wanted.sjlids,]


