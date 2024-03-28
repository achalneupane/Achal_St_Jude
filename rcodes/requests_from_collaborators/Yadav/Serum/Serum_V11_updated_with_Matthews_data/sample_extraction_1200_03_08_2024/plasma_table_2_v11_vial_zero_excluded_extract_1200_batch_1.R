library(haven)
library(dplyr)
table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v10_output/plasma_table_2_v10.txt", header = T, sep = "\t") # after removing num vial zero and with 18 or older
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
# 112   6
## select by younger age
# Find the 10 sjlid with the oldest sample_age of vials
sorted_data <- zhaoming.control.from.mathew %>%
  arrange(Sample_age, )
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

#######################################################
## Batch 1, subset 1: Now, extract 120 sample subset ##
#######################################################
# This extract should be based on the proportion of all.wanted.df.1200$selection_group
# Define the counts for each group
# Define the counts for each group
counts <- table(all.wanted.df.1200$selection_group)
groups <- names(counts)

# Calculate the proportion of each group
proportions <- counts / sum(counts)

# Sample 120 rows based on the proportions
set.seed(54321) # Set seed for reproducibility
sampled_rows <- lapply(groups, function(group) {
  n <- round(proportions[group] * 120)
  sample(which(all.wanted.df.1200$selection_group == group), n)
})

# Combine the sampled rows
sampled_rows <- unlist(sampled_rows)

# Extract the sampled rows
sampled_data.120 <- all.wanted.df.1200[sampled_rows, ]


table_2.keep <- table_2[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status")]
## add community control for the record
table_2.keep <- rbind.data.frame(table_2.keep, first_ageatsample.zhaoming.100[c("tb_number", "sjlid",  "num_vials", "ageevent", "Sample_age", "vitalstatus", "diaggrp", "CMP_status")])
table_2.keep$CMP_grade <- table_2$grade[match(table_2.keep$tb_number, table_2$tb_number)]
table_2.keep$new_event_number <- table_2$new_event_number[match(table_2.keep$tb_number, table_2$tb_number)]

table_2.keep$Batch1.1200.selecion <- all.wanted.df.1200$tb_number[match(table_2.keep$tb_number, all.wanted.df.1200$tb_number)]
table_2.keep$Batch1.1200.selecion_group <- all.wanted.df.1200$selection_group[match(table_2.keep$tb_number, all.wanted.df.1200$tb_number)]

table_2.keep$Batch1.1200.subset.120.selecion <- sampled_data.120$tb_number[match(table_2.keep$tb_number, sampled_data.120$tb_number)]
table_2.keep$Batch1.1200.subset.120.selecion_group <- sampled_data.120$selection_group[match(table_2.keep$tb_number, sampled_data.120$tb_number)]

write.table(table_2.keep, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v10_output/plasma_data_complete_list_of_table_2_v10.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = "")
write.table(all.wanted.df.1200, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v10_output/plasma_data_batch1_1200_samples.txt", col.names = T, row.names = F, quote = F, sep = "\t", na="")
write.table(sampled_data.120, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v10_output/plasma_data_batch1_1200_samples_subset1_120_samples.txt", col.names = T, row.names = F, quote = F, sep = "\t", na="")

# Note from Yadav: I think you should only provide the necessary information
# when you send these files. You would only need tb_number, sjlid, num_vials and
# the group to ECC people. When you send to the proteomics core, you should only
# send them tb_number and the group so that they can include all 4 groups of
# survivors in each experiment.



# Yadav on 3/26/2024: Achal, can you please provide sex, age at sample and race of these samples?
demographic <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
all.wanted.df.1200.to.update <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v10_output/plasma_data_batch1_1200_samples.txt", header = T, sep = "\t")
all.wanted.df.1200.to.update$Sex <- demographic$gender[match(all.wanted.df.1200.to.update$sjlid, demographic$sjlid)]
all.wanted.df.1200.to.update$racegrp <- demographic$racegrp[match(all.wanted.df.1200.to.update$sjlid, demographic$sjlid)]

write.table(all.wanted.df.1200.to.update, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v10_output/plasma_data_batch1_1200_samples_to_proteomics_core.txt", col.names = T, row.names = F, sep = "\t", quote = F)
table(all.wanted.df.1200.to.update$num_vials, all.wanted.df.1200.to.update$selection_group)






## Extra Checks
cases <- all.wanted.df.1200.to.update[all.wanted.df.1200.to.update$selection_group == "171_CMP_cases" & all.wanted.df.1200.to.update$num_vials == 1 &all.wanted.df.1200.to.update$vitalstatus == "Deceased",]
cbind.data.frame(sjlid=cases$sjlid, diaggrp = cases$diaggrp, racegrp = cases$racegrp)

cases.num.vials.1 <- all.wanted.df.1200.to.update[all.wanted.df.1200.to.update$selection_group == "171_CMP_cases" & all.wanted.df.1200.to.update$num_vials == 1,]
cases.num.vials.1 <- cbind.data.frame(sjlid=cases.num.vials.1$sjlid, diaggrp = cases.num.vials.1$diaggrp, racegrp = cases.num.vials.1$racegrp)

cases.num.vials.1 <- table_2[table_2$sjlid %in% cases.num.vials.1$sjlid,]

Num_samples_before_CMP <- cases.num.vials.1 %>%
  filter(grade_2_or_higher == "grade_0") %>%
  group_by(sjlid) %>%
  summarise(grade_2_or_higher_count = n())
cases.num.vials.1$Num_samples_before_CMP <- Num_samples_before_CMP$grade_2_or_higher_count[match(cases.num.vials.1$sjlid, Num_samples_before_CMP$sjlid)]
## select vial 1
cases.num.vials.1 <- cases.num.vials.1[which(cases.num.vials.1$num_vials == 1),]
cases.num.vials.1$TB_to_exclude[cases.num.vials.1$Num_samples_before_CMP >=2] <- cases.num.vials.1$tb_number[cases.num.vials.1$Num_samples_before_CMP >=2]
# cases.num.vials.1 <- cases.num.vials.1[match(cc$V1, cases.num.vials.1$sjlid), ] ## cc is the list of 41 samples to order
cases.num.vials.1 <- cases.num.vials.1[c("sjlid", "diaggrp", "racegrp", "vitalstatus", "Num_samples_before_CMP","TB_to_exclude", "Sample_age")]
TB <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/ForAchal_Survivors.txt", header = T, sep = "\t")
df <- TB
df$ageatsample <- floor(df$ageatsample * 10) / 10
df$KEY <- paste0(df$sjlid, ":", df$ageatsample)
df.plasma <- df[df$aliquot_type == "Plasma",]
df.serum <- df[df$aliquot_type == "Serum",]
cases.num.vials.1$KEY <- paste0(cases.num.vials.1$sjlid, ":",cases.num.vials.1$Sample_age)

cases.num.vials.1$Plasma_vials <- df.plasma$num_vials[match(cases.num.vials.1$KEY, df.plasma$KEY)]
cases.num.vials.1$Serum_vials <- df.serum$num_vials[match(cases.num.vials.1$KEY, df.serum$KEY)]
cases.num.vials.1.remove <- na.omit(cases.num.vials.1$TB_to_exclude)

## controls
controls <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/plasma_table_2_v9.txt", header = T, sep = "\t") # after removing num vial zero and with 18 or older
remove <- unique(all.wanted.df.1200.to.update$sjlid[all.wanted.df.1200.to.update$selection_group == "171_CMP_cases"])
controls <- controls [!controls$sjlid %in% remove,]
controls$Num_samples <- ave(controls$sjlid, controls$sjlid, FUN = length)
controls.num.vial.1 <- controls[controls$num_vials == 1 & controls$Num_samples ==1,]

controls.num.vial.1$KEY <- paste0(controls.num.vial.1$sjlid, ":",controls.num.vial.1$Sample_age)
controls.num.vial.1$Plasma_vials <- df.plasma$num_vials[match(controls.num.vial.1$KEY, df.plasma$KEY)]
controls.num.vial.1$Serum_vials <- df.serum$num_vials[match(controls.num.vial.1$KEY, df.serum$KEY)]
## exclude last vial available 
controls.num.vial.1.last.vial <- controls.num.vial.1[(controls.num.vial.1$Plasma_vials == 1 & controls.num.vial.1$Serum_vials == 0)|(controls.num.vial.1$Plasma_vials == 1 & is.na(controls.num.vial.1$Serum_vials)),]
controls.num.vial.1.last.vial.to.remove <- controls.num.vial.1.last.vial$tb_number

## exclude control TB if subsequent TBs are available
controls.to.remove <- controls[which(controls$num_vials == 1),]
controls.to.remove$TB_to_exclude[controls.to.remove$Num_samples >=2] <- controls.to.remove$tb_number[controls.to.remove$Num_samples >=2]
controls.to.remove <- c(controls.to.remove$TB_to_exclude[controls.to.remove$Num_samples >=2], controls.num.vial.1.last.vial.to.remove)
all.to.remove <- c(cases.num.vials.1.remove, controls.to.remove)
write.table(as.data.frame(all.to.remove), "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v10_output/all.TB.to.removefrom_V10.txt", sep = "\t", col.names = F)
# ## Test
# library("blockrand")
# randomized_samples <- block_ra(sample_data, n = 14, id_col = "Sample", block_col = "selection_group", strata_cols = c("Sample_age", "Sex", "racegrp"))
# 
# # View the randomized samples
# randomized_samples