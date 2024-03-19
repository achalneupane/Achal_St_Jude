library(haven)
library(dplyr)
table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/plasma_table_2_v9.txt", header = T, sep = "\t") # after removing num vial zero and with 18 or older
#########################
## Zhaoming's controls ##
#########################
zhaoming.control <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/controls.selected.110.txt", header = F, sep = "\t")
# zhaoming.control$V1 %in% SERUM$sjlid
zhaoming.control.from.mathew <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/data_from_Matthew/forachal_controls.txt", header = T, sep = "\t")
zhaoming.control.from.mathew <- zhaoming.control.from.mathew[grepl("Plasma", zhaoming.control.from.mathew$aliquot_type, ignore.case = T),]
## remove 10 samples with the lowest number of vials
grouped_data <- zhaoming.control.from.mathew %>%
  group_by(sjlid) %>%
  summarize(num_vials = min(num_vials))

# Find the 10 sjlid with the lowest number of vials
sorted_data <- zhaoming.control.from.mathew %>%
  arrange(num_vials)
remove.10.samples <- head(unique(sorted_data$sjlid),10)

zhaoming.control.from.mathew <- zhaoming.control.from.mathew[!zhaoming.control.from.mathew$sjlid %in% remove.10.samples,]
length(unique(zhaoming.control.from.mathew$sjlid))
# 100
zhaoming.controls.100 <- unique(zhaoming.control.from.mathew$sjlid)
###########################
## Add primary diagnosis ##
###########################
demog <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/demog.sas7bdat")
table_2$primdx <- demog$diag[match(table_2$mrn, demog$MRN)]

###################
## Get 171 cases ##
###################
CA.171 <- table_2[table_2$grade_2_or_higher =="grade_2_or_higher",]
CA.171 <- unique(CA.171$sjlid)

#################################################################################################
## randomly select 200 Hodgkin lymphoma survivors from all eligible Hodgkin lymphoma survivors ##
#################################################################################################
table_3 <- table_2[!table_2$sjlid %in% CA.171,] # exclude those in CA.171 samples from the original table
all.hodgkin <- table_3[grepl("\\bHodgkin's\\b", table_3$primdx, ignore.case = T),]
all.hodgkin <- all.hodgkin[!grepl("non", all.hodgkin$primdx, ignore.case = T),]

set.seed(54321)
# select 200
all.hodgkin.extract <- sample(unique(all.hodgkin$sjlid), 200)

####################################################
## randomly select the remaining 1100-171-200=733 ##
####################################################
table_3 <- table_2[!table_2$sjlid %in% c(CA.171,all.hodgkin.extract),] # exclude those in CA.171 + 200 hodgkins samples from the original table
# Randomly select 733 non- Hodgkin lymphoma survivors from all eligible non- Hodgkin lymphoma survivors.
all.non.hodgkin <- table_3[!table_3$sjlid %in% all.hodgkin$sjlid,]

set.seed(54321)
# select 733
all.non.hodgkin.extract <- sample(unique(all.non.hodgkin$sjlid), 729)

#################################
## concatenate all wanted ones ##
#################################
all.wanted.sjlids <- c(CA.171, all.hodgkin.extract, all.non.hodgkin.extract, zhaoming.controls.100)
length(all.wanted.sjlids)

#########################
## Now, extract TB IDs ##
#########################
extract.1200.samples <- table_2[table_2$sjlid %in% all.wanted.sjlids,]
