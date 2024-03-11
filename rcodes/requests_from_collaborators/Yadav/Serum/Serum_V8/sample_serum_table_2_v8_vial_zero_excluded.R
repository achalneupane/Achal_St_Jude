table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_v8.txt", header = T, sep = "\t") # after removing num vial zero and with 18 or older
## Add primary diagnosis
demog <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/demog.sas7bdat")
table_2$primdx <- demog$diag[match(table_2$mrn, demog$MRN)]

# ## Add TB ID
# TBID <- read.table("All_serum_samples_for_R01_22Aug2023_for_Achal_edited.txt", header = T, sep = "\t")
# TBID <- TBID[TBID$sjlid %in% table_2$sjlid,]
# table_2$TBID_YN  <- ifelse(table_2$sjlid %in% TBID$sjlid, "Yes", "No")

table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number == 2] <- "CA91"
table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number >=3] <- "CA38"

CA91.samples <-  table_2$sjlid[which(table_2$CAcount == "CA91")]
CA38.samples <- table_2$sjlid[which(table_2$CAcount == "CA38")]

#####################
## Get 600 samples ##
#####################
visit2.grade.0 <- table_2[table_2$grade == 0 & table_2$new_event_number >= 2,]
nrow(visit2.grade.0)
# 1558 # v6
# 1521 # v8
# 1033 After removing vial zero
## exclude 100 and 67 from visit2.grade.0
to.exclude <- c(CA91.samples, CA38.samples)
length(to.exclude)
# 177
# 167
# 129
sum(visit2.grade.0$sjlid %in% to.exclude)
# 82
# 81
# 46
visit2.grade.0 <- visit2.grade.0 [!visit2.grade.0$sjlid %in% to.exclude,]

## Get 600-67=532 sjlids from this randomly
## Get 600-38=532 sjlids from this randomly After removin vial zero

set.seed(54321)
random_samples.600 <- sample(unique(visit2.grade.0$sjlid), 562)
## add 67 to this:
## add 38 to this: # after removing vial zero
random_samples.600 <- c(random_samples.600, CA38.samples)
length(random_samples.600)
# 600
#####################
## Get 800 samples ##
#####################
## Now remove 600 and 100 from 2437 and get 800-100=700 samples
## Now remove 600 and 100 from 2228 and get 800-91=709 samples # after removing vial zero
CO.2228 <- table_2$sjlid [table_2$grade==0 & table_2$new_event_number ==1]
length(CO.2228)
# 2799
# 2437
# 2228 # after removing vial zero

## remove 600 samples from this
CO.2228 <- CO.2228[!CO.2228 %in% random_samples.600]
## also remove 100 samples from this
## also remove 91 samples from this # after removing vial zero
CO.2228 <- CO.2228[!CO.2228 %in% CA91.samples]

## Now randomly select 800-91=709 samples from this
random_samples.800 <- sample(unique(CO.2228), 709)
## add 100 to this
random_samples.800 <- c(random_samples.800, CA91.samples)

extract.1400 <- c(random_samples.600, random_samples.800)
length(extract.1400)
# 1400

extract.1400.df <- as.data.frame(extract.1400)
extract.1400.df$extracted_group <- ifelse(extract.1400.df$extract.1400 %in% random_samples.600, "600.with.grade0.at.first2visits", "800.with.grade0.at.first.visit")

## classify CA and CO
extract.1400.df$status <- ifelse(extract.1400.df$extract.1400 %in% c(CA38.samples, CA91.samples), "CA", "CO")


# write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples_v8.txt", col.names = T, row.names = F, sep = "\t")


# Yadav on 02/15/2024: Can you please add age at serum samples at baseline for
# all the 1400 samples, along with age at first follow-up after baseline for the
# 600 samples?
# extract.1400.df <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples.txt", sep = "\t", header = T)

get.baseline <- table_2[table_2$new_event_number == 1,]
age.at.first.followup <- table_2[table_2$new_event_number == 2,]

extract.1400.df$serum.sample.age.at.baseline <- get.baseline$Sample_age[match(extract.1400.df$extract.1400, get.baseline$sjlid)]
extract.1400.df$age.at.first.follow.up <- age.at.first.followup$Sample_age[match(extract.1400.df$extract.1400, age.at.first.followup$sjlid)]

extract.1400.df$age.at.first.follow.up[extract.1400.df$extracted_group == "800.with.grade0.at.first.visit"] <- NA
# write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples_v8_vial_zero_excluded.txt", col.names = T, row.names = F, sep = "\t")

################################
## Submitted; round 1; N=1100 ##
################################
cases.167 <- extract.1400.df[extract.1400.df$status == "CA",]

table_2.hodgkin <- table_2[grepl("Hodgkin", table_2$primdx, ignore.case = T),]
