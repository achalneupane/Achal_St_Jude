# table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2.txt", header = T, sep = "\t") # all age and vial zero included
# table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_v8_with_vial_zero.txt", header = T, sep = "\t") # 18 or older only but vial zero included

# table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_v8.txt", header = T, sep = "\t") # 18 or older only and vial zero excluded 
table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/plasma_table_2_v8.txt", header = T, sep = "\t") # 18 or older only and vial zero excluded

table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number == 2] <- "CA107"
table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number >=3] <- "CA64"

CA107.samples <-  table_2$sjlid[which(table_2$CAcount == "CA107")]
CA64.samples <- table_2$sjlid[which(table_2$CAcount == "CA64")]

#####################
## Get 600 samples ##
#####################
visit2.grade.0 <- table_2[table_2$grade == 0 & table_2$new_event_number >= 2,]
nrow(visit2.grade.0)
# 1503 # plasma
## exclude 107 and 64 from visit2.grade.0
to.exclude <- c(CA107.samples, CA64.samples)
length(to.exclude)
# 171
sum(visit2.grade.0$sjlid %in% to.exclude)
# 77

visit2.grade.0 <- visit2.grade.0 [!visit2.grade.0$sjlid %in% to.exclude,]

## Get 600-64=537 sjlids from this randomly

set.seed(54321)
random_samples.600 <- sample(unique(visit2.grade.0$sjlid), 536)
## add 63 to this:
random_samples.600 <- c(random_samples.600, CA64.samples)
length(random_samples.600)
# 600
#####################
## Get 800 samples ##
#####################
## Now remove 600 and 107 from 2432 and get 800-107=693 samples
CO.2432 <- table_2$sjlid [table_2$grade==0 & table_2$new_event_number ==1]
length(CO.2432)
# 2799
# 2432

## remove 600 samples from this
CO.2432 <- CO.2432[!CO.2432 %in% random_samples.600]
## also remove 100 samples from this
CO.2432 <- CO.2432[!CO.2432 %in% CA107.samples]

## Now randomly select 800-107=693 samples from this
random_samples.800 <- sample(unique(CO.2432), 693)
## add 100 to this
random_samples.800 <- c(random_samples.800, CA107.samples)

extract.1400 <- c(random_samples.600, random_samples.800)
length(extract.1400)
# 1400

extract.1400.df <- as.data.frame(extract.1400)
extract.1400.df$extracted_group <- ifelse(extract.1400.df$extract.1400 %in% random_samples.600, "600.with.grade0.at.first2visits", "800.with.grade0.at.first.visit")

## classify CA and CO
extract.1400.df$status <- ifelse(extract.1400.df$extract.1400 %in% c(CA64.samples, CA107.samples), "CA", "CO")


# write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples_v8.txt", col.names = T, row.names = F, sep = "\t")
# write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/plasma_table_2_extracted_batch1_1400_samples_v8.txt", col.names = T, row.names = F, sep = "\t")

# Yadav on 02/15/2024: Can you please add age at serum samples at baseline for
# all the 1400 samples, along with age at first follow-up after baseline for the
# 600 samples?
# extract.1400.df <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples.txt", sep = "\t", header = T)

get.baseline <- table_2[table_2$new_event_number == 1,]
age.at.first.followup <- table_2[table_2$new_event_number == 2,]

extract.1400.df$serum.sample.age.at.baseline <- get.baseline$Sample_age[match(extract.1400.df$extract.1400, get.baseline$sjlid)]
extract.1400.df$age.at.first.follow.up <- age.at.first.followup$Sample_age[match(extract.1400.df$extract.1400, age.at.first.followup$sjlid)]

extract.1400.df$age.at.first.follow.up[extract.1400.df$extracted_group == "800.with.grade0.at.first.visit"] <- NA
# write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples_v8.txt", col.names = T, row.names = F, sep = "\t")
write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/plasma_table_2_extracted_batch1_1400_samples_v8.txt", col.names = T, row.names = F, sep = "\t")
