# table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2.txt", header = T, sep = "\t")
table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_v8.txt", header = T, sep = "\t")


table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number == 2] <- "CA100"
table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number >=3] <- "CA67"

CA100.samples <-  table_2$sjlid[which(table_2$CAcount == "CA100")]
CA67.samples <- table_2$sjlid[which(table_2$CAcount == "CA67")]

#####################
## Get 600 samples ##
#####################
visit2.grade.0 <- table_2[table_2$grade == 0 & table_2$new_event_number >= 2,]
nrow(visit2.grade.0)
# 1558 # v6
# 1521 # v8
## exclude 100 and 67 from visit2.grade.0
to.exclude <- c(CA100.samples, CA67.samples)
length(to.exclude)
# 177
# 167
sum(visit2.grade.0$sjlid %in% to.exclude)
# 82
# 81
visit2.grade.0 <- visit2.grade.0 [!visit2.grade.0$sjlid %in% to.exclude,]

## Get 600-67=532 sjlids from this randomly

set.seed(54321)
random_samples.600 <- sample(unique(visit2.grade.0$sjlid), 533)
## add 67 to this:
random_samples.600 <- c(random_samples.600, CA67.samples)
length(random_samples.600)
# 600
#####################
## Get 800 samples ##
#####################
## Now remove 600 and 100 from 2737 and get 800-100=700 samples
CO.2437 <- table_2$sjlid [table_2$grade==0 & table_2$new_event_number ==1]
length(CO.2437)
# 2799
# 2437

## remove 600 samples from this
CO.2437 <- CO.2437[!CO.2437 %in% random_samples.600]
## also remove 100 samples from this
CO.2437 <- CO.2437[!CO.2437 %in% CA100.samples]

## Now randomly select 800-100=700 samples from this
random_samples.800 <- sample(unique(CO.2437), 700)
## add 100 to this
random_samples.800 <- c(random_samples.800, CA100.samples)

extract.1400 <- c(random_samples.600, random_samples.800)
length(extract.1400)
# 1400

extract.1400.df <- as.data.frame(extract.1400)
extract.1400.df$extracted_group <- ifelse(extract.1400.df$extract.1400 %in% random_samples.600, "600.with.grade0.at.first2visits", "800.with.grade0.at.first.visit")

## classify CA and CO
extract.1400.df$status <- ifelse(extract.1400.df$extract.1400 %in% c(CA67.samples, CA100.samples), "CA", "CO")


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
write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples_v8.txt", col.names = T, row.names = F, sep = "\t")
