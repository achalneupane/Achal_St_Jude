table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2.txt", header = T, sep = "\t")

table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number == 2] <- "CA109"
table_2$CAcount[table_2$grade >= 2 & table_2$new_event_number >=3] <- "CA68"

CA109.samples <-  table_2$sjlid[which(table_2$CAcount == "CA109")]
CA68.samples <- table_2$sjlid[which(table_2$CAcount == "CA68")]

#####################
## Get 600 samples ##
#####################
visit2.grade.0 <- table_2[table_2$grade == 0 & table_2$new_event_number >= 2,]
nrow(visit2.grade.0)
# 1558

## exclude 109 and 68 from visit2.grade.0
to.exclude <- c(CA109.samples, CA68.samples)
length(to.exclude)
# 177
sum(visit2.grade.0$sjlid %in% to.exclude)
# 82
visit2.grade.0 <- visit2.grade.0 [!visit2.grade.0$sjlid %in% to.exclude,]

## Get 600-68=532 sjlids from this randomly

set.seed(54321)
random_samples.600 <- sample(unique(visit2.grade.0$sjlid), 532)
## add 68 to this:
random_samples.600 <- c(random_samples.600, CA68.samples)
length(random_samples.600)
# 600
#####################
## Get 800 samples ##
#####################
## Now remove 600 and 109 from 2799 and get 800-109=691 samples
CO.2799 <- table_2$sjlid [table_2$grade==0 & table_2$new_event_number ==1]
length(CO.2799)
# 2799

## remove 600 samples from this
CO.2799 <- CO.2799[!CO.2799 %in% random_samples.600]
## also remove 109 samples from this
CO.2799 <- CO.2799[!CO.2799 %in% CA109.samples]

## Now randomly select 800-109=691 samples from this
random_samples.800 <- sample(unique(CO.2799), 691)
## add 109 to this
random_samples.800 <- c(random_samples.800, CA109.samples)

extract.1400 <- c(random_samples.600, random_samples.800)
length(extract.1400)

extract.1400.df <- as.data.frame(extract.1400)
extract.1400.df$extracted_group <- ifelse(extract.1400.df$extract.1400 %in% random_samples.600, "600.with.grade0.at.first2visits", "800.with.grade0.at.first.visit")

## classify CA and CO
extract.1400.df$status <- ifelse(extract.1400.df$extract.1400 %in% c(CA68.samples, CA109.samples), "CA", "CO")
write.table(extract.1400.df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_extracted_batch1_1400_samples.txt", col.names = T, row.names = F, sep = "\t")
