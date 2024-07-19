library(dplyr)
library(tidyr)
library(haven)

df <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma_table_2_extracted_batch1_1400_samples_v8_ejection_fraction.txt", header = T, sep = "\t", stringsAsFactors = F)
df.600 <- df[df$extracted_group == "600.with.grade0.at.first2visits",]
dim(df.600)

# gg <- EF.600[grepl("sjlid|Ejec", colnames(EF.600), ignore.case = T)]
EF <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/echo_machine.sas7bdat")

# We have EF data from multiple sources and in that case, we used a hierarchy
# rule using the first non-missing EF values from 4D, 3D, Biplane and 2D (in
# this order).

## ********* This probably need to be revised. start with 3D, biplane, 2D... since 4D is not available
EF <- EF[grepl("sjlid|date|Ejec", colnames(EF), ignore.case = T)]
EF$ejection_fraction <- EF$LV_Ejection_Fraction_3D
EF$ejection_fraction[is.na(EF$ejection_fraction)] <- EF$LV_Ejection_Fraction_MOD_BP[is.na(EF$ejection_fraction)]  
EF$ejection_fraction[is.na(EF$ejection_fraction)] <- EF$LV_Ejection_Fraction_2D_Teich[is.na(EF$ejection_fraction)]  
EF$ejection_fraction[is.na(EF$ejection_fraction)] <- EF$LV_Ejection_Fraction_2D_Cubed[is.na(EF$ejection_fraction)]  
EF$ejection_fraction[is.na(EF$ejection_fraction)] <- EF$LV_Ejection_Fraction_MOD_4C[is.na(EF$ejection_fraction)]
EF$ejection_fraction[is.na(EF$ejection_fraction)] <- EF$LV_Ejection_Fraction_4C_AL[is.na(EF$ejection_fraction)]
EF$ejection_fraction[is.na(EF$ejection_fraction)] <- EF$LV_Ejection_Fraction_MOD_2C[is.na(EF$ejection_fraction)]  
EF$ejection_fraction[is.na(EF$ejection_fraction)] <- EF$LV_Ejection_Fraction_MM_Teich[is.na(EF$ejection_fraction)]  

# EF.res <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/echo_research.sas7bdat")
# EF.600.res <- EF[EF$sjlid %in% df.600$SJLID ,]


EF.600.rest <- EF[EF$sjlid %in% df.600$SJLID ,]
EF.600.rest <- EF.600.rest[grepl("sjlid|date|Date|ejection_fraction", colnames(EF.600.rest), ignore.case = F)]
EF.600.rest <- EF.600.rest[!is.na(EF.600.rest$ejection_fraction),]
pp <- EF.600.rest

table_2 <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/plasma_table_2_v8.txt", header = T, sep = "\t") # 18 or older only and vial zero excluded
table_2.600 <- table_2[table_2$sjlid %in% df.600$SJLID,]

table_2.600$gradedt <- as.Date(table_2.600$gradedt)


# Initialize an empty vector to store the matched ejection_fraction values
table_2.600$ejection_fraction <- NA

# Step 1: Match sjlid with sjlid and handle duplicates
for (i in 1:nrow(table_2.600)) {
  print(paste0("Doing i: ", i))
  # Get the current sjlid and gradedt
  sjlid_current <- table_2.600$sjlid[i]
  gradedt_current <- table_2.600$gradedt[i]
  # sub.pp <- pp[pp$sjlid == sjlid_current,]
  # Find rows in pp that match the current sjlid
  matching_rows <- which(pp$sjlid == sjlid_current)
  
  if (length(matching_rows) > 0) {
    # Step 2: Within duplicated sjlid rows, match gradedt with DateVisitStart
    matching_date <- abs(as.numeric(difftime(pp$studydatetime[matching_rows], gradedt_current, units = "days"))) <= 7
    
    if (any(matching_date)) {
      # Extract the ejection_fraction for the matched row
      table_2.600$ejection_fraction[i] <- pp$ejection_fraction[matching_rows[matching_date]][1]
    } else {
      # Handle cases where no exact date match is found
      table_2.600$ejection_fraction[i] <- NA
    }
  } else {
    # Handle cases where no sjlid match is found
    table_2.600$ejection_fraction[i] <- NA
  }
}

table_2.600$ejection_fraction[which(table_2.600$ejection_fraction > 1)] <- table_2.600$ejection_fraction[which(table_2.600$ejection_fraction > 1)]/100
# cc <- cbind.data.frame(table_2.600$sjlid, table_2.600$ejection_fraction)

table(table_2.600$new_event_number, table_2.600$grade_2_or_higher)

table_2.600 <- table_2.600[table_2.600$grade_2_or_higher == "grade_0",]

table(table_2.600$new_event_number, table_2.600$grade_2_or_higher)

table_2.600.event1 <- table_2.600[table_2.600$new_event_number == 1,]
table_2.600.event2 <- table_2.600[table_2.600$new_event_number == 2,]


table_2.600.event1$ejection_fraction_Baseline <- table_2.600.event1$ejection_fraction
table_2.600.event1$gradedtBase <- table_2.600.event1$gradedt

table_2.600.event1$ejection_fraction_FU <- table_2.600.event2$ejection_fraction[match(table_2.600.event1$sjlid, table_2.600.event2$sjlid)]
table_2.600.event1$gradedtFU <- table_2.600.event2$gradedt[match(table_2.600.event1$sjlid, table_2.600.event2$sjlid)]

table_2.600.event <- table_2.600.event1[c("sjlid",	"gender",	"condition", "gradedtBase", "gradedtFU", "ejection_fraction_Baseline", "ejection_fraction_FU")]


write.table(table_2.600.event, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_plasma_600_ejection_fraction_updated.txt", col.names = T, row.names = F, quote = F, sep = "\t")

table_2.600.event <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/table_2_plasma_600_ejection_fraction_updated.txt", header = T, sep = "\t") # 18 or older only and vial zero excluded

# calculate the change in two EF and provide a summary of the change
table_2.600.event$change_in_ef <- table_2.600.event$ejection_fraction_Baseline-table_2.600.event$ejection_fraction_FU



# Remove rows with NA values in the change_in_ef
cleaned_data <- na.omit(table_2.600.event)

# Summary of the change in ejection fraction
summary_change_ef <- summary(cleaned_data$change_in_ef)

# Print the summary
print(summary_change_ef)


# Perform T-test
data <- cbind.data.frame(ejection_fraction_Baseline=table_2.600.event$ejection_fraction_Baseline, ejection_fraction_FU=table_2.600.event$ejection_fraction_FU)
data_clean <- na.omit(data)

# Perform paired t-test
t_test_result <- t.test(data_clean$ejection_fraction_Baseline, data_clean$ejection_fraction_FU, paired = TRUE)

# Print the results
print(t_test_result)

## Results:

# The analysis of the change in ejection fraction (EF) from baseline to follow-up reveals a range of variations among the individuals. The minimum observed change in EF was a decrease of -0.467, while the maximum observed change was an increase of 0.291. The first quartile indicates that the bottom 25% of the sample experienced a decrease of at least -0.060 in EF. The median change in EF is -0.008, suggesting that half of the sample had a change in EF less than or equal to this value, indicating a slight overall decrease. The mean change in EF across the sample is -0.017, further suggesting a modest overall decline in EF. Finally, the third quartile shows that the top 25% of the sample experienced an increase in EF up to 0.037. These findings highlight the variability in EF changes within the sample, with some individuals experiencing significant decreases, some experiencing increases, and most having modest changes.
# 
# Table 1. Summary of the change in ejection fraction from baseline to follow-up
# Min	1st Quartile	Median	Mean	3rd Quartile	Max
# -0.467	-0.060	-0.008	-0.017	0.037	0.291
# 
# Additionally, the paired t-test conducted to assess the change in ejection fraction from baseline to follow-up revealed a statistically significant difference (P=1.16×10⁻⁵). The mean difference between baseline and follow-up ejection fractions was -0.017 (95% CI=-0.024 to -0.0095), which highlights a meaningful decrease in ejection fraction.

