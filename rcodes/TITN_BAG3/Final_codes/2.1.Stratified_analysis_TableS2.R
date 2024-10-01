rm(list=ls())

## Perform stratified association analysis with top TTN and BAG3 common missense variants

## From Kendrick
# sjlife_cmp_data <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Re_A_new_manuscript_on_TTN_BAG3_for_your_review//sjlife_ccm_analysis_dat.rds") # N = 3686
# ccss_cmp_data <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Re_A_new_manuscript_on_TTN_BAG3_for_your_review/ccss_ccm_analysis_dat.rds") 
EUR_common_Kendrick <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis//EUR_CMP_analysis_data.rds")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/EUR_common_variants.RData")
## EUR
EUR_common_Kendrick$CMP <- EUR_common_variants$CMP[match(EUR_common_Kendrick$iid, EUR_common_variants$IID)]
EUR_common_Kendrick$CMP[EUR_common_Kendrick$CMP==1] <- 0; EUR_common_Kendrick$CMP[EUR_common_Kendrick$CMP==2] <- 1
EUR_common_Kendrick$anthra_jco_dose_any <- EUR_common_Kendrick$anthracyclines_dose_5
EUR_common_Kendrick$hrtavg <- EUR_common_Kendrick$heartavg
EUR_common_Kendrick$hrtavg <- EUR_common_Kendrick$hrtavg/100
colnames(EUR_common_Kendrick) <- gsub("^pc", "PC", colnames(EUR_common_Kendrick), ignore.case = TRUE)
EUR_common_Kendrick <- EUR_common_Kendrick[!(is.na(EUR_common_Kendrick$anthracyclines_dose_5)|is.na(EUR_common_Kendrick$heartavg)),]
colnames(EUR_common_Kendrick)[colnames(EUR_common_Kendrick) == "fid"] <- "FID"
colnames(EUR_common_Kendrick)[colnames(EUR_common_Kendrick) == "iid"] <- "IID"
EUR_common_Kendrick.pheno <- EUR_common_Kendrick[,c("FID", "IID", "CMP", "agelstcontact", "gender", "agevent", "ejection_fraction_hrt", "PC1", "PC2",
                                                   "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "cohort", "cohort_two", "agedx", "anthra_jco_dose_any", "hrtavg")]

EUR_common_Kendrick.pheno$CMP[EUR_common_Kendrick.pheno$CMP==1] <- 2; EUR_common_Kendrick.pheno$CMP[EUR_common_Kendrick.pheno$CMP==0] <- 1
table(EUR_common_Kendrick.pheno$CMP)
# 1    2 
# 5729  453 
# write.table(EUR_common_Kendrick.pheno, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", sep ="\t", col.names = T, row.names = F, quote = F)

AFR_common_Kendrick <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis//AFR_CMP_analysis_data.rds")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/AFR_common_variants.RData")
## AFR
AFR_common_Kendrick$CMP <- AFR_common_variants$CMP[match(AFR_common_Kendrick$iid, AFR_common_variants$IID)]
# AFR_common_Kendrick$CMP[AFR_common_Kendrick$CMP==1] <- 0; AFR_common_Kendrick$CMP[AFR_common_Kendrick$CMP==2] <- 1
AFR_common_Kendrick$anthra_jco_dose_any <- AFR_common_Kendrick$anthracyclines_dose_5
AFR_common_Kendrick$hrtavg <- AFR_common_Kendrick$heartavg
AFR_common_Kendrick$hrtavg <- AFR_common_Kendrick$hrtavg/100
colnames(AFR_common_Kendrick) <- gsub("^pc", "PC", colnames(AFR_common_Kendrick), ignore.case = TRUE)
AFR_common_Kendrick <- AFR_common_Kendrick[!(is.na(AFR_common_Kendrick$anthracyclines_dose_5)|is.na(AFR_common_Kendrick$heartavg)),]
colnames(AFR_common_Kendrick)[colnames(AFR_common_Kendrick) == "fid"] <- "FID"
colnames(AFR_common_Kendrick)[colnames(AFR_common_Kendrick) == "iid"] <- "IID"
AFR_common_Kendrick.pheno <- AFR_common_Kendrick[,c("FID", "IID", "CMP", "agelstcontact", "gender", "ejection_fraction_hrt", "PC1", "PC2",
                                                   "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "cohort", "agedx", "anthra_jco_dose_any", "hrtavg")]

AFR_common_Kendrick.pheno$CMP[AFR_common_Kendrick.pheno$CMP==1] <- 2; AFR_common_Kendrick.pheno$CMP[AFR_common_Kendrick.pheno$CMP==0] <- 1
table(AFR_common_Kendrick.pheno$CMP)
# 1   2 
# 201  37
# write.table(AFR_common_Kendrick.pheno, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3_afr_kendrick.pheno", sep ="\t", col.names = T, row.names = F, quote = F)



# EUR_common_variants <- EUR_common_variants[EUR_common_variants$IID %in% EUR_common_Kendrick$iid,]
# EUR_common_variants$CMP <- factor(ifelse(EUR_common_variants$CMP == 2, 1, 0))
# EUR_common_variants$CMP_kendrick <- EUR_common_Kendrick$pre_diag_5_cmp[match(EUR_common_variants$IID, EUR_common_Kendrick$iid)]
# table(Our=EUR_common_variants$CMP, Kendrick=EUR_common_variants$CMP_kendrick)
# EUR_common_variants$anthra_jco_dose_any.kendrick <- EUR_common_Kendrick$anthracyclines_dose_5[match(EUR_common_variants$IID, EUR_common_Kendrick$iid)]
# table(EUR_common_variants$anthra_jco_dose_any == EUR_common_variants$anthra_jco_dose_any.kendrick)
# EUR_common_variants$hrtavg.kendrick <- EUR_common_Kendrick$heartavg[match(EUR_common_variants$IID, EUR_common_Kendrick$iid)]
# table(EUR_common_variants$hrtavg==EUR_common_variants$hrtavg.kendrick)
# EUR_common_variants$hrtavg <- EUR_common_variants$hrtavg.kendrick
# EUR_common_variants$anthra_jco_dose_any <- EUR_common_variants$anthra_jco_dose_any.kendrick
# EUR_common_variants$agelstcontact <- EUR_common_variants$agelstcontact
# # EUR_common_variants$CMP <- EUR_common_variants$CMP_kendrick

dat_final <- EUR_common_Kendrick
# dat_final$CMP <- factor(ifelse(dat_final$CMP == 2, 1, 0))


# Function to calculate odds ratio and confidence interval
calculate_odds_ratio <- function(coef, se) {
  odds_ratio <- exp(coef)
  ci_lower <- exp(coef - 1.96 * se)
  ci_upper <- exp(coef + 1.96 * se)
  return(list(odds_ratio = odds_ratio, ci_lower = ci_lower, ci_upper = ci_upper))
}


# Function to combine OR and CI and format decimal places
combine_or_ci <- function(estimate, ci_lower, ci_upper) {
  combined <- paste0(round(estimate, 3), " (", round(ci_lower, 3), "-", round(ci_upper, 3), ")")
  return(combined)
}

# Extract Gender, P_Value, and OR_CI
extract_info <- function(df) {
  gender <- df$Group
  p_value <- df$P_Value
  or_ci <- df$OR_CI
  N = df$n
  return(data.frame(Gender = gender, P_Value = p_value, OR_CI = or_ci, n = N))
}

#########
## TTN ##
#########
#------------------1. ALL



# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c("All"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)

# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.all <- table_data



# Overall Analysis (carrier and gender in the same model):


# Based on https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(23)00012-8/fulltext

# High Risk: Patients who have received a cumulative anthracycline dose of ≥250 mg/m² and/or chest-directed radiotherapy of ≥30 Gy, or a combination of ≥100 mg/m² anthracycline and ≥15 Gy chest-directed radiotherapy.
# Moderate Risk: Patients who have received a cumulative anthracycline dose between 100 and <250 mg/m² and/or chest-directed radiotherapy between 15 and <30 Gy. The table has "NA" (Not Applicable) for the combination of anthracycline and radiotherapy in this risk group.
# Low Risk: Patients who have received a cumulative anthracycline dose >0 and <100 mg/m² and/or chest-directed radiotherapy >0 and <15 Gy.

#------------------2. High risk treatment
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250)| 
                        (dat_final$hrtavg >=30)| 
                        (dat_final$anthra_jco_dose_any >=100 &  dat_final$hrtavg >= 15),]


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary


# Create the table
table_data <- data.frame(
  Group = c("High-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.high.risk.trt <- table_data


#------------------3. moderate treatment risk
all_data.moderate.risk <- dat_final[
  ((dat_final$anthra_jco_dose_any >= 100 & dat_final$anthra_jco_dose_any < 250 & dat_final$hrtavg == 0)| 
     (dat_final$hrtavg >= 15 & dat_final$hrtavg < 30 & dat_final$anthra_jco_dose_any == 0)),
]



all_data <- all_data.moderate.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Moderate-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.moderate.risk.trt <- table_data

#------------------3. Low treatment risk

all_data.low.risk <- dat_final[
  (dat_final$anthra_jco_dose_any > 0 & dat_final$anthra_jco_dose_any < 100 & 
     dat_final$hrtavg ==0)|
    (dat_final$hrtavg > 0 & dat_final$hrtavg < 15 & dat_final$anthra_jco_dose_any == 0),
]

all_data <- all_data.low.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Low-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.low.risk.trt <- table_data


#------------------4. Male
all_data <- dat_final[dat_final$gender == 0,]


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary


# Create the table
table_data <- data.frame(
  Group = c("Male"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.male <- table_data


#------------------5. Female
all_data <- dat_final[dat_final$gender == 1,]


overall_model <- glm(CMP ~ rs3829746 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary


# Create the table
table_data <- data.frame(
  Group = c("Female"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.female <- table_data
##########
## BAG3 ##
##########
#------------------1. All

# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c("All"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.all <- table_data


#------------------2. High risk treatment
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250)| 
                        (dat_final$hrtavg >=30)| 
                        (dat_final$anthra_jco_dose_any >=100 &  dat_final$hrtavg >= 15),]

cc <- cbind.data.frame(dat_final$iid, dat_final$hrtavg, dat_final$anthra_jco_dose_any, dat_final$risk_group)


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("High-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.high.risk.trt <- table_data


#------------------3. moderate treatment risk
all_data.moderate.risk <- dat_final[
  ((dat_final$anthra_jco_dose_any >= 100 & dat_final$anthra_jco_dose_any < 250 & dat_final$hrtavg == 0)| 
     (dat_final$hrtavg >= 15 & dat_final$hrtavg < 30 & dat_final$anthra_jco_dose_any == 0)),
]



all_data <- all_data.moderate.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Moderate-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.moderate.risk.trt <- table_data



#------------------3. Low treatment risk

all_data.low.risk <- dat_final[
  (dat_final$anthra_jco_dose_any > 0 & dat_final$anthra_jco_dose_any < 100 & 
     dat_final$hrtavg ==0)|
    (dat_final$hrtavg > 0 & dat_final$hrtavg < 15 & dat_final$anthra_jco_dose_any == 0),
]

all_data <- all_data.low.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Low-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.low.risk.trt <- table_data



#------------------4. Male
all_data <- dat_final[dat_final$gender == 0,]


overall_model <- glm(CMP ~ rs2234962 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Male"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.male <- table_data


#------------------5. Female
all_data <- dat_final[dat_final$gender == 1,]


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Female"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.female <- table_data


# Create a list of your objects
objects_list <- list(ttn.all, ttn.male, ttn.female, ttn.high.risk.trt, ttn.moderate.risk.trt, ttn.low.risk.trt, bag3.all, bag3.male, bag3.female, bag3.high.risk.trt, bag3.moderate.risk.trt, bag3.low.risk.trt)

# Extract information and combine into a data frame
combined_info <- do.call(rbind, lapply(objects_list, extract_info))

# Print the combined information
print(combined_info)


######################
## African ancestry ##
######################
######################
## African ancestry ##
######################
######################
## African ancestry ##
######################
######################
## African ancestry ##
######################
dat_final <- AFR_common_Kendrick
table(dat_final$CMP)


#########
## TTN ##
#########
#------------------1. ALL

# Overall Analysis
all_data <- dat_final

overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c("All"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)

# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.all <- table_data



# Overall Analysis (carrier and gender in the same model):


# Based on https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(23)00012-8/fulltext

# High Risk: Patients who have received a cumulative anthracycline dose of ≥250 mg/m² and/or chest-directed radiotherapy of ≥30 Gy, or a combination of ≥100 mg/m² anthracycline and ≥15 Gy chest-directed radiotherapy.
# Moderate Risk: Patients who have received a cumulative anthracycline dose between 100 and <250 mg/m² and/or chest-directed radiotherapy between 15 and <30 Gy. The table has "NA" (Not Applicable) for the combination of anthracycline and radiotherapy in this risk group.
# Low Risk: Patients who have received a cumulative anthracycline dose >0 and <100 mg/m² and/or chest-directed radiotherapy >0 and <15 Gy.

#------------------2. High risk treatment
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250)| 
                        (dat_final$hrtavg >=30)| 
                        (dat_final$anthra_jco_dose_any >=100 &  dat_final$hrtavg >= 15),]


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary


# Create the table
table_data <- data.frame(
  Group = c("High-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.high.risk.trt <- table_data


#------------------3. moderate treatment risk
all_data.moderate.risk <- dat_final[
  ((dat_final$anthra_jco_dose_any >= 100 & dat_final$anthra_jco_dose_any < 250 & dat_final$hrtavg == 0)| 
     (dat_final$hrtavg >= 15 & dat_final$hrtavg < 30 & dat_final$anthra_jco_dose_any == 0)),
]



all_data <- all_data.moderate.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Moderate-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.moderate.risk.trt <- table_data

#------------------3. Low treatment risk

all_data.low.risk <- dat_final[
  (dat_final$anthra_jco_dose_any > 0 & dat_final$anthra_jco_dose_any < 100 & 
     dat_final$hrtavg ==0)|
    (dat_final$hrtavg > 0 & dat_final$hrtavg < 15 & dat_final$anthra_jco_dose_any == 0),
]

all_data <- all_data.low.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Low-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.low.risk.trt <- table_data


#------------------4. Male
all_data <- dat_final[dat_final$gender == 0,]


# Overall Analysis
overall_model <- glm(CMP ~ rs3829746 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary


# Create the table
table_data <- data.frame(
  Group = c("Male"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.male <- table_data


#------------------5. Female
all_data <- dat_final[dat_final$gender == 1,]


overall_model <- glm(CMP ~ rs3829746 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary


# Create the table
table_data <- data.frame(
  Group = c("Female"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

ttn.female <- table_data
##########
## BAG3 ##
##########
#------------------1. All

# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c("All"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.all <- table_data


#------------------2. High risk treatment
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250)| 
                        (dat_final$hrtavg >=30)| 
                        (dat_final$anthra_jco_dose_any >=100 &  dat_final$hrtavg >= 15),]


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("High-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.high.risk.trt <- table_data


#------------------3. moderate treatment risk
all_data.moderate.risk <- dat_final[
  ((dat_final$anthra_jco_dose_any >= 100 & dat_final$anthra_jco_dose_any < 250 & dat_final$hrtavg == 0)| 
     (dat_final$hrtavg >= 15 & dat_final$hrtavg < 30 & dat_final$anthra_jco_dose_any == 0)),
]



all_data <- all_data.moderate.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Moderate-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.moderate.risk.trt <- table_data



#------------------3. Low treatment risk

all_data.low.risk <- dat_final[
  (dat_final$anthra_jco_dose_any > 0 & dat_final$anthra_jco_dose_any < 100 & 
     dat_final$hrtavg ==0)|
    (dat_final$hrtavg > 0 & dat_final$hrtavg < 15 & dat_final$anthra_jco_dose_any == 0),
]

all_data <- all_data.low.risk


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Low-risk"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.low.risk.trt <- table_data



#------------------4. Male
all_data <- dat_final[dat_final$gender == 0,]


overall_model <- glm(CMP ~ rs2234962 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Male"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.male <- table_data


#------------------5. Female
all_data <- dat_final[dat_final$gender == 1,]


# Overall Analysis
overall_model <- glm(CMP ~ rs2234962 + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary



# Create the table
table_data <- data.frame(
  Group = c("Female"),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)


# Apply function to data frame
table_data$OR_CI <- combine_or_ci(table_data$Odds_Ratio, table_data$CI_Lower, table_data$CI_Upper)

# Print the updated data frame
print(table_data)

bag3.female <- table_data


# Create a list of your objects
objects_list <- list(ttn.all, ttn.male, ttn.female, ttn.high.risk.trt, ttn.moderate.risk.trt, ttn.low.risk.trt, bag3.all, bag3.male, bag3.female, bag3.high.risk.trt, bag3.low.risk.trt)

# Extract information and combine into a data frame
combined_info <- do.call(rbind, lapply(objects_list, extract_info))

# Print the combined information
print(combined_info)
