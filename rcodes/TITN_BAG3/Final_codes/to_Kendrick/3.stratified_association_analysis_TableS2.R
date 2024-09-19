rm(list=ls())

## Perform stratified association analysis with top TTN and BAG3 common missense variants

load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/EUR_common_variants.RData")
dat_final <- EUR_common_variants
dat_final$CMP <- factor(ifelse(dat_final$CMP == 2, 1, 0))


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
  gender <- df$Gender
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
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250 & dat_final$hrtavg ==0)| 
                  (dat_final$hrtavg >=30 & dat_final$anthra_jco_dose_any ==0)| 
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
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250 & dat_final$hrtavg ==0)| 
                        (dat_final$hrtavg >=30 & dat_final$anthra_jco_dose_any ==0)| 
                        (dat_final$anthra_jco_dose_any >=100 &  dat_final$hrtavg >= 15),]


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
objects_list <- list(ttn.all, ttn.male, ttn.female, ttn.high.risk.trt, ttn.moderate.risk.trt, ttn.low.risk.trt, bag3.all, bag3.male, bag3.female, bag3.high.risk.trt, bag3.low.risk.trt)

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
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/AFR_common_variants.RData")
dat_final = AFR_common_variants
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
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250 & dat_final$hrtavg ==0)| 
                        (dat_final$hrtavg >=30 & dat_final$anthra_jco_dose_any ==0)| 
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
all_data <- dat_final[(dat_final$anthra_jco_dose_any >=250 & dat_final$hrtavg ==0)| 
                        (dat_final$hrtavg >=30 & dat_final$anthra_jco_dose_any ==0)| 
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
