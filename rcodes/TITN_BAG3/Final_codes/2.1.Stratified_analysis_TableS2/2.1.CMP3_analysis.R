rm(list=ls())
diagDT <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/diagDT.rds")
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


## CMP3 plus
ccss_cmp3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_eur_cardiotoxic_exposed.pheno", header =T)
ccss_cmp3 <- ccss_cmp3[!is.na(ccss_cmp3$CMP3plus),]
ccss_cmp3 <- ccss_cmp3[c("IID", "CMP3plus")]
sjlife_cmp3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3.pheno", header =T)
sjlife_cmp3 <- sjlife_cmp3[sjlife_cmp3$CMP_EF_HEIRARCHY !=2,]
sjlife_cmp3$CMP3plus <- ifelse(sjlife_cmp3$CMP_EF_HEIRARCHY == 0, 1,2)
sjlife_cmp3 <- sjlife_cmp3[c("IID", "CMP3plus")]
sjlife_afr_cmp3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3_afr.pheno", header =T)
sjlife_afr_cmp3 <- sjlife_afr_cmp3[sjlife_afr_cmp3$CMP_EF_HEIRARCHY !=2,]
sjlife_afr_cmp3$CMP3plus <- ifelse(sjlife_afr_cmp3$CMP_EF_HEIRARCHY == 0, 1,2)
sjlife_afr_cmp3 <- sjlife_afr_cmp3[c("IID", "CMP3plus")]

all.cmp3plus <- rbind.data.frame(ccss_cmp3, sjlife_cmp3, sjlife_afr_cmp3)

all.cmp3plus$CMP3plus <- ifelse(all.cmp3plus$CMP3plus == 1, 0,1)

table(EUR_common_Kendrick$IID %in% all.cmp3plus$IID)
# FALSE  TRUE 
# 216  5966 

EUR_common_Kendrick <- EUR_common_Kendrick[EUR_common_Kendrick$IID %in% all.cmp3plus$IID,]
EUR_common_Kendrick$CMP <- all.cmp3plus$CMP3plus[match(EUR_common_Kendrick$IID, all.cmp3plus$IID)]

table(AFR_common_Kendrick$IID %in% all.cmp3plus$IID)
# FALSE  TRUE 
# 19   219 
AFR_common_Kendrick <- AFR_common_Kendrick[AFR_common_Kendrick$IID %in% all.cmp3plus$IID,]
AFR_common_Kendrick$CMP <- all.cmp3plus$CMP3plus[match(AFR_common_Kendrick$IID, all.cmp3plus$IID)]


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
  Estimate=df$Estimate 
  Std_Error=df$Std_Error
  p_value <- df$P_Value
  or_ci <- df$OR_CI
  N = df$n
  return(data.frame(Gender = gender, Estimate=Estimate, Std_Error=Std_Error, P_Value = p_value, OR_CI = or_ci, n = N))
}

EUR_common_Kendrick$era_numeric <- diagDT$era_numeric[match(EUR_common_Kendrick$IID, diagDT$IID)]
AFR_common_Kendrick$era_numeric <- diagDT$era_numeric[match(AFR_common_Kendrick$IID, diagDT$IID)]

##############################################  SJLIFE + CCSS EUR ##################################################
dat_final <- EUR_common_Kendrick
analysis.group <- "TTN SJLIFE+CCSS EUR"
#########
## TTN ##
#########
rsID <- "rs3829746"
#------------------1. ALL



# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric +cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
  stringsAsFactors = FALSE
)
# n=c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")"))

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

ttn.sjlife.ccss.eur <- table_data




##########
## BAG3 ##
##########
rsID <- "rs2234962"
#------------------1. All

analysis.group <- "BAG3 SJLIFE+CCSS EUR"

# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric +cohort_two, family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
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

bag3.sjlife.ccss.eur <- table_data


############################################## SJLIFE EUR ##########################################################

dat_final <- EUR_common_Kendrick[EUR_common_Kendrick$cohort_two == 1,]

#########
## TTN ##
#########
rsID <- "rs3829746"

#------------------1. ALL
analysis.group <- "TTN SJLIFE EUR"


# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric , family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
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

ttn.sjlife.eur <- table_data




##########
## BAG3 ##
##########
#------------------1. All
rsID <- "rs2234962"

analysis.group <- "BAG3 SJLIFE EUR"

# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric , family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
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

bag3.sjlife.eur <- table_data


############################################## CCSS EUR ##########################################################

dat_final <- EUR_common_Kendrick[EUR_common_Kendrick$cohort_two == 2,]

#########
## TTN ##
#########
rsID <- "rs3829746"

#------------------1. ALL
analysis.group <- "TTN CCSS EUR"


# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric , family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
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

ttn.ccss.eur <- table_data




##########
## BAG3 ##
##########
rsID <- "rs2234962"

#------------------1. All

analysis.group <- "BAG3 CCSS EUR"

# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric , family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
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

bag3.ccss.eur <- table_data


############################################## SJLIFE AFR ##########################################################

dat_final <- AFR_common_Kendrick

#########
## TTN ##
#########
rsID <- "rs3829746"

#------------------1. ALL

analysis.group <- "TTN SJLIFE AFR"

# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs3829746 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric , family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
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

ttn.sjlife.afr <- table_data




##########
## BAG3 ##
##########
rsID <- "rs2234962"

#------------------1. All

analysis.group <- "BAG3 SJLIFE AFR"

# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ rs2234962 + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + era_numeric , family = binomial, data = all_data)
overall_summary <- summary(overall_model)


# Create the table
table_data <- data.frame(
  Group = c(analysis.group),
  Estimate = c(overall_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4]),
  n = c(paste0(sum(!is.na(all_data[rsID])), " (", sum(all_data$CMP==1), ")")),
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

bag3.sjlife.afr <- table_data



objects_list <- list(ttn.sjlife.eur, ttn.ccss.eur, ttn.sjlife.ccss.eur, ttn.sjlife.afr, bag3.sjlife.eur, bag3.ccss.eur, bag3.sjlife.ccss.eur, bag3.sjlife.afr)

# Extract information and combine into a data frame
combined_info <- do.call(rbind, lapply(objects_list, extract_info))

# Print the combined information
print(combined_info)

## metal file
data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/cmp3_meta_analysis_results_1.tbl", header = TRUE, sep = "\t")

# Calculate OR_CI as effect size and standard error
data$OR_CI <- paste(round(exp(data$Effect), 3), " (", round(exp(data$Effect - 1.96 * data$StdErr), 3), "-", round(exp(data$Effect + 1.96 * data$StdErr), 3), ")", sep = "")

# Select and display relevant columns
result <- data[, c("MarkerName", "P.value", "OR_CI")]
result
