rm(list=ls())
# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/")
load("EUR.dat.PLP_with_ccss_org.RData")
# EUR_common_Kendrick <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis//EUR_CMP_analysis_data.rds")
# EUR.dat.PLP <- EUR.dat.PLP[EUR.dat.PLP$IID %in% EUR_common_Kendrick$iid,]
dim(EUR.dat.PLP)
# # 3086   34 (SJLIFE+CCSS exp)
# 6063 36 (SJLIFE+CCSS exp + CCSS exp)
# ## Using the same number of samples that Kendric used in common variant analysis (in this EUR_CMP_analysis_data.rds file Kendric has dropped some samples with missing anthracycline and missing heart radiation)
# EUR_common_Kendrick <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", header = T, sep = "\t")
## Using the same from KL and CCSS orginal from MIngjuan
EUR_common_Kendrick <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ccss_org_KL_Mingjuan.pheno", header = T, sep = "\t")
EUR.dat.PLP <- EUR.dat.PLP[EUR.dat.PLP$IID %in% EUR_common_Kendrick$IID,]
dim(EUR.dat.PLP)
# 6063
# saveRDS(EUR.dat.PLP, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis/EUR_dat_PLP_v2.rds")

# 3062   34
genes <- colnames(EUR.dat.PLP)[grepl("carrier", colnames(EUR.dat.PLP))]
## EUR analysis



#########################################################################
#########################################################################
#########################################################################
#########################################################################


#################
## SJLIFE only ##
#################

EUR.dat.PLP.sjlife <- EUR.dat.PLP[EUR.dat.PLP$cohort_two == 1,]

# Empty dataframe
results <- data.frame(
  Gene = character(),
  TotalCarriers = numeric(),
  TotalCases = numeric(),
  CarriersCases = numeric(),
  TotalControls = numeric(),
  CarriersControls = numeric(),
  OR_CI = character(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene set/carrier variables
genes <- colnames(EUR.dat.PLP.sjlife)[grepl("carrier", colnames(EUR.dat.PLP.sjlife))]

for (i in 1:length(genes)){
  EUR.dat.PLP.sjlife$carriers <- EUR.dat.PLP.sjlife[, genes[i]]
  
  # Try to run Fisher's test, handle errors if no carriers
  tryCatch({
    gene.test <- fisher.test(table(EUR.dat.PLP.sjlife$CMP, EUR.dat.PLP.sjlife$carriers))
    pvalue <- gene.test$p.value
    OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  }, error = function(e) {
    pvalue <- NA
    OR.CI <- NA
  })
  
  total.carriers <- sum(EUR.dat.PLP.sjlife$carriers == 1, na.rm = TRUE)
  total.cases <- sum(EUR.dat.PLP.sjlife$CMP == 2, na.rm = TRUE)
  carriers.cases <- sum(EUR.dat.PLP.sjlife$CMP == 2 & EUR.dat.PLP.sjlife$carriers == 1, na.rm = TRUE)
  total.controls <- sum(EUR.dat.PLP.sjlife$CMP == 1, na.rm = TRUE)
  carriers.controls <- sum(EUR.dat.PLP.sjlife$CMP == 1 & EUR.dat.PLP.sjlife$carriers == 1, na.rm = TRUE)
  
  # Append the results 
  results <- rbind(results, data.frame(
    Gene = genes[i],
    TotalCarriers = total.carriers,
    TotalCases = total.cases,
    CarriersCases = carriers.cases,
    TotalControls = total.controls,
    CarriersControls = carriers.controls,
    OR_CI = OR.CI,
    PValue = pvalue
  ))
}

# View the results
print(results)
View(results)

## Check prevalence of rare variants EUR
EUR.dat.PLP.sjlife.CCM <- EUR.dat.PLP.sjlife[EUR.dat.PLP.sjlife$CMP ==2,]
EUR.dat.PLP.sjlife.without.CCM <- EUR.dat.PLP.sjlife[EUR.dat.PLP.sjlife$CMP ==1,]
sum(EUR.dat.PLP.sjlife.without.CCM$LMNA.PLP.carrier==1)/nrow(EUR.dat.PLP.sjlife.without.CCM)*100


# data = EUR.dat.PLP.sjlife
# carrier_column="BAG3.PLP.carrier"
## Function to calculate the proportion and its 95% CI for each carrier gene
calculate_proportion_and_ci <- function(data, carrier_column) {
  # Number of carriers
  num_carriers <- sum(data[[carrier_column]] == 1, na.rm = TRUE)
  
  # Total sample size
  total_samples <- nrow(data)
  
  # Calculate binomial test for exact CI
  binom_result <- binom.test(num_carriers, total_samples, conf.level = 0.95)
  
  # Calculate proportion and confidence intervals
  proportion <- round((num_carriers / total_samples) * 100, 2)  # Convert to percentage and round
  ci_lower <- round(binom_result$conf.int[1] * 100, 2)          # Lower CI in percentage
  ci_upper <- round(binom_result$conf.int[2] * 100, 2)          # Upper CI in percentage
  
  # Create a list to store results
  result <- list(
    Gene = carrier_column,
    Proportion = proportion,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    WANTED_VAR = paste0(num_carriers, " (", proportion, " [", ci_lower, "-", ci_upper, "]", ")")
  )
  
  return(result)
}

# 1. EUR Cases 
EUR.dat.PLP.sjlife.subset <- EUR.dat.PLP.sjlife[EUR.dat.PLP.sjlife$CMP == 2, ]

# Get the carrier column names
carrier_columns <- colnames(EUR.dat.PLP.sjlife.subset)[grepl("carrier", colnames(EUR.dat.PLP.sjlife.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(EUR.dat.PLP.sjlife.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)


# 1. EUR controls
EUR.dat.PLP.sjlife.subset <- EUR.dat.PLP.sjlife[EUR.dat.PLP.sjlife$CMP == 1, ]

# Get the carrier column names
carrier_columns <- colnames(EUR.dat.PLP.sjlife.subset)[grepl("carrier", colnames(EUR.dat.PLP.sjlife.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(EUR.dat.PLP.sjlife.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)


###############
## CCSS only ##
###############

EUR.dat.PLP.ccss <- EUR.dat.PLP[EUR.dat.PLP$cohort_two == 2,]

# Empty dataframe
results <- data.frame(
  Gene = character(),
  TotalCarriers = numeric(),
  TotalCases = numeric(),
  CarriersCases = numeric(),
  TotalControls = numeric(),
  CarriersControls = numeric(),
  OR_CI = character(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene set/carrier variables
genes <- colnames(EUR.dat.PLP.ccss)[grepl("carrier", colnames(EUR.dat.PLP.ccss))]

for (i in 1:length(genes)){
  EUR.dat.PLP.ccss$carriers <- EUR.dat.PLP.ccss[, genes[i]]
  
  # Try to run Fisher's test, handle errors if no carriers
  tryCatch({
    gene.test <- fisher.test(table(EUR.dat.PLP.ccss$CMP, EUR.dat.PLP.ccss$carriers))
    pvalue <- gene.test$p.value
    OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  }, error = function(e) {
    pvalue <- NA
    OR.CI <- NA
  })
  
  total.carriers <- sum(EUR.dat.PLP.ccss$carriers == 1, na.rm = TRUE)
  total.cases <- sum(EUR.dat.PLP.ccss$CMP == 2, na.rm = TRUE)
  carriers.cases <- sum(EUR.dat.PLP.ccss$CMP == 2 & EUR.dat.PLP.ccss$carriers == 1, na.rm = TRUE)
  total.controls <- sum(EUR.dat.PLP.ccss$CMP == 1, na.rm = TRUE)
  carriers.controls <- sum(EUR.dat.PLP.ccss$CMP == 1 & EUR.dat.PLP.ccss$carriers == 1, na.rm = TRUE)
  
  # Append the results 
  results <- rbind(results, data.frame(
    Gene = genes[i],
    TotalCarriers = total.carriers,
    TotalCases = total.cases,
    CarriersCases = carriers.cases,
    TotalControls = total.controls,
    CarriersControls = carriers.controls,
    OR_CI = OR.CI,
    PValue = pvalue
  ))
}

# View the results
print(results)
View(results)

## Check prevalence of rare variants EUR
EUR.dat.PLP.ccss.CCM <- EUR.dat.PLP.ccss[EUR.dat.PLP.ccss$CMP ==2,]
EUR.dat.PLP.ccss.without.CCM <- EUR.dat.PLP.ccss[EUR.dat.PLP.ccss$CMP ==1,]
sum(EUR.dat.PLP.ccss.CCM$MYH7.PLP.carrier==1)/nrow(EUR.dat.PLP.ccss.CCM)*100
sum(EUR.dat.PLP.ccss.without.CCM$MYH7.PLP.carrier==1)/nrow(EUR.dat.PLP.ccss.without.CCM)*100

# data = EUR.dat.PLP.ccss
# carrier_column="BAG3.PLP.carrier"
## Function to calculate the proportion and its 95% CI for each carrier gene
calculate_proportion_and_ci <- function(data, carrier_column) {
  # Number of carriers
  num_carriers <- sum(data[[carrier_column]] == 1, na.rm = TRUE)
  
  # Total sample size
  total_samples <- nrow(data)
  
  # Calculate binomial test for exact CI
  binom_result <- binom.test(num_carriers, total_samples, conf.level = 0.95)
  
  # Calculate proportion and confidence intervals
  proportion <- round((num_carriers / total_samples) * 100, 2)  # Convert to percentage and round
  ci_lower <- round(binom_result$conf.int[1] * 100, 2)          # Lower CI in percentage
  ci_upper <- round(binom_result$conf.int[2] * 100, 2)          # Upper CI in percentage
  
  # Create a list to store results
  result <- list(
    Gene = carrier_column,
    Proportion = proportion,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    WANTED_VAR = paste0(num_carriers, " (", proportion, " [", ci_lower, "-", ci_upper, "]", ")")
  )
  
  return(result)
}

# 1. EUR Cases 
EUR.dat.PLP.ccss.subset <- EUR.dat.PLP.ccss[EUR.dat.PLP.ccss$CMP == 2, ]

# Get the carrier column names
carrier_columns <- colnames(EUR.dat.PLP.ccss.subset)[grepl("carrier", colnames(EUR.dat.PLP.ccss.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(EUR.dat.PLP.ccss.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)


# 1. EUR controls
EUR.dat.PLP.ccss.subset <- EUR.dat.PLP.ccss[EUR.dat.PLP.ccss$CMP == 1, ]

# Get the carrier column names
carrier_columns <- colnames(EUR.dat.PLP.ccss.subset)[grepl("carrier", colnames(EUR.dat.PLP.ccss.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(EUR.dat.PLP.ccss.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)

##################
## AFR Analysis ##
##################

rm(list=ls())
## AFR analysis
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/AFR.dat.PLP.RData")
# AFR_common_Kendrick <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis/AFR_CMP_analysis_data.rds")
# AFR.dat.PLP <- AFR.dat.PLP[AFR.dat.PLP$IID %in% AFR_common_Kendrick$iid,]
# dim(AFR.dat.PLP)
# # 244   34

## Using the same number of samples that Kendric used in common variant analysis (in this EUR_CMP_analysis_data.rds file Kendric has dropped some samples with missing anthracycline and missing heart radiation)
AFR_common_Kendrick <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3_afr_kendrick.pheno", header = T, sep = "\t")
AFR.dat.PLP <- AFR.dat.PLP[AFR.dat.PLP$IID %in% AFR_common_Kendrick$IID,]
dim(AFR.dat.PLP)
# 238  34

genes <- colnames(AFR.dat.PLP)[grepl("carrier", colnames(AFR.dat.PLP))]
# saveRDS(AFR.dat.PLP, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis/AFR_dat_PLP_v2.rds")
## AFR analysis

# Empty dataframe
results <- data.frame(
  Gene = character(),
  TotalCarriers = numeric(),
  TotalCases = numeric(),
  CarriersCases = numeric(),
  TotalControls = numeric(),
  CarriersControls = numeric(),
  OR_CI = character(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene set/carrier variables
genes <- colnames(AFR.dat.PLP)[grepl("carrier", colnames(AFR.dat.PLP))]

for (i in 1:length(genes)){
  AFR.dat.PLP$carriers <- AFR.dat.PLP[, genes[i]]
  
  # Try to run Fisher's test, handle errors if no carriers
  tryCatch({
    gene.test <- fisher.test(table(AFR.dat.PLP$CMP, AFR.dat.PLP$carriers))
    pvalue <- gene.test$p.value
    OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  }, error = function(e) {
    pvalue <- NA
    OR.CI <- NA
  })
  
  total.carriers <- sum(AFR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.cases <- sum(AFR.dat.PLP$CMP == 2, na.rm = TRUE)
  carriers.cases <- sum(AFR.dat.PLP$CMP == 2 & AFR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.controls <- sum(AFR.dat.PLP$CMP == 1, na.rm = TRUE)
  carriers.controls <- sum(AFR.dat.PLP$CMP == 1 & AFR.dat.PLP$carriers == 1, na.rm = TRUE)
  
  # Append the results 
  results <- rbind(results, data.frame(
    Gene = genes[i],
    TotalCarriers = total.carriers,
    TotalCases = total.cases,
    CarriersCases = carriers.cases,
    TotalControls = total.controls,
    CarriersControls = carriers.controls,
    OR_CI = OR.CI,
    PValue = pvalue
  ))
}

# View the results
print(results)

## Check prevalence of rare variants in AFR
sum(AFR.dat.PLP$BAG3.PLP.carrier==1)/nrow(AFR.dat.PLP)*100
# 0.42
(sum(AFR.dat.PLP$DSP.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 0
(sum(AFR.dat.PLP$LMNA.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 0.84
(sum(AFR.dat.PLP$MYH7.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 2.1
(sum(AFR.dat.PLP$SCN5A.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 0
(sum(AFR.dat.PLP$TCAP.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 0
(sum(AFR.dat.PLP$TNNT2.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 0.84
(sum(AFR.dat.PLP$TTN_PSI.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 0.84
(sum(AFR.dat.PLP$TTN_PSI_A_Band.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 0
(sum(AFR.dat.PLP$ALL_PLP_With_TTN_PSI_A_Band.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 4.2
(sum(AFR.dat.PLP$ALL_PLP_With_TTN_PSI.PLP.carrier==1)/nrow(AFR.dat.PLP))*100
# 5.04


# 1. AFR cases
AFR.dat.PLP.subset <- AFR.dat.PLP[AFR.dat.PLP$CMP == 2, ]

# Get the carrier column names
carrier_columns <- colnames(AFR.dat.PLP.subset)[grepl("carrier", colnames(AFR.dat.PLP.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(AFR.dat.PLP.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)


# 1. AFR controls
AFR.dat.PLP.subset <- AFR.dat.PLP[AFR.dat.PLP$CMP == 1, ]

# Get the carrier column names
carrier_columns <- colnames(AFR.dat.PLP.subset)[grepl("carrier", colnames(AFR.dat.PLP.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(AFR.dat.PLP.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)


#######################################
## Meta analysis with Cochran–Mantel ##
#######################################

## New analysis with 2 X 2 matrix
EUR.dat.PLP <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis/EUR_dat_PLP_v2.rds")

library(dplyr)

data <- EUR.dat.PLP

# List of carrier status columns
carrier_cols <- grep("\\.PLP\\.carrier$", colnames(data), value = TRUE)

# Function to calculate Mantel-Haenszel OR, p-value, and detailed carrier counts
analyze_carrier <- function(carrier_col) {
  # Creating a contingency table stratified by cohort_two
  table_list <- lapply(split(data, data$cohort_two), function(subdata) {
    table(subdata$CMP, subdata[[carrier_col]])
  })
  
  # combining tables into a 3D array
  data_array <- array(unlist(table_list), dim = c(2, 2, length(table_list)))
  
  # Mantel-Haenszel test
  test_result <- mantelhaen.test(data_array)
  
  # extracting OR, CI, and p-value
  OR <- test_result$estimate
  CI <- test_result$conf.int
  p_value <- test_result$p.value
  
  # calculating carrier counts
  total_carriers <- sum(data[[carrier_col]] == 1)
  cohort1_carriers <- sum(data[[carrier_col]] == 1 & data$cohort_two == 1)
  cohort2_carriers <- sum(data[[carrier_col]] == 1 & data$cohort_two == 2)
  cmp_cohort1_carriers <- sum(data[[carrier_col]] == 1 & data$CMP == 1 & data$cohort_two == 1)
  no_cmp_cohort1_carriers <- sum(data[[carrier_col]] == 1 & data$CMP == 2 & data$cohort_two == 1)
  cmp_cohort2_carriers <- sum(data[[carrier_col]] == 1 & data$CMP == 1 & data$cohort_two == 2)
  no_cmp_cohort2_carriers <- sum(data[[carrier_col]] == 1 & data$CMP == 2 & data$cohort_two == 2)
  
  # return results
  return(data.frame(
    Carrier = carrier_col,
    Total_Carriers = total_carriers,
    Cohort1_Carriers = cohort1_carriers,
    Cohort2_Carriers = cohort2_carriers,
    CMP_Cohort1_Carriers = cmp_cohort1_carriers,
    No_CMP_Cohort1_Carriers = no_cmp_cohort1_carriers,
    CMP_Cohort2_Carriers = cmp_cohort2_carriers,
    No_CMP_Cohort2_Carriers = no_cmp_cohort2_carriers,
    OR = paste0(round(OR, 2), paste0(" (", round(CI[1], 2), ", ", round(CI[2], 2), ")")), # OR and CI
    # CI = paste0("(", round(CI[1], 2), ", ", round(CI[2], 2), ")"),
    P_value = signif(p_value, 3)
  ))
}

# running analysis for all carrier statuses
results <- bind_rows(lapply(carrier_cols, analyze_carrier))

print(results)
