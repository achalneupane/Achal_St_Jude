rm(list=ls())
# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/")
load("EUR.dat.PLP.RData")
# EUR_common_Kendrick <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis//EUR_CMP_analysis_data.rds")
# EUR.dat.PLP <- EUR.dat.PLP[EUR.dat.PLP$IID %in% EUR_common_Kendrick$iid,]
# dim(EUR.dat.PLP)
# # 3086   34

## Using the same number of samples that Kendric used in common variant analysis (in this EUR_CMP_analysis_data.rds file Kendric has dropped some samples with missing anthracycline and missing heart radiation)
EUR_common_Kendrick <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", header = T, sep = "\t")
EUR.dat.PLP <- EUR.dat.PLP[EUR.dat.PLP$IID %in% EUR_common_Kendrick$IID,]
dim(EUR.dat.PLP)

# saveRDS(EUR.dat.PLP, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/analysis_from_Kendrick/Final_analysis/EUR_dat_PLP_v2.rds")

# 3062   34
genes <- colnames(EUR.dat.PLP)[grepl("carrier", colnames(EUR.dat.PLP))]
## EUR analysis

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
genes <- colnames(EUR.dat.PLP)[grepl("carrier", colnames(EUR.dat.PLP))]

for (i in 1:length(genes)){
  EUR.dat.PLP$carriers <- EUR.dat.PLP[, genes[i]]
  
  # Try to run Fisher's test, handle errors if no carriers
  tryCatch({
    gene.test <- fisher.test(table(EUR.dat.PLP$CMP, EUR.dat.PLP$carriers))
    pvalue <- gene.test$p.value
    OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  }, error = function(e) {
    pvalue <- NA
    OR.CI <- NA
  })
  
  total.carriers <- sum(EUR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.cases <- sum(EUR.dat.PLP$CMP == 2, na.rm = TRUE)
  carriers.cases <- sum(EUR.dat.PLP$CMP == 2 & EUR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.controls <- sum(EUR.dat.PLP$CMP == 1, na.rm = TRUE)
  carriers.controls <- sum(EUR.dat.PLP$CMP == 1 & EUR.dat.PLP$carriers == 1, na.rm = TRUE)
  
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


## Check prevalence of rare variants EUR
sum(EUR.dat.PLP$BAG3.PLP.carrier==1)/nrow(EUR.dat.PLP)*100
# 0.58
(sum(EUR.dat.PLP$DSP.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.03
(sum(EUR.dat.PLP$LMNA.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.91
(sum(EUR.dat.PLP$MYH7.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.75
(sum(EUR.dat.PLP$SCN5A.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0
(sum(EUR.dat.PLP$TCAP.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.23
(sum(EUR.dat.PLP$TNNC1.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.03
(sum(EUR.dat.PLP$TNNT2.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.23
(sum(EUR.dat.PLP$TTN_PSI.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.20
(sum(EUR.dat.PLP$TTN_PSI_A_Band.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 0.03
(sum(EUR.dat.PLP$ALL_PLP_With_TTN_PSI_A_Band.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 2.71
(sum(EUR.dat.PLP$ALL_PLP_With_TTN_PSI.PLP.carrier==1)/nrow(EUR.dat.PLP))*100
# 2.84

EUR.dat.PLP.CCM <- EUR.dat.PLP[EUR.dat.PLP$CMP ==2,]
EUR.dat.PLP.without.CCM <- EUR.dat.PLP[EUR.dat.PLP$CMP ==1,]

(sum(EUR.dat.PLP.CCM$TTN_PSI.PLP.carrier==1)/nrow(EUR.dat.PLP.CCM))*100
# 0
(sum(EUR.dat.PLP.without.CCM$TTN_PSI.PLP.carrier==1)/nrow(EUR.dat.PLP.without.CCM))*100
# 0.21

(sum(EUR.dat.PLP.CCM$TTN_PSI_A_Band.PLP.carrier==1)/nrow(EUR.dat.PLP.CCM))*100
# 0
(sum(EUR.dat.PLP.without.CCM$TTN_PSI_A_Band.PLP.carrier==1)/nrow(EUR.dat.PLP.without.CCM))*100
# 0.04

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


data = AFR.dat.PLP.subset
carrier_column="BAG3.PLP.carrier"
# Function to calculate the proportion and its 95% CI for each carrier gene
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
EUR.dat.PLP.subset <- EUR.dat.PLP[EUR.dat.PLP$CMP == 2, ]

# Get the carrier column names
carrier_columns <- colnames(EUR.dat.PLP.subset)[grepl("carrier", colnames(EUR.dat.PLP.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(EUR.dat.PLP.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)


# 1. EUR controls
EUR.dat.PLP.subset <- EUR.dat.PLP[EUR.dat.PLP$CMP == 1, ]

# Get the carrier column names
carrier_columns <- colnames(EUR.dat.PLP.subset)[grepl("carrier", colnames(EUR.dat.PLP.subset))]

# Initialize a list to store results for all genes
results_list <- list()

# Loop through each carrier column and calculate proportions and CIs
for (carrier in carrier_columns) {
  result <- calculate_proportion_and_ci(EUR.dat.PLP.subset, carrier)
  results_list[[carrier]] <- result
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Print the results
print(results_df)


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
