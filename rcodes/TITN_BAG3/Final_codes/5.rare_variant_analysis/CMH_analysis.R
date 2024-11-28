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



