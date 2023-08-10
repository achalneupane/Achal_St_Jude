## Distribution of overlapping "Unknown" lifestle factors
# df <- PHENO.ANY_SN
# Columns to check
# columns_to_check <- c("PhysicalActivity_yn", "Current_smoker_yn", "RiskyHeavyDrink_yn", "Obese_yn")

get_missing_combinations <- function(df, cols){
  # df <- df[grepl("ccssid|sjlid|ANY_SNs|SMNs|NMSCs|BREASTcancer|THYROIDcancer|MENINGIOMA|SARCOMA|^Current_smoker_yn$|^PhysicalActivity_yn$|^RiskyHeavyDrink_yn$|^Obese_yn$", colnames(df))]
  colnames(df)[grepl("ANY_SNs|SMNs|NMSCs|BREASTcancer|THYROIDcancer|MENINGIOMA|SARCOMA", colnames(df))] <- "CACO"
  # Create a matrix to store the binary representation of missing values for each sample
  # Create subsets for cases and controls
  cases <- df[df$CACO == "1", ]
  controls <- df[df$CACO == "0", ]
  
  # Columns to check
  columns_to_check <- cols
  
  # Function to calculate unknown percentages
  calculate_unknown_percentages <- function(data) {
    unknown_percentages <- sapply(1:4, function(k) {
      any_unknown <- rowSums(data[columns_to_check] == "Unknown") >= k
      percentage <- sum(any_unknown) / nrow(data) * 100
      return(round(percentage, 2))
    })
    return(unknown_percentages)
  }
  
  # Calculate unknown percentages for cases and controls
  unknown_percentages_cases <- calculate_unknown_percentages(cases)
  unknown_percentages_controls <- calculate_unknown_percentages(controls)
  
  combinations <- c(unknown_percentages_cases, unknown_percentages_controls)
  names(combinations) <- c("atleast_1_missing_CA", "atleast_2_missing_CA", "atleast_3_missing_CA", "all_missing_CA", "atleast_1_missing_CO", "atleast_2_missing_CO", "atleast_3_missing_CO", "all_missing_CO")
  return(combinations)
}

# get_missing_combinations(df)


## Yutaka on 08/10/2023: Could you breakdown the "any 1 missing" to each item missing so that I can see what variables are missing more
library(dplyr)

calculate_missing_counts <- function(data) {
  colnames(data)[grepl("ANY_SNs|SMNs|NMSCs|BREASTcancer|THYROIDcancer|MENINGIOMA|SARCOMA", colnames(data))] <- "CACO"
  result <- data %>%
    group_by(CACO) %>%
    summarize(
      Total_CA = sum(CACO == "1"),
      Total_CO = sum(CACO == "0"),
      
      Current_smoker_missing_CA = sum(Current_smoker_yn == "Unknown" & CACO == "1"),
      RiskyHeavyDrink_missing_CA = sum(RiskyHeavyDrink_yn == "Unknown" & CACO == "1"),
      PhysicalActivity_missing_CA = sum(PhysicalActivity_yn == "Unknown" & CACO == "1"),
      Obese_missing_CA = sum(Obese_yn == "Unknown" & CACO == "1"),
      
      Current_smoker_missing_CO = sum(Current_smoker_yn == "Unknown" & CACO == "0"),
      RiskyHeavyDrink_missing_CO = sum(RiskyHeavyDrink_yn == "Unknown" & CACO == "0"),
      PhysicalActivity_missing_CO = sum(PhysicalActivity_yn == "Unknown" & CACO == "0"),
      Obese_missing_CO = sum(Obese_yn == "Unknown" & CACO == "0")
    )
  
  return(result)
}

# # Call the function with your dataframe
# missing_counts <- calculate_missing_counts(df)
# 
# # Print the result
# print(missing_counts)



library(dplyr)

calculate_missing_percentages <- function(data) {
  colnames(data)[grepl("ANY_SNs|SMNs|NMSCs|BREASTcancer|THYROIDcancer|MENINGIOMA|SARCOMA", colnames(data))] <- "CACO"
  result <- data %>%
    group_by(CACO) %>%
    summarize(
      Total_CA = sum(CACO == "1"),
      Total_CO = sum(CACO == "0"),
      PhysicalActivity_missing_CA = (sum(PhysicalActivity_yn == "Unknown" & CACO == "1") / Total_CA) * 100,
      PhysicalActivity_missing_CO = (sum(PhysicalActivity_yn == "Unknown" & CACO == "0") / Total_CO) * 100,
      Current_smoker_missing_CA = (sum(Current_smoker_yn == "Unknown" & CACO == "1") / Total_CA) * 100,
      Current_smoker_missing_CO = (sum(Current_smoker_yn == "Unknown" & CACO == "0") / Total_CO) * 100,
      RiskyHeavyDrink_missing_CA = (sum(RiskyHeavyDrink_yn == "Unknown" & CACO == "1") / Total_CA) * 100,
      RiskyHeavyDrink_missing_CO = (sum(RiskyHeavyDrink_yn == "Unknown" & CACO == "0") / Total_CO) * 100,
      Obese_missing_CA = (sum(Obese_yn == "Unknown" & CACO == "1") / Total_CA) * 100,
      Obese_missing_CO = (sum(Obese_yn == "Unknown" & CACO == "0") / Total_CO) * 100
    )
  
  return(result)
}

# # Call the function with your dataframe
# missing_percentages <- calculate_missing_percentages(df)
# 
# # Print the result
# print(missing_percentages)

