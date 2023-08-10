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


