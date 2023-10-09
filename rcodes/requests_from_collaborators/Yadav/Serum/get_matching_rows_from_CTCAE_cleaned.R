## 1. Get serum age within X-days of ageevent; also extracts vials and vital status
get_rows_with_smaller_sample_age.all <- function(CTCAE, SERUM, days){
  CTCAE$Sample_age <- NA  # Initialize ageevent with NA values
  CTCAE$num_vials <- NA
  CTCAE$vitalstatus <- NA
  for (i in 1:nrow(CTCAE)) {
    # Subset SERUM and CTCAE data for the current SJLID
    CTCAE_subset <- CTCAE[i, ]
    # Find matching rows in CTCAE based on sjlid and age difference
    matching_rows <- which(
      (SERUM$sjlid == CTCAE_subset$sjlid) &
        # (abs(CTCAE_subset$ageevent - SERUM$ageatsample) <= days/365.25)
        (abs(CTCAE_subset$ageevent - SERUM$ageatsample) <= days/365.25)
    )
    # Check if there are matching rows
    if (length(matching_rows) > 0) {
      # If matching rows exist, update the ageevent and grade using the first matching row
      CTCAE$Sample_age[i] <- SERUM$ageatsample[matching_rows[1]]
      CTCAE$num_vials[i] <- SERUM$num_vials[matching_rows[1]]
      CTCAE$vitalstatus[i] <- SERUM$vitalstatus[matching_rows[1]]
    }
  }
  return(CTCAE)
}


## 2. Remove rows once grades 2 or higher are seen in ordered df ordered by event_number
filter_rows_by_condition <- function(data, group_column, condition_column) {
  # Split the dataframe by the group_column
  split_data <- split(data, data[[group_column]])
  # Define a function to filter and retain rows within each group
  filter_within_group <- function(group_df) {
    first_condition_true_row <- which(group_df[[condition_column]] >= 2)[1]
    if (is.na(first_condition_true_row)) {
      return(group_df)
    } else {
      return(group_df[1:first_condition_true_row, ])
    }
  }
  # Apply the filtering function to each group and store the results in a list
  filtered_data_list <- lapply(split_data, filter_within_group)
  # Combine the filtered data frames into a single dataframe
  result <- do.call(rbind, filtered_data_list)
  # Reset row names
  rownames(result) <- NULL
  return(result)
}


# ## 3. To calculate age difference from the last visit with grade 0 to the first CMP
# # data = CTCAE.2
# calculate_ageevent_difference <- function(data) {
#   # Sort the data by sjlid and gradedt
#   sorted_data <- data[order(data$sjlid, data$gradedt), ]
#   # Initialize variables to store results
#   result <- data.frame(sjlid = character(0), ageevent_difference = double(0))
#   # Loop through unique sjlid values
#   unique_sjlids <- unique(sorted_data$sjlid)
#   for (sjlid in unique_sjlids) {
#     # Subset data for the current sjlid
#     sjlid_data <- sorted_data[sorted_data$sjlid == sjlid, ]
#     # Find the index of the first occurrence where grade > 0
#     first_gt_0_index <- which(sjlid_data$grade > 0)[1]
#     # Find the index of the last occurrence where grade == 0 before grade > 0
#     last_eq_0_index <- max(which(sjlid_data$grade == 0 & seq_along(sjlid_data$grade) < first_gt_0_index))
#     # Calculate the ageevent difference
#     ageevent_difference <- sjlid_data$ageevent[first_gt_0_index] - sjlid_data$ageevent[last_eq_0_index]
#     # Append the result to the data frame
#     result <- rbind(result, data.frame(sjlid = sjlid, ageevent_difference = ageevent_difference))
#   }
#   # Remove rows with NA ageevent_difference
#   result <- result[!is.na(result$ageevent_difference), ]
#   return(result)
# }


## 3. To calculate age difference from the first visit to the first CMP
calculate_ageevent_difference <- function(data) {
  # Sort the data by sjlid and gradedt
  sorted_data <- data[order(data$sjlid, data$gradedt), ]
  # Initialize variables to store results
  result <- data.frame(sjlid = character(0), ageevent_difference = double(0))
  # Loop through unique sjlid values
  unique_sjlids <- unique(sorted_data$sjlid)
  for (sjlid in unique_sjlids) {
    # Subset data for the current sjlid
    sjlid_data <- sorted_data[sorted_data$sjlid == sjlid, ]
    # Find the index of the last occurrence where grade == 0 before grade > 0
    first_eq_0_index <- max(which(sjlid_data$grade == 0))
    # Find the index of the first occurrence where grade > 0
    first_gt_0_index <- which(sjlid_data$grade > 0)[1]
    # Calculate the ageevent difference
    ageevent_difference <- sjlid_data$ageevent[first_gt_0_index] - sjlid_data$ageevent[first_eq_0_index]
    # Append the result to the data frame
    result <- rbind(result, data.frame(sjlid = sjlid, ageevent_difference = ageevent_difference))
  }
  # Remove rows with NA ageevent_difference
  result <- result[!is.na(result$ageevent_difference), ]
  return(result)
}


