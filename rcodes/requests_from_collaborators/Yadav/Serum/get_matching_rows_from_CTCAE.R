get_matching_rows <- function(PLASMA, CTCAE, days){
PLASMA$ageevent <- NA  # Initialize ageevent with NA values
PLASMA$grade <- NA     # Initialize grade with NA values

for (i in 1:nrow(PLASMA)) {
  # Subset PLASMA and CTCAE data for the current SJLID
  PLASMA_subset <- PLASMA[i, ]
  # Find matching rows in CTCAE based on sjlid and age difference
  matching_rows <- which(
    (CTCAE$sjlid == PLASMA_subset$sjlid) &
      (abs(CTCAE$ageevent - PLASMA_subset$ageatsample) <= days/365.25)
  )
  # Check if there are matching rows
  if (length(matching_rows) > 0) {
    # If matching rows exist, update the ageevent and grade using the first matching row
    PLASMA$ageevent[i] <- CTCAE$ageevent[matching_rows[1]]
    PLASMA$grade[i] <- CTCAE$grade[matching_rows[1]]
  }
}

return(PLASMA)
}


# cc <- get_matching_rows(PLASMA.cc, CTCAE.cc, 7/365.25)


# CTCAE <- first_CMP_event
# SERUM <- SERUM.original
get_rows_with_smaller_sample_age <- function(CTCAE, SERUM, days){
  CTCAE$Sample_age <- NA  # Initialize ageevent with NA values
  
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
    }
  }
  return(CTCAE)
}

## also include vials
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




## Remove rows once grades 2 or higher are seen in ordered df ordered by event_number
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


# filter_rows_by_condition(small.CTCAE, "sjlid", "grade")
