count_zeros_ones <- function(data) {
  # Create an empty data frame to store the results
  result_df <- data.frame(Column = character(),
                          Zeros = integer(),
                          Ones = integer(),
                          stringsAsFactors = FALSE)
  
  # Loop through each column of the data
  for (i in 1:ncol(data)) {
    # Get the column name
    col_name <- colnames(data)[i]
    
    # Count the number of 0s and 1s in the column
    zeros_count <- sum(data[, i] == 0, na.rm = TRUE)
    ones_count <- sum(data[, i] == 1, na.rm = TRUE)
    
    # Add the results to the data frame
    result_df <- rbind(result_df, data.frame(Column = col_name,
                                             Zeros = zeros_count,
                                             Ones = ones_count))
  }
  
  return(result_df)
}

