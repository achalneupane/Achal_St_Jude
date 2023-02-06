## Distribution of overlapping "Unknown" lifestle factors
get_missing_combinations <- function(df){
  # Create a matrix to store the binary representation of missing values for each sample
  column_names <- colnames(df)[1:ncol(df)]
  missing_matrix <- matrix(nrow = nrow(df), ncol = ncol(df))
  colnames(missing_matrix) <- column_names
  
  # Loop through each column and create the binary representation of missing values for each sample
  for (i in 1:ncol(df)) {
    missing_matrix[, i] <- ifelse(df[, i] == "Unknown", 1, 0)
  }
  head(missing_matrix)
  # Find the unique combinations of missing values across all columns
  combinations <- apply(missing_matrix, 1, function(x) paste(column_names[x == 1], collapse = ","))
  combinations <- table(combinations)
  return(combinations)
}

# get_missing_combinations(df)