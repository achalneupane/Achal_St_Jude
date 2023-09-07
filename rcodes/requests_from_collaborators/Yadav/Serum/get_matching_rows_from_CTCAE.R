get_matching_rows <- function(PLASMA, CTCAE, days){
PLASMA$ageevent <- NA  # Initialize ageevent with NA values
PLASMA$grade <- NA     # Initialize grade with NA values

for (i in 1:nrow(PLASMA)) {
  # Subset PLASMA and CTCAE data for the current SJLID
  PLASMA_subset <- PLASMA[i, ]
  # Find matching rows in CTCAE based on sjlid and age difference
  matching_rows <- which(
    (CTCAE$sjlid == PLASMA_subset$sjlid) &
      (abs(CTCAE$ageevent - PLASMA_subset$ageatsample) <= days)
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
