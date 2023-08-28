library(dplyr)
check_grades_eq_or_higher_than <- function(data, y){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique tb_number with grades two or higher
  tb_number_with_higher_grades <- data %>%
    group_by(tb_number) %>%
    filter(grade >= y) %>%
    distinct(tb_number) %>%
    pull(tb_number)
  
  # Count the number of unique tb_number with grades two or higher
  num_tb_number_with_higher_grades <- length(tb_number_with_higher_grades)
  
  cat("Number of tb_number with grades", y, "or higher:", num_tb_number_with_higher_grades, "\n")
  # if (num_tb_number_with_higher_grades >= y) {
  cat("tb_number with grades", y, "or higher:", paste(tb_number_with_higher_grades, collapse = ", "), "\n")
  # }
  
  return(num_tb_number_with_higher_grades)
}


check_grades_transition <- function(data, lagging_grade, leading_grade){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique tb_number with grades two or higher
  tb_number_with_higher_grades <- data %>%
    group_by(tb_number) %>%
    filter(
      any(grade == lagging_grade) &  # Check if any grade is lagging_grade
        any(grade == leading_grade) # Check if any grade is leading_grade
    ) %>%
    distinct(tb_number) %>%
    pull(tb_number)
  
  # Count the number of unique tb_number with grades two or higher
  num_tb_number_with_higher_grades <- length(tb_number_with_higher_grades)
  
  cat("Number of tb_number with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", num_tb_number_with_higher_grades, "\n")
  cat("Number of tb_number with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", paste(tb_number_with_higher_grades, collapse = ", "), "\n")
  return(num_tb_number_with_higher_grades)
}



check_grades_eq_or_higher_than.min.agevent.grade.0 <- function(data, y){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique tb_number with grades two or higher
  tb_number_with_higher_grades <- data %>%
    group_by(tb_number) %>%
    filter(!is.na(ageevent) & grade >= y & grade[which.min(ageevent)] == 0) %>%
    distinct(tb_number) %>%
    pull(tb_number)
  
  # Count the number of unique tb_number with grades two or higher
  num_tb_number_with_higher_grades <- length(tb_number_with_higher_grades)
  
  cat("Number of tb_number with grades", y, "or higher:", num_tb_number_with_higher_grades, "\n")
  # if (num_tb_number_with_higher_grades >= y) {
  cat("tb_number with grades", y, "or higher:", paste(tb_number_with_higher_grades, collapse = ", "), "\n")
  # }
  
  return(num_tb_number_with_higher_grades)
}


check_grades_transition.agevent.grade.0 <- function(data, lagging_grade, leading_grade){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique tb_number with grades two or higher
  tb_number_with_higher_grades <- data %>%
    group_by(tb_number) %>%
    filter(
      any(grade == lagging_grade) &  # Check if any grade is lagging_grade
        any(grade == leading_grade) &   # Check if any grade is leading_grade
        grade[which.min(ageevent)] == lagging_grade  # Check if grade at minimum ageevent is lagging_grade
    ) %>%
    distinct(tb_number) %>%
    pull(tb_number)
  
  # Count the number of unique tb_number with grades two or higher
  num_tb_number_with_higher_grades <- length(tb_number_with_higher_grades)
  
  cat("Number of tb_number with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", num_tb_number_with_higher_grades, "\n")
  cat("Number of tb_number with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", paste(tb_number_with_higher_grades, collapse = ", "), "\n")
  return(num_tb_number_with_higher_grades)
}