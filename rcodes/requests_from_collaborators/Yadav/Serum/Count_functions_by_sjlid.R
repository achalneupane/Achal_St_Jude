library(dplyr)
check_grades_eq_or_higher_than <- function(data, y){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique sjlid with grades two or higher
  sjlid_with_higher_grades <- data %>%
    group_by(sjlid) %>%
    filter(grade >= y) %>%
    distinct(sjlid) %>%
    pull(sjlid)
  
  # Count the number of unique sjlid with grades two or higher
  num_sjlid_with_higher_grades <- length(sjlid_with_higher_grades)
  
  cat("Number of sjlid with grades", y, "or higher:", num_sjlid_with_higher_grades, "\n")
  # if (num_sjlid_with_higher_grades >= y) {
  cat("sjlid with grades", y, "or higher:", paste(sjlid_with_higher_grades, collapse = ", "), "\n")
  # }
  
  return(num_sjlid_with_higher_grades)
}


check_grades_transition <- function(data, lagging_grade, leading_grade){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique sjlid with grades two or higher
  sjlid_with_higher_grades <- data %>%
    group_by(sjlid) %>%
    filter(
      any(grade == lagging_grade) &  # Check if any grade is lagging_grade
        any(grade == leading_grade) # Check if any grade is leading_grade
    ) %>%
    distinct(sjlid) %>%
    pull(sjlid)
  
  # Count the number of unique sjlid with grades two or higher
  num_sjlid_with_higher_grades <- length(sjlid_with_higher_grades)
  
  cat("Number of sjlid with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", num_sjlid_with_higher_grades, "\n")
  cat("Number of sjlid with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", paste(sjlid_with_higher_grades, collapse = ", "), "\n")
  return(num_sjlid_with_higher_grades)
}



check_grades_eq_or_higher_than.min.agevent.grade.0 <- function(data, y){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique sjlid with grades two or higher
  sjlid_with_higher_grades <- data %>%
    group_by(sjlid) %>%
    filter(!is.na(ageevent) & grade >= y & grade[which.min(ageevent)] == 0) %>%
    distinct(sjlid) %>%
    pull(sjlid)
  
  # Count the number of unique sjlid with grades two or higher
  num_sjlid_with_higher_grades <- length(sjlid_with_higher_grades)
  
  cat("Number of sjlid with grades", y, "or higher:", num_sjlid_with_higher_grades, "\n")
  # if (num_sjlid_with_higher_grades >= y) {
  cat("sjlid with grades", y, "or higher:", paste(sjlid_with_higher_grades, collapse = ", "), "\n")
  # }
  
  return(num_sjlid_with_higher_grades)
}


check_grades_transition.agevent.grade.0 <- function(data, lagging_grade, leading_grade){
  
  # Convert "grade" column to numeric
  data$grade <- as.numeric(data$grade)
  
  # Find unique sjlid with grades two or higher
  sjlid_with_higher_grades <- data %>%
    group_by(sjlid) %>%
    filter(
      any(grade == lagging_grade) &  # Check if any grade is lagging_grade
        any(grade == leading_grade) &   # Check if any grade is leading_grade
        grade[which.min(ageevent)] == lagging_grade  # Check if grade at minimum ageevent is lagging_grade
    ) %>%
    distinct(sjlid) %>%
    pull(sjlid)
  
  # Count the number of unique sjlid with grades two or higher
  num_sjlid_with_higher_grades <- length(sjlid_with_higher_grades)
  
  cat("Number of sjlid with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", num_sjlid_with_higher_grades, "\n")
  cat("Number of sjlid with grades transitioning from ", lagging_grade, "to leading", leading_grade, " : ", paste(sjlid_with_higher_grades, collapse = ", "), "\n")
  return(num_sjlid_with_higher_grades)
}