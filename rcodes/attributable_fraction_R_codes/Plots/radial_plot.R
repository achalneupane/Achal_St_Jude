library(ggplot2)
library(dplyr)

# Example data
data <- data.frame(
  Ethnicity = c("European", "African", "Asian"),
  Breakfast = c(300, 250, 200),
  Lunch = c(500, 400, 350),
  Dinner = c(600, 550, 500)
)

# Reshape the data from wide to long format
data_long <- data %>%
  tidyr::pivot_longer(cols = -Ethnicity, names_to = "Meal", values_to = "Weight")


# Plot
ggplot(data_long, aes(x = Meal, y = Weight, fill = Ethnicity, label = Weight)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(y = Weight , angle = 0), position = position_dodge(width = 1), vjust = 0) +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), color = "black", linetype = "dashed") +
  coord_polar() +
  theme_minimal() +
  labs(title = "Food Weight by Meal and Ethnicity",
       x = NULL, y = "Weight (grams)",
       fill = "Ethnicity")

